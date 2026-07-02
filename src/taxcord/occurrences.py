"""Count regional occurrences for each taxon in GBIF and BOLD.

For every query lineage this annotates the finest resolved rank with the
number of records GBIF and BOLD hold for the target region (by default the
Iberian Peninsula: Portugal and Spain). GBIF is queried through its public
occurrence API; BOLD is read from its public taxon page. When GBIF returns no
records at a rank the search moves one rank up.

This step needs network access and depends on the live structure of the BOLD
taxon page, so its output can change as those services change.
"""

from __future__ import annotations

import ast
import csv
import re
import sys
import time
from dataclasses import dataclass, field

import requests

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
GBIF_OCCURRENCE_URL = "https://api.gbif.org/v1/occurrence/search"
BOLD_TAXON_URL = "http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage"

DEFAULT_COUNTRIES = ("PT", "ES")
BOLD_COUNTRY_NAMES = ("Portugal", "Spain")

# Ranks reported in the output, coarsest first. GBIF can match only these ranks.
REPORTED_RANKS = ["order", "family", "genus", "species"]
GBIF_KEY_FIELD = {
    "species": ("name", "speciesKey"),
    "genus": ("genus", "genusKey"),
    "family": ("family", "familyKey"),
    "order": ("order", "orderKey"),
}

# Values that mark a rank as unresolved and not worth querying.
SKIP_TOKENS = {"uncultured", "fungus", "no-match", "IncompleteTaxonomy"}

REQUEST_TIMEOUT = 30
HEADERS = {"User-Agent": "taxcord/1.0 (taxonomy occurrence lookup)"}

# Transient HTTP statuses worth retrying (rate limiting and server errors).
RETRYABLE_STATUS = {429, 500, 502, 503, 504}
MAX_RETRIES = 2  # extra attempts after the first, for transient failures
RETRY_BACKOFF = 1.0  # seconds before the first retry, doubled each time
# Consecutive failures before a service is declared down and skipped for the
# rest of the run, so a collapsed endpoint is not hammered 1700 times over.
COLLAPSE_THRESHOLD = 8


class ServiceError(RuntimeError):
    """A network/service failure annotated with which service failed."""

    def __init__(self, service, detail):
        super().__init__(f"{service} {detail}")
        self.service = service
        self.detail = detail


@dataclass
class _Health:
    """Running failure state for one external service."""

    consecutive: int = 0
    failures: int = 0
    down: bool = False


@dataclass
class _Caches:
    """Per-run memoisation (one query per unique taxon) and service health."""

    gbif: dict = field(default_factory=dict)
    keys: dict = field(default_factory=dict)
    bold: dict = field(default_factory=dict)
    health: dict = field(default_factory=lambda: {"GBIF": _Health(), "BOLD": _Health()})


def _fetch(session, service, url, params, caches):
    """GET ``url`` with retry/backoff, returning the response.

    Transient failures (timeouts, dropped connections, 429/5xx) are retried a
    few times. After :data:`COLLAPSE_THRESHOLD` consecutive failures the
    service is marked down and further calls short-circuit immediately. Raises
    :class:`ServiceError` (naming the service and reason) on any failure.
    """
    health = caches.health[service]
    if health.down:
        raise ServiceError(service, "skipped (service marked down this run)")

    detail = "unknown error"
    for attempt in range(MAX_RETRIES + 1):
        try:
            response = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
        except requests.RequestException as exc:
            detail = f"unreachable ({exc.__class__.__name__})"
        else:
            if response.status_code in RETRYABLE_STATUS:
                detail = f"HTTP {response.status_code} {response.reason}"
            elif not response.ok:
                # Service answered but rejected this request (e.g. 400/404):
                # it is up, so do not count this toward a collapse.
                health.consecutive = 0
                raise ServiceError(
                    service, f"HTTP {response.status_code} {response.reason}"
                )
            else:
                health.consecutive = 0
                return response
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_BACKOFF * (2**attempt))

    health.consecutive += 1
    health.failures += 1
    if health.consecutive >= COLLAPSE_THRESHOLD:
        health.down = True
    raise ServiceError(service, detail)


def gbif_taxon_key(session, rank, name, caches):
    """Resolve a taxon name to its GBIF key for the given rank, or None."""
    if rank not in GBIF_KEY_FIELD:
        return None
    param, key_field = GBIF_KEY_FIELD[rank]
    cache_key = (param, name)
    if cache_key not in caches.keys:
        response = _fetch(session, "GBIF", GBIF_MATCH_URL, {param: name}, caches)
        try:
            data = response.json()
        except ValueError:
            raise ServiceError("GBIF", "invalid JSON in match response")
        caches.keys[cache_key] = data.get(key_field)
    return caches.keys[cache_key]


def gbif_occurrences(session, taxon_key, countries, caches):
    """Return the GBIF occurrence count for a taxon key in the given countries."""
    if taxon_key not in caches.gbif:
        params = [("taxonKey", taxon_key), *(("country", c) for c in countries)]
        response = _fetch(session, "GBIF", GBIF_OCCURRENCE_URL, params, caches)
        try:
            caches.gbif[taxon_key] = response.json().get("count", 0)
        except ValueError:
            raise ServiceError("GBIF", "invalid JSON in occurrence response")
    return caches.gbif[taxon_key]


def bold_occurrences(session, name, country_names, caches):
    """Return the BOLD record count for a taxon across the given countries.

    BOLD embeds per-country counts in a JavaScript ``allCountriesData`` object
    on the taxon page; this finds that object and sums the requested countries.
    """
    if name in caches.bold:
        return caches.bold[name]

    count = 0
    response = _fetch(session, "BOLD", BOLD_TAXON_URL, {"taxon": name}, caches)
    match = re.search(r"allCountriesData\s*=\s*(\{.*?\})", response.text, re.DOTALL)
    if match:
        try:
            per_country = ast.literal_eval(match.group(1))
            count = sum(per_country.get(c, 0) for c in country_names)
        except (ValueError, SyntaxError):
            count = 0

    caches.bold[name] = count
    return count


def resolved_rank_index(lineage_values, ranks):
    """Return the index of the finest resolved, queryable rank, or None."""
    for index in range(len(ranks) - 1, -1, -1):
        value = lineage_values[index]
        if not value or value == "NA":
            continue
        if value.split()[0] in SKIP_TOKENS:
            return None
        return index
    return None


def annotate_row(session, values, ranks, countries, country_names, caches):
    """Append GBIF and BOLD counts for the finest rank with GBIF support.

    Walks up from the resolved rank until GBIF reports records, then reports
    that count under its rank and the BOLD count for the original taxon. The
    two services are queried independently: if one fails (or has collapsed)
    the other's counts are still kept, and the failure is returned alongside
    the cells rather than discarding the whole row.

    Returns ``(cells, errors)`` — the appended count cells (GBIF per reported
    rank, then BOLD) and a list of :class:`ServiceError`/request failures hit
    while building them. Returns ``(None, [])`` for an unresolvable lineage.
    """
    index = resolved_rank_index(values, ranks)
    if index is None:
        return None, []

    errors = []
    bold_counts = {}
    gbif_rank = None
    gbif_count = 0
    gbif_failed = bold_failed = False
    while index >= 0:
        rank, name = ranks[index], values[index]
        if rank in REPORTED_RANKS and rank not in bold_counts and not bold_failed:
            try:
                bold_counts[rank] = bold_occurrences(
                    session, name, country_names, caches
                )
            except (ServiceError, requests.RequestException, ValueError) as exc:
                errors.append(exc)
                bold_failed = True

        if gbif_rank is None and not gbif_failed:
            try:
                taxon_key = gbif_taxon_key(session, rank, name, caches)
                if taxon_key is not None:
                    count = gbif_occurrences(session, taxon_key, countries, caches)
                    if count > 0:
                        gbif_rank, gbif_count = rank, count
            except (ServiceError, requests.RequestException, ValueError) as exc:
                errors.append(exc)
                gbif_failed = True

        if gbif_rank is not None:
            break
        index -= 1

    gbif_cells = [
        str(gbif_count) if rank == gbif_rank else "-" for rank in REPORTED_RANKS
    ]
    bold_cells = [
        str(bold_counts[rank]) if bold_counts.get(rank) else "-"
        for rank in REPORTED_RANKS
    ]
    return gbif_cells + bold_cells, errors


def _count_rows(in_path):
    """Return the number of data rows (lines minus the header)."""
    with open(in_path, encoding="utf-8") as handle:
        total = sum(1 for _ in handle)
    return max(total - 1, 0)


def _row_label(row, ranks):
    """A short human label for a row: its id and finest resolved rank value."""
    name = next(
        (v for v in reversed(row[1:]) if v and v != "NA"),
        "",
    )
    return f"{row[0]} {name}".strip()


def _report_progress(stream, done, total, errors, label):
    """Render a single-line progress bar to ``stream``.

    On a TTY this overwrites one line with a carriage return; in a log file
    (e.g. SLURM stderr) it prints discrete lines so the log stays readable.
    """
    width = len(str(total))
    pct = done * 100 // total if total else 100
    bar = "#" * (pct // 5) + "." * (20 - pct // 5)
    msg = f"[{bar}] {pct:3d}%  {done:>{width}}/{total}  {label[:45]:<45}"
    if errors:
        msg += f"  errors: {errors}"
    if stream.isatty():
        stream.write("\r" + msg)
    else:
        stream.write(msg + "\n")
    stream.flush()


def annotate_file(in_path, out_path, countries, country_names, progress=sys.stderr):
    """Annotate every lineage in a tab-delimited file with occurrence counts.

    Progress and per-row errors are written to ``progress`` (stderr by
    default); a failed lookup is reported and that row is written with empty
    counts rather than aborting the whole run. Pass ``progress=None`` to
    silence it.
    """
    caches = _Caches()
    extra_header = [f"GBIF.{rank}" for rank in REPORTED_RANKS] + [
        f"BOLD.{rank}" for rank in REPORTED_RANKS
    ]
    placeholder = ["-"] * len(extra_header)
    total = _count_rows(in_path)
    is_tty = bool(progress) and progress.isatty()

    with (
        requests.Session() as session,
        open(in_path, newline="", encoding="utf-8") as source,
        open(out_path, "w", newline="", encoding="utf-8") as out,
    ):
        session.headers.update(HEADERS)
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(out, delimiter="\t")

        header = next(reader)
        ranks = [name.lower() for name in header[1:]]
        writer.writerow(header + extra_header)

        done = errors = 0
        warned = set()
        for row in reader:
            label = _row_label(row, ranks)
            counts, row_errors = annotate_row(
                session, row[1:], ranks, countries, country_names, caches
            )
            if row_errors:
                errors += 1
                if progress:
                    end = "\n" if is_tty else ""
                    for exc in row_errors:
                        progress.write(f"{end}  ! {label}: {exc}\n")
                        end = ""
                    for service, health in caches.health.items():
                        if health.down and service not in warned:
                            warned.add(service)
                            progress.write(
                                f"  !! {service} appears to be down after "
                                f"{COLLAPSE_THRESHOLD} consecutive failures; "
                                f"skipping {service} for the rest of this run. "
                                f"Re-run later to fill the gaps.\n"
                            )
            writer.writerow(row + (counts if counts else placeholder))
            out.flush()  # write each row through so the file grows live and
            # partial results survive an interrupted run rather than sitting
            # in the buffer until close.

            done += 1
            if progress and (is_tty or done % 25 == 0 or done == total):
                _report_progress(progress, done, total, errors, label)

        if progress and is_tty:
            progress.write("\n")
        return done, errors, caches.health


def _parse_countries(specs):
    """Turn ``CODE:NAME`` strings into aligned GBIF-code and BOLD-name tuples.

    GBIF filters on the ISO 2-letter code while BOLD keys on the full country
    name, so each region is given as one ``CODE:NAME`` token (e.g.
    ``PT:Portugal``) to keep the two lists in step.
    """
    countries, names = [], []
    for spec in specs:
        code, sep, name = spec.partition(":")
        code, name = code.strip(), name.strip()
        if not sep or not code or not name:
            raise SystemExit(
                f"taxcord occurrences: --country expects CODE:NAME, got {spec!r} "
                f"(e.g. PT:Portugal)"
            )
        countries.append(code.upper())
        names.append(name)
    return tuple(countries), tuple(names)


def configure(parser):
    parser.add_argument("input", help="tab-delimited lineage table (id + ranks)")
    parser.add_argument("output", help="annotated tab-delimited output path")
    default = " ".join(
        f"{c}:{n}" for c, n in zip(DEFAULT_COUNTRIES, BOLD_COUNTRY_NAMES)
    )
    parser.add_argument(
        "--country",
        dest="countries",
        action="append",
        metavar="CODE:NAME",
        help=(
            "region to count records for, as a GBIF 2-letter code and BOLD "
            "country name, e.g. PT:Portugal. Repeat for several regions. "
            f"Defaults to: {default}"
        ),
    )


def execute(args):
    if args.countries:
        countries, country_names = _parse_countries(args.countries)
    else:
        countries, country_names = DEFAULT_COUNTRIES, BOLD_COUNTRY_NAMES
    done, errors, health = annotate_file(
        args.input, args.output, countries, country_names
    )
    summary = f"Annotated {done} lineages -> {args.output}"
    if errors:
        summary += f" ({errors} row(s) failed and were left blank)"
    breakdown = [
        f"{service}: {h.failures} failed{' — DOWN' if h.down else ''}"
        for service, h in health.items()
        if h.failures
    ]
    if breakdown:
        summary += "  [" + "; ".join(breakdown) + "]"
    print(summary, file=sys.stderr)
