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
from dataclasses import dataclass, field

import requests

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
GBIF_OCCURRENCE_URL = "https://api.gbif.org/v1/occurrence/search"
BOLD_TAXON_URL = "http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage"

DEFAULT_COUNTRIES = ("PT", "ES")
BOLD_COUNTRY_NAMES = ("Portugal", "Spain")

# Ranks reported in the output, finest first. GBIF can match only these ranks.
REPORTED_RANKS = ["species", "genus", "family", "order"]
GBIF_KEY_FIELD = {
    "species": ("name", "speciesKey"),
    "genus": ("genus", "genusKey"),
    "family": ("family", "familyKey"),
    "order": ("order", "orderKey"),
}

# Values that mark a rank as unresolved and not worth querying.
SKIP_TOKENS = {"uncultured", "fungus", "no-match", "IncompleteTaxonomy"}

REQUEST_TIMEOUT = 30
HEADERS = {"User-Agent": "metatax/1.0 (taxonomy occurrence lookup)"}


@dataclass
class _Caches:
    """Per-run memoisation so each unique taxon is queried only once."""

    gbif: dict = field(default_factory=dict)
    keys: dict = field(default_factory=dict)
    bold: dict = field(default_factory=dict)


def gbif_taxon_key(session, rank, name, cache):
    """Resolve a taxon name to its GBIF key for the given rank, or None."""
    if rank not in GBIF_KEY_FIELD:
        return None
    param, key_field = GBIF_KEY_FIELD[rank]
    cache_key = (param, name)
    if cache_key not in cache:
        response = session.get(
            GBIF_MATCH_URL, params={param: name}, timeout=REQUEST_TIMEOUT
        )
        data = response.json()
        cache[cache_key] = data.get(key_field)
    return cache[cache_key]


def gbif_occurrences(session, taxon_key, countries, cache):
    """Return the GBIF occurrence count for a taxon key in the given countries."""
    if taxon_key not in cache:
        params = [("taxonKey", taxon_key), *(("country", c) for c in countries)]
        response = session.get(
            GBIF_OCCURRENCE_URL, params=params, timeout=REQUEST_TIMEOUT
        )
        cache[taxon_key] = response.json().get("count", 0)
    return cache[taxon_key]


def bold_occurrences(session, name, country_names, cache):
    """Return the BOLD record count for a taxon across the given countries.

    BOLD embeds per-country counts in a JavaScript ``allCountriesData`` object
    on the taxon page; this finds that object and sums the requested countries.
    """
    if name in cache:
        return cache[name]

    count = 0
    response = session.get(
        BOLD_TAXON_URL, params={"taxon": name}, timeout=REQUEST_TIMEOUT
    )
    match = re.search(r"allCountriesData\s*=\s*(\{.*?\})", response.text, re.DOTALL)
    if match:
        try:
            per_country = ast.literal_eval(match.group(1))
            count = sum(per_country.get(c, 0) for c in country_names)
        except (ValueError, SyntaxError):
            count = 0

    cache[name] = count
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
    that count under its rank and the BOLD count for the original taxon.
    Returns the appended count cells (GBIF for each reported rank, then BOLD).
    """
    index = resolved_rank_index(values, ranks)
    if index is None:
        return None

    bold_counts = {}
    gbif_rank = None
    gbif_count = 0
    while index >= 0:
        rank, name = ranks[index], values[index]
        if rank in REPORTED_RANKS:
            bold_counts[rank] = bold_occurrences(
                session, name, country_names, caches.bold
            )

        taxon_key = gbif_taxon_key(session, rank, name, caches.keys)
        if taxon_key is not None:
            count = gbif_occurrences(session, taxon_key, countries, caches.gbif)
            if count > 0:
                gbif_rank, gbif_count = rank, count
                break
        index -= 1

    gbif_cells = [
        str(gbif_count) if rank == gbif_rank else "-" for rank in REPORTED_RANKS
    ]
    bold_cells = [
        str(bold_counts[rank]) if bold_counts.get(rank) else "-"
        for rank in REPORTED_RANKS
    ]
    return gbif_cells + bold_cells


def annotate_file(in_path, out_path, countries, country_names):
    """Annotate every lineage in a tab-delimited file with occurrence counts."""
    caches = _Caches()
    extra_header = [f"IP.{rank}" for rank in REPORTED_RANKS] + [
        f"BOLD.{rank}" for rank in REPORTED_RANKS
    ]

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

        for row in reader:
            counts = annotate_row(
                session, row[1:], ranks, countries, country_names, caches
            )
            writer.writerow(row + (counts if counts else ["-"] * len(extra_header)))


def configure(parser):
    parser.add_argument("input", help="tab-delimited lineage table (id + ranks)")
    parser.add_argument("output", help="annotated tab-delimited output path")


def execute(args):
    annotate_file(args.input, args.output, DEFAULT_COUNTRIES, BOLD_COUNTRY_NAMES)
    print(f"Annotated occurrences -> {args.output}", file=sys.stderr)
