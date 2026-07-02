"""Merge NCBI and BOLD taxonomy tables into one consensus lineage per query.

The two sources are joined on query id. At each rank the call is the value
they agree on; a disagreement clears that rank. Once a rank is empty the lower
ranks are cleared too, so a lineage never claims more detail than it supports.

That rule alone discards good calls: when the sources agree on the genus or
species but place it in differently-named families (synonyms or different
classifications), the family conflict cascades down and erases the agreed
species. With ``--gbif-backbone`` such a row is reconciled instead — the agreed
taxon is resolved against GBIF's name backbone, the agreed ranks are kept, and
only the unresolved coarser rank(s) are filled from GBIF's accepted
classification, so the species is kept.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import requests

LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Columns of the reconciliation report written alongside the merged table.
REPORT_COLUMNS = [
    "id",
    "conflicted_rank",
    "ncbi",
    "bold",
    "gbif_filled",
    "resolved_from",
    "status",
]

GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
REQUEST_TIMEOUT = 30
HEADERS = {"User-Agent": "taxcord/1.0 (taxonomy merge backbone)"}

# Map each level to the GBIF match response field and rank parameter.
_GBIF_FIELD = {level: level.lower() for level in LEVELS}
_GBIF_RANK = {level: level.upper() for level in LEVELS}


def _consensus(ncbi_value, bold_value):
    """Return the agreed value, or NaN when the sources conflict."""
    ncbi_present, bold_present = pd.notna(ncbi_value), pd.notna(bold_value)
    if ncbi_present and bold_present:
        return ncbi_value if ncbi_value == bold_value else np.nan
    if ncbi_present:
        return ncbi_value
    if bold_present:
        return bold_value
    return np.nan


def _gbif_classification(session, name, level):
    """Return GBIF's accepted classification for ``name`` at ``level``, or None.

    None means the name could not be resolved at that rank, so the caller
    should fall back rather than trust an empty lineage.
    """
    try:
        response = session.get(
            GBIF_MATCH_URL,
            params={"name": name, "rank": _GBIF_RANK[level]},
            timeout=REQUEST_TIMEOUT,
        )
        data = response.json()
    except (requests.RequestException, ValueError):
        return None
    if data.get("matchType") in (None, "NONE") or not data.get("usageKey"):
        return None
    classification = {lvl: data.get(field) for lvl, field in _GBIF_FIELD.items()}
    # The queried rank must come back populated, else the anchor would vanish.
    if not classification.get(level):
        return None
    return classification


def gbif_resolver(session):
    """Build a memoised resolver: (name, level) -> accepted classification/None."""
    cache = {}

    def resolve(name, level):
        key = (level, name)
        if key not in cache:
            cache[key] = _gbif_classification(session, name, level)
        return cache[key]

    return resolve


def _reconcile_row(row, resolver, report=None):
    """Build the consensus lineage for one joined row.

    Without ``resolver`` (or when a row needs no reconciliation) this is the
    plain rank-wise consensus with downward gap propagation. When the sources
    agree at a finer rank but a coarser rank is unresolved, and a resolver is
    given, that coarser rank is filled from GBIF's backbone for the agreed
    taxon while the agreed ranks are kept; if GBIF cannot resolve it, the
    coarser rank is left blank rather than cascading down and losing the agreed
    anchor. Each such fill is appended to ``report`` when one is supplied.
    """
    out = {"id": row["id"]}
    pairs = {lvl: (row.get(f"{lvl}_ncbi"), row.get(f"{lvl}_bold")) for lvl in LEVELS}
    base = {lvl: _consensus(n, b) for lvl, (n, b) in pairs.items()}
    agreed = {
        lvl for lvl, (n, b) in pairs.items() if pd.notna(n) and pd.notna(b) and n == b
    }
    anchor_idx = max((LEVELS.index(lvl) for lvl in agreed), default=None)
    coarser_gap = anchor_idx is not None and any(
        pd.isna(base[LEVELS[i]]) for i in range(anchor_idx)
    )

    if resolver is not None and coarser_gap:
        anchor = LEVELS[anchor_idx]
        anchor_name = pairs[anchor][0]
        backbone = resolver(anchor_name, anchor)
        for i, lvl in enumerate(LEVELS):
            if i > anchor_idx:
                # Finer than the agreed anchor: not jointly supported.
                out[lvl] = np.nan
            elif pd.notna(base[lvl]):
                # Keep what the sources agree on (or a lone source value).
                out[lvl] = base[lvl]
            else:
                # Coarser rank the sources did not resolve: fill from GBIF if it
                # can, else leave blank (never cascade onto the agreed anchor).
                filled = backbone.get(lvl) if backbone else None
                out[lvl] = filled if filled is not None else np.nan
                if report is not None:
                    ncbi_value, bold_value = pairs[lvl]
                    report.append(
                        {
                            "id": row["id"],
                            "conflicted_rank": lvl,
                            "ncbi": "" if pd.isna(ncbi_value) else ncbi_value,
                            "bold": "" if pd.isna(bold_value) else bold_value,
                            "gbif_filled": "" if filled is None else filled,
                            "resolved_from": f"{anchor}={anchor_name}",
                            "status": "filled" if filled is not None else "unresolved",
                        }
                    )
        return out

    cleared = False
    for lvl in LEVELS:
        if pd.isna(base[lvl]):
            cleared = True
        out[lvl] = np.nan if cleared else base[lvl]
    return out


def merge_taxonomy(ncbi_df, bold_df, resolver=None, report=None):
    """Merge two condensed taxonomy tables into one consensus table.

    Pass ``resolver`` (see :func:`gbif_resolver`) to reconcile finer-rank
    agreements that a coarser-rank conflict would otherwise discard. Pass a
    list as ``report`` to collect one record per GBIF-filled rank (the columns
    in :data:`REPORT_COLUMNS`).
    """
    merged = pd.merge(
        ncbi_df, bold_df, on="id", how="outer", suffixes=("_ncbi", "_bold")
    )
    records = [_reconcile_row(row, resolver, report) for _, row in merged.iterrows()]
    return pd.DataFrame(records, columns=["id", *LEVELS])


def _report_path(output):
    """Path for the reconciliation report, alongside the merged output."""
    stem, _ = os.path.splitext(output)
    return f"{stem}.reconciliation.tsv"


def configure(parser):
    parser.add_argument("ncbi", help="NCBI condensed taxonomy CSV")
    parser.add_argument("bold", help="BOLD condensed taxonomy CSV")
    parser.add_argument("output", help="output CSV path")
    parser.add_argument(
        "--gbif-backbone",
        action="store_true",
        help="when the sources agree at a finer rank but conflict at a coarser "
        "one, resolve the agreed taxon via GBIF and fill the higher lineage "
        "from GBIF's accepted classification instead of discarding it (needs "
        "network access)",
    )


def execute(args):
    ncbi_df = pd.read_csv(args.ncbi)
    bold_df = pd.read_csv(args.bold)
    resolver = None
    report = None
    if args.gbif_backbone:
        session = requests.Session()
        session.headers.update(HEADERS)
        resolver = gbif_resolver(session)
        report = []
    result = merge_taxonomy(ncbi_df, bold_df, resolver, report)
    result.to_csv(args.output, index=False)

    message = f"Merged {len(result)} queries -> {args.output}"
    if report:
        report_path = _report_path(args.output)
        pd.DataFrame(report, columns=REPORT_COLUMNS).to_csv(
            report_path, sep="\t", index=False
        )
        message += f"; reconciled {len(report)} rank(s) via GBIF -> {report_path}"
    elif args.gbif_backbone:
        message += "; no ranks needed GBIF reconciliation"
    print(message)
