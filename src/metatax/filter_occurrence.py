"""Trim taxonomic assignments to the finest rank with occurrence support.

Given a table that carries occurrence counts from GBIF (``IP.<rank>``) and
BOLD (``BOLD.<rank>``), keep each lineage only down to the finest rank that
either source actually recorded in the target region.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TAX_LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Ranks carrying occurrence counts, finest first.
OCCURRENCE_RANKS = ["species", "genus", "family", "order"]

# Maps the NCBI condensed-file headers onto the canonical level names.
NCBI_RENAME = {
    "Taxonomy.kingdom": "Kingdom",
    "Taxonomy.phylum": "Phylum",
    "Taxonomy.class": "Class",
    "Taxonomy.order": "Order",
    "Taxonomy.family": "Family",
    "Taxonomy.genus": "Genus",
    "Taxonomy.species": "Species",
}


def _supported_rank(row):
    """Return the finest rank with a recorded occurrence count, or None."""
    for rank in OCCURRENCE_RANKS:
        for prefix in ("IP", "BOLD"):
            # isdigit() is true only for a real count: it rejects "", "-" and NaN.
            if str(row.get(f"{prefix}.{rank}")).strip().isdigit():
                return rank
    return None


def filter_by_occurrence(df, source):
    """Keep rows with occurrence support and trim them to that rank.

    Args:
        df: Table with taxonomy columns plus ``IP.<rank>`` / ``BOLD.<rank>``.
        source: ``"ncbi"`` renames ``Taxonomy.<rank>`` headers first;
            ``"bold"`` expects canonical level names already.

    Returns:
        DataFrame with ``id`` and the six taxonomy levels, trimmed per row.
    """
    if source == "ncbi":
        df = df.rename(columns=NCBI_RENAME)
    df = df.copy()

    lowest = df.apply(_supported_rank, axis=1)
    supported = lowest.notna()
    df, lowest = df[supported].copy(), lowest[supported]

    cutoff = lowest.map(lambda rank: TAX_LEVELS.index(rank.capitalize()))
    for depth, level in enumerate(TAX_LEVELS):
        df[level] = [
            value if depth <= limit else np.nan
            for value, limit in zip(df[level], cutoff)
        ]

    return df[["id", *TAX_LEVELS]]


def configure(parser):
    parser.add_argument("input", help="tab-delimited table with occurrence counts")
    parser.add_argument("output", help="output CSV path")
    parser.add_argument(
        "--source",
        choices=("ncbi", "bold"),
        required=True,
        help="header style of the input table",
    )


def execute(args):
    df = pd.read_csv(args.input, sep="\t")
    result = filter_by_occurrence(df, args.source)
    result.to_csv(args.output, index=False)
    print(f"Kept {len(result)} supported queries -> {args.output}")
