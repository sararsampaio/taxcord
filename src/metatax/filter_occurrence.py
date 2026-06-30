"""Trim taxonomic assignments to the finest rank with occurrence support.

Given a table that carries occurrence counts from GBIF (``GBIF.<rank>``) and
BOLD (``BOLD.<rank>``), keep each lineage only down to the finest rank that
either source actually recorded in the target region.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TAX_LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Ranks carrying occurrence counts, finest first.
OCCURRENCE_RANKS = ["species", "genus", "family", "order"]

# Canonical capitalised name for each taxonomic level we might receive, keyed by
# lower case. Upstream steps (condense/occurrences) emit lower-case ranks
# (``phylum``); raw NCBI annotation prefixes them with ``Taxonomy.``; a BOLD
# table may already be canonical (``Phylum``). All three are accepted.
_CANONICAL_LEVELS = {name.lower(): name for name in ["Kingdom", *TAX_LEVELS]}


def _rank_renames(columns):
    """Map each taxonomy column to its canonical name.

    Handles any letter case and an optional ``Taxonomy.`` prefix; leaves
    non-rank columns (``id``, ``GBIF.*``, ``BOLD.*``) untouched.
    """
    renames = {}
    for col in columns:
        base = col.split(".", 1)[1] if col.startswith("Taxonomy.") else col
        canonical = _CANONICAL_LEVELS.get(base.lower())
        if canonical and canonical != col:
            renames[col] = canonical
    return renames


def _supported_rank(row):
    """Return the finest rank with a recorded occurrence count, or None."""
    for rank in OCCURRENCE_RANKS:
        for prefix in ("GBIF", "BOLD"):
            # isdigit() is true only for a real count: it rejects "", "-" and NaN.
            if str(row.get(f"{prefix}.{rank}")).strip().isdigit():
                return rank
    return None


def filter_by_occurrence(df):
    """Keep rows with occurrence support and trim them to that rank.

    Args:
        df: Table with taxonomy columns plus ``GBIF.<rank>`` / ``BOLD.<rank>``.
            Rank columns are accepted in any case and with or without a
            ``Taxonomy.`` prefix, regardless of which source produced them.

    Returns:
        DataFrame with ``id`` and the six taxonomy levels, trimmed per row.
    """
    df = df.rename(columns=_rank_renames(df.columns)).copy()

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


def execute(args):
    df = pd.read_csv(args.input, sep="\t")
    result = filter_by_occurrence(df)
    result.to_csv(args.output, index=False)
    print(f"Kept {len(result)} supported queries -> {args.output}")
