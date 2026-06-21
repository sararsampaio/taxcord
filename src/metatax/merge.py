"""Merge NCBI and BOLD taxonomy tables into one consensus lineage per query.

The two sources are joined on query id. At each rank the call is the value
they agree on; a disagreement clears that rank. Once a rank is empty the lower
ranks are cleared too, so a lineage never claims more detail than it supports.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


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


def _propagate_gaps(row):
    """Clear every rank below the first empty one."""
    cleared = False
    for level in LEVELS:
        if pd.isna(row[level]):
            cleared = True
        if cleared:
            row[level] = np.nan
    return row


def merge_taxonomy(ncbi_df, bold_df):
    """Merge two condensed taxonomy tables into one consensus table."""
    merged = pd.merge(
        ncbi_df, bold_df, on="id", how="outer", suffixes=("_ncbi", "_bold")
    )
    for level in LEVELS:
        merged[level] = [
            _consensus(ncbi_value, bold_value)
            for ncbi_value, bold_value in zip(
                merged[f"{level}_ncbi"], merged[f"{level}_bold"]
            )
        ]
    merged = merged.apply(_propagate_gaps, axis=1)
    return merged[["id", *LEVELS]]


def configure(parser):
    parser.add_argument("ncbi", help="NCBI condensed taxonomy CSV")
    parser.add_argument("bold", help="BOLD condensed taxonomy CSV")
    parser.add_argument("output", help="output CSV path")


def execute(args):
    ncbi_df = pd.read_csv(args.ncbi)
    bold_df = pd.read_csv(args.bold)
    result = merge_taxonomy(ncbi_df, bold_df)
    result.to_csv(args.output, index=False)
    print(f"Merged {len(result)} queries -> {args.output}")
