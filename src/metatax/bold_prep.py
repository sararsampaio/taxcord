"""Reshape a BOLDigger identification table into a condensed lineage file.

BOLDigger returns one identification per query (OTU), so it is the BOLD-branch
counterpart of the ``condense`` step rather than an input to it: there are no
multiple hits to collapse. This selects the ``id`` and taxonomy columns,
normalises BOLDigger's ``no-match`` sentinel and blank sub-rank cells to
``NA``, and writes the tab-delimited table that ``occurrences`` consumes.
"""

from __future__ import annotations

import sys

import pandas as pd

# Taxonomy columns expected from BOLDigger, coarsest first.
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

# BOLDigger marks an unidentified query with this sentinel across every rank.
NO_MATCH = "no-match"


def reshape(df):
    """Return ``id`` + taxonomy columns with no-match/blank cells set to ``NA``.

    Extra BOLDigger columns (pct_identity, status, records, ...) are dropped so
    ``occurrences`` does not mistake them for ranks.
    """
    missing = [name for name in ["id", *RANKS] if name not in df.columns]
    if missing:
        raise SystemExit(
            f"metatax bold-prep: input is missing expected column(s) "
            f"{', '.join(missing)}. Found: {', '.join(map(str, df.columns))}"
        )
    out = df[["id", *RANKS]].copy()
    out[RANKS] = out[RANKS].replace(NO_MATCH, "NA").fillna("NA")
    return out


def _read_table(path):
    """Read a BOLDigger results table from Excel or delimited text."""
    lowered = path.lower()
    if lowered.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    sep = "\t" if lowered.endswith((".txt", ".tsv")) else ","
    return pd.read_csv(path, sep=sep)


def configure(parser):
    parser.add_argument(
        "input", help="BOLDigger identification table (.xlsx, .csv or .tsv)"
    )
    parser.add_argument("output", help="tab-delimited condensed lineage file")


def execute(args):
    out = reshape(_read_table(args.input))

    duplicates = int(out["id"].duplicated().sum())
    if duplicates:
        raise SystemExit(
            f"metatax bold-prep: {duplicates} duplicate id(s) in input; "
            f"expected one row per OTU"
        )

    out.to_csv(args.output, sep="\t", index=False)
    resolved = int(out[RANKS].ne("NA").any(axis=1).sum())
    print(
        f"Prepared {len(out)} OTUs ({resolved} with a taxon) -> {args.output}",
        file=sys.stderr,
    )
