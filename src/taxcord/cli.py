"""Command-line entry point dispatching to the pipeline steps."""

from __future__ import annotations

import argparse

from . import bold_prep, condense, filter_occurrence, merge, occurrences

COMMANDS = {
    "condense": (condense, "collapse BLAST hits into one lineage per query"),
    "bold-prep": (bold_prep, "reshape a BOLDigger table into a condensed lineage file"),
    "occurrences": (occurrences, "annotate lineages with GBIF/BOLD record counts"),
    "filter": (filter_occurrence, "trim lineages to the occurrence-supported rank"),
    "merge": (merge, "merge NCBI and BOLD tables into a consensus lineage"),
}


def build_parser():
    parser = argparse.ArgumentParser(
        prog="taxcord",
        description="Taxonomic assignment and refinement for metabarcoding data.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, (module, help_text) in COMMANDS.items():
        subparser = subparsers.add_parser(name, help=help_text)
        module.configure(subparser)
        subparser.set_defaults(handler=module.execute)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.handler(args)


if __name__ == "__main__":
    main()
