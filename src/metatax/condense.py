"""Collapse BLAST hits into one taxonomic lineage per query.

A single query sequence (OTU) typically produces many BLAST hits spanning
several taxa. This module walks up the rank hierarchy and assigns the finest
rank at which the hits agree, given their percent identity.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

EMPTY_VALUES = {"", " ", "NA", None}

# Taxonomic ranks ordered from most general to most specific. The actual ranks
# used at runtime are read from the input header; this is the conventional set.
DEFAULT_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Fraction of voting hits that must share a value for a rank to be assigned.
AGREEMENT = 0.8

DEFAULT_QCOVS_MIN = 91.0


@dataclass(frozen=True)
class Tier:
    """Identity band that decides how finely a query can be resolved.

    Attributes:
        trigger: Apply this tier when the best hit's identity is at least this.
        floor: Only hits with identity at least this value vote on the call.
        rank: Finest rank this tier attempts to resolve.
    """

    trigger: float
    floor: float
    rank: str


# Default identity bands. The best hit's identity selects a tier; the tier's
# floor selects the voting hits and its rank is where resolution starts.
DEFAULT_TIERS = [
    Tier(trigger=99, floor=97, rank="species"),
    Tier(trigger=97, floor=97, rank="genus"),
    Tier(trigger=90, floor=90, rank="family"),
    Tier(trigger=85, floor=85, rank="order"),
    Tier(trigger=0, floor=0, rank="class"),
]


@dataclass
class Hit:
    """One BLAST hit reduced to the fields needed for assignment."""

    identity: float
    coverage: float
    lineage: dict  # rank name -> value, with missing values normalised to None


def _clean(value):
    """Normalise a raw cell to a value or None."""
    if isinstance(value, str):
        value = value.strip()
    return None if value in EMPTY_VALUES else value


def select_tier(top_identity, tiers):
    """Return the first tier whose trigger the best hit meets."""
    for tier in tiers:
        if top_identity >= tier.trigger:
            return tier
    return tiers[-1]


def _missing_agree_upward(missing, representative, rank, ranks):
    """Check that hits lacking a value at `rank` still agree higher up.

    For each hit missing `rank`, find its nearest populated higher rank and
    compare it with the representative hit. Returns True only when more than
    AGREEMENT of the missing hits agree.
    """
    higher = ranks[: ranks.index(rank)]  # general -> specific
    agree = 0
    for hit in missing:
        for higher_rank in reversed(higher):  # nearest higher rank first
            value = hit.lineage.get(higher_rank)
            if value is not None:
                if value == representative.lineage.get(higher_rank):
                    agree += 1
                break
    return bool(missing) and agree / len(missing) > AGREEMENT


def _resolve_rank(hits, rank, ranks):
    """Return a hit representing the agreed value at `rank`, or None.

    A rank resolves when one value covers at least AGREEMENT of the hits that
    have any value there. When most hits lack a value, the few that have one
    must also be backed by agreement at a higher rank.
    """
    populated = [hit for hit in hits if hit.lineage.get(rank) is not None]
    if not populated:
        return None

    winner, count = Counter(hit.lineage[rank] for hit in populated).most_common(1)[0]
    if count / len(populated) < AGREEMENT:
        return None

    representative = next(hit for hit in populated if hit.lineage[rank] == winner)

    missing = [hit for hit in hits if hit.lineage.get(rank) is None]
    if len(missing) / len(hits) >= 0.5 and not _missing_agree_upward(
        missing, representative, rank, ranks
    ):
        return None

    return representative


def assign_lineage(hits, ranks=DEFAULT_RANKS, tiers=DEFAULT_TIERS):
    """Assign one lineage to a query from its BLAST hits.

    Returns a dict mapping every rank to a value or None. Ranks below the
    resolved level are left as None.
    """
    hits = sorted(hits, key=lambda hit: hit.identity, reverse=True)
    top_identity = hits[0].identity
    tier = select_tier(top_identity, tiers)

    voting = [hit for hit in hits if hit.identity >= tier.floor]
    if len(voting) <= 1:
        # Too few hits to vote: widen the band to 2% below the best hit.
        voting = [hit for hit in hits if hit.identity >= top_identity - 2]

    start = ranks.index(tier.rank) if tier.rank in ranks else len(ranks) - 1
    for idx in range(start, -1, -1):
        representative = _resolve_rank(voting, ranks[idx], ranks)
        if representative is not None:
            return {
                rank: (representative.lineage.get(rank) if i <= idx else None)
                for i, rank in enumerate(ranks)
            }

    return {rank: None for rank in ranks}


def read_blast_hits(path, qcovs_min=DEFAULT_QCOVS_MIN):
    """Read a pipe-delimited annotated BLAST file into hits grouped by query.

    The file is the output of the R annotation step: columns separated by
    ``|`` with ``id``, ``pident``, ``qcovs`` and one ``Taxonomy.<rank>`` column
    per rank. R's ``write.table`` prefixes each data row with a row number,
    which is detected and dropped.

    Returns (otus, ranks) where otus maps query id -> list[Hit].
    """
    otus = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("|")
        column = {name: i for i, name in enumerate(header)}
        rank_columns = {
            name.split(".", 1)[1]: i
            for name, i in column.items()
            if name.startswith("Taxonomy.")
        }
        ranks = list(rank_columns)
        id_col, pident_col, qcovs_col = column["id"], column["pident"], column["qcovs"]

        for line in handle:
            fields = line.rstrip("\n").split("|")
            if len(fields) == len(header) + 1:
                fields = fields[1:]  # drop R row-number column
            if len(fields) != len(header):
                continue  # skip malformed rows

            coverage = float(fields[qcovs_col])
            if coverage < qcovs_min:
                continue

            lineage = {rank: _clean(fields[i]) for rank, i in rank_columns.items()}
            hit = Hit(float(fields[pident_col]), coverage, lineage)
            otus.setdefault(fields[id_col], []).append(hit)

    return otus, ranks


def write_assignments(otus, ranks, out_path, tiers=DEFAULT_TIERS):
    """Assign a lineage to every query and write a tab-delimited table."""
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(["id", *ranks]) + "\n")
        for query_id, hits in otus.items():
            lineage = assign_lineage(hits, ranks, tiers)
            out.write(
                "\t".join([query_id] + [lineage[rank] or "NA" for rank in ranks]) + "\n"
            )


def configure(parser):
    parser.add_argument("input", help="pipe-delimited annotated BLAST file")
    parser.add_argument("output", help="tab-delimited condensed taxonomy file")
    parser.add_argument(
        "--qcovs-min",
        type=float,
        default=DEFAULT_QCOVS_MIN,
        help="minimum query coverage to keep a hit (default: %(default)s)",
    )


def execute(args):
    otus, ranks = read_blast_hits(args.input, args.qcovs_min)
    write_assignments(otus, ranks, args.output)
    print(f"Condensed {len(otus)} queries -> {args.output}")
