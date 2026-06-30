from pathlib import Path

import pytest

from metatax.condense import (
    DEFAULT_RANKS,
    DEFAULT_TIERS,
    Hit,
    _resolve_tiers,
    assign_lineage,
    read_blast_hits,
    write_assignments,
)

FIXTURE = Path(__file__).parent / "fixtures" / "annotated_blast.txt"


def make_hit(identity, **lineage):
    full = {rank: lineage.get(rank) for rank in DEFAULT_RANKS}
    return Hit(identity=identity, coverage=100.0, lineage=full)


def test_unanimous_hits_resolve_to_species():
    hits = [make_hit(99.0, genus="Mockia", species="Mockia alpha") for _ in range(3)]
    result = assign_lineage(hits)
    assert result["species"] == "Mockia alpha"


def test_species_disagreement_falls_back_to_genus():
    hits = [
        make_hit(99.0, genus="Mockia", species="Mockia beta"),
        make_hit(98.5, genus="Mockia", species="Mockia beta"),
        make_hit(97.5, genus="Mockia", species="Mockia gamma"),
    ]
    result = assign_lineage(hits)
    assert result["genus"] == "Mockia"
    assert result["species"] is None


def test_no_agreement_returns_empty_lineage():
    hits = [
        make_hit(99.0, kingdom="Metazoa", phylum="Arthropoda"),
        make_hit(99.0, kingdom="Plantae", phylum="Streptophyta"),
    ]
    # Two equally common phyla -> neither reaches the agreement threshold, and
    # the conflicting kingdoms above it also fail.
    result = assign_lineage(hits)
    assert all(value is None for value in result.values())


def test_sparse_lower_rank_requires_higher_agreement():
    # Only one of four hits names a genus; the other three are blank there but
    # all agree on family, so the genus call is allowed to stand.
    hits = [
        make_hit(99.0, family="Mockidae", genus="Mockia"),
        make_hit(99.0, family="Mockidae"),
        make_hit(99.0, family="Mockidae"),
        make_hit(99.0, family="Mockidae"),
    ]
    result = assign_lineage(hits)
    assert result["genus"] == "Mockia"


def test_read_blast_drops_row_numbers_and_filters_coverage():
    otus, ranks = read_blast_hits(FIXTURE)
    assert ranks == DEFAULT_RANKS
    assert set(otus) == {"OTU_A", "OTU_B"}  # OTU_C dropped: qcovs 80 < 91
    assert len(otus["OTU_A"]) == 3


def test_write_assignments_round_trip(tmp_path):
    otus, ranks = read_blast_hits(FIXTURE)
    out = tmp_path / "condensed.txt"
    write_assignments(otus, ranks, out)

    rows = {
        line.split("\t")[0]: line.rstrip("\n").split("\t")
        for line in out.read_text(encoding="utf-8").splitlines()[1:]
    }
    assert rows["OTU_A"][-1] == "Mockia alpha"
    assert rows["OTU_B"][ranks.index("genus") + 1] == "Mockia"
    assert rows["OTU_B"][-1] == "NA"


def test_resolve_tiers_defaults_when_no_override():
    assert _resolve_tiers(None) == DEFAULT_TIERS
    assert _resolve_tiers([]) == DEFAULT_TIERS


def test_resolve_tiers_overrides_only_named_rank():
    tiers = _resolve_tiers(["species:98:95"])
    by_rank = {tier.rank: tier for tier in tiers}
    # the overridden rank takes the new thresholds...
    assert by_rank["species"].trigger == 98
    assert by_rank["species"].floor == 95
    # ...while the others keep their defaults
    assert by_rank["family"].trigger == 90
    # result stays ordered by trigger, highest first
    assert [tier.trigger for tier in tiers] == sorted(
        (tier.trigger for tier in tiers), reverse=True
    )


def test_resolve_tiers_floor_defaults_to_trigger():
    (tier,) = [t for t in _resolve_tiers(["genus:96"]) if t.rank == "genus"]
    assert tier.trigger == 96
    assert tier.floor == 96


def test_resolve_tiers_rejects_bad_spec():
    with pytest.raises(SystemExit):
        _resolve_tiers(["species"])
    with pytest.raises(SystemExit):
        _resolve_tiers(["species:abc"])
