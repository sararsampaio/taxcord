import pandas as pd

from metatax.merge import merge_taxonomy

LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _frame(id_, values):
    return pd.DataFrame([{"id": id_, **dict(zip(LEVELS, values))}])


def test_agreement_is_kept_and_conflict_clears_rank():
    ncbi = _frame(
        "OTU_A",
        ["Arthropoda", "Insecta", "Diptera", "Mockidae", "Mockia", "Mockia alpha"],
    )
    bold = _frame(
        "OTU_A",
        ["Arthropoda", "Insecta", "Diptera", "Mockidae", "Mockia", "Mockia beta"],
    )
    result = merge_taxonomy(ncbi, bold).set_index("id")

    assert result.loc["OTU_A", "Genus"] == "Mockia"
    # Species conflicts -> cleared, and nothing sits below it to propagate.
    assert pd.isna(result.loc["OTU_A", "Species"])


def test_gap_propagates_to_lower_ranks():
    ncbi = _frame(
        "OTU_B", ["Arthropoda", "Insecta", "Diptera", None, "Mockia", "Mockia beta"]
    )
    bold = _frame(
        "OTU_B", ["Arthropoda", "Insecta", "Diptera", None, "Mockia", "Mockia beta"]
    )
    result = merge_taxonomy(ncbi, bold).set_index("id")

    assert result.loc["OTU_B", "Order"] == "Diptera"
    # Family is empty, so genus and species are cleared even though they agreed.
    assert pd.isna(result.loc["OTU_B", "Family"])
    assert pd.isna(result.loc["OTU_B", "Genus"])
    assert pd.isna(result.loc["OTU_B", "Species"])


def test_single_source_value_is_used():
    ncbi = _frame(
        "OTU_C",
        ["Arthropoda", "Insecta", "Diptera", "Mockidae", "Mockia", "Mockia gamma"],
    )
    bold = _frame(
        "OTU_C", ["Arthropoda", "Insecta", "Diptera", "Mockidae", "Mockia", None]
    )
    result = merge_taxonomy(ncbi, bold).set_index("id")

    assert result.loc["OTU_C", "Species"] == "Mockia gamma"


# A backbone for the OTU8-style case: sources agree on the species but disagree
# on family; GBIF places the species in Helotiaceae.
_BACKBONE = {
    "Articulospora tetracladia": {
        "Phylum": "Ascomycota",
        "Class": "Leotiomycetes",
        "Order": "Helotiales",
        "Family": "Helotiaceae",
        "Genus": "Articulospora",
        "Species": "Articulospora tetracladia",
    }
}


def _fake_resolver(name, level):
    return _BACKBONE.get(name)


def test_family_conflict_loses_species_without_backbone():
    # Same species, different family -> family conflict cascades and clears species.
    ncbi = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Discinellaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    bold = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Helotiaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    result = merge_taxonomy(ncbi, bold).set_index("id")
    assert pd.isna(result.loc["OTU8", "Family"])
    assert pd.isna(result.loc["OTU8", "Species"])


def test_gbif_backbone_keeps_species_and_fills_family():
    ncbi = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Discinellaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    bold = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Helotiaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    result = merge_taxonomy(ncbi, bold, resolver=_fake_resolver).set_index("id")
    assert result.loc["OTU8", "Species"] == "Articulospora tetracladia"
    assert result.loc["OTU8", "Family"] == "Helotiaceae"  # from GBIF backbone


def test_backbone_fallback_keeps_species_blank_family_when_unresolved():
    ncbi = _frame(
        "OTU9",
        ["Ascomycota", "Leotiomycetes", "Helotiales", "Famone", "Genusx", "Genusx sp"],
    )
    bold = _frame(
        "OTU9",
        ["Ascomycota", "Leotiomycetes", "Helotiales", "Famtwo", "Genusx", "Genusx sp"],
    )
    # resolver returns None (GBIF cannot resolve) -> graceful fallback
    result = merge_taxonomy(ncbi, bold, resolver=lambda name, level: None).set_index(
        "id"
    )
    assert result.loc["OTU9", "Species"] == "Genusx sp"  # kept
    assert pd.isna(result.loc["OTU9", "Family"])  # conflicted family left blank


def test_report_records_each_reconciled_rank():
    ncbi = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Discinellaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    bold = _frame(
        "OTU8",
        [
            "Ascomycota",
            "Leotiomycetes",
            "Helotiales",
            "Helotiaceae",
            "Articulospora",
            "Articulospora tetracladia",
        ],
    )
    report = []
    merge_taxonomy(ncbi, bold, resolver=_fake_resolver, report=report)

    assert len(report) == 1
    entry = report[0]
    assert entry["id"] == "OTU8"
    assert entry["conflicted_rank"] == "Family"
    assert entry["ncbi"] == "Discinellaceae"
    assert entry["bold"] == "Helotiaceae"
    assert entry["gbif_filled"] == "Helotiaceae"
    assert entry["status"] == "filled"
    assert entry["resolved_from"] == "Species=Articulospora tetracladia"
