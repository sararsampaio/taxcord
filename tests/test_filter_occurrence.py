import pandas as pd

from metatax.filter_occurrence import filter_by_occurrence


def test_bold_source_trims_to_supported_rank():
    df = pd.DataFrame(
        {
            "id": ["OTU_A", "OTU_B"],
            "Phylum": ["Arthropoda", "Arthropoda"],
            "Class": ["Insecta", "Insecta"],
            "Order": ["Diptera", "Diptera"],
            "Family": ["Mockidae", "Mockidae"],
            "Genus": ["Mockia", "Mockia"],
            "Species": ["Mockia alpha", "Mockia beta"],
            "GBIF.species": ["12", "-"],
            "GBIF.genus": ["30", "-"],
            "GBIF.family": ["-", "5"],
            "GBIF.order": ["-", "-"],
            "BOLD.species": ["-", "-"],
            "BOLD.genus": ["-", "-"],
            "BOLD.family": ["-", "-"],
            "BOLD.order": ["-", "-"],
        }
    )
    result = filter_by_occurrence(df).set_index("id")

    # OTU_A is supported to species; OTU_B only to family, so finer ranks clear.
    assert result.loc["OTU_A", "Species"] == "Mockia alpha"
    assert result.loc["OTU_B", "Family"] == "Mockidae"
    assert pd.isna(result.loc["OTU_B", "Genus"])
    assert pd.isna(result.loc["OTU_B", "Species"])


def test_unsupported_rows_are_dropped():
    df = pd.DataFrame(
        {
            "id": ["OTU_C"],
            "Phylum": ["Arthropoda"],
            "Class": ["Insecta"],
            "Order": ["Diptera"],
            "Family": ["Mockidae"],
            "Genus": ["Mockia"],
            "Species": ["Mockia gamma"],
            "GBIF.species": ["-"],
            "GBIF.genus": ["-"],
            "GBIF.family": ["-"],
            "GBIF.order": ["-"],
            "BOLD.species": ["-"],
            "BOLD.genus": ["-"],
            "BOLD.family": ["-"],
            "BOLD.order": ["-"],
        }
    )
    assert filter_by_occurrence(df).empty


def test_ncbi_source_renames_headers():
    df = pd.DataFrame(
        {
            "id": ["OTU_A"],
            "Taxonomy.kingdom": ["Metazoa"],
            "Taxonomy.phylum": ["Arthropoda"],
            "Taxonomy.class": ["Insecta"],
            "Taxonomy.order": ["Diptera"],
            "Taxonomy.family": ["Mockidae"],
            "Taxonomy.genus": ["Mockia"],
            "Taxonomy.species": ["Mockia alpha"],
            "GBIF.species": ["7"],
            "GBIF.genus": ["-"],
            "GBIF.family": ["-"],
            "GBIF.order": ["-"],
            "BOLD.species": ["-"],
            "BOLD.genus": ["-"],
            "BOLD.family": ["-"],
            "BOLD.order": ["-"],
        }
    )
    result = filter_by_occurrence(df)
    assert list(result.columns) == [
        "id",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    assert result.iloc[0]["Species"] == "Mockia alpha"
