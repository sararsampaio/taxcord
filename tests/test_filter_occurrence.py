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
            "IP.species": ["12", "-"],
            "IP.genus": ["30", "-"],
            "IP.family": ["-", "5"],
            "IP.order": ["-", "-"],
            "BOLD.species": ["-", "-"],
            "BOLD.genus": ["-", "-"],
            "BOLD.family": ["-", "-"],
            "BOLD.order": ["-", "-"],
        }
    )
    result = filter_by_occurrence(df, source="bold").set_index("id")

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
            "IP.species": ["-"],
            "IP.genus": ["-"],
            "IP.family": ["-"],
            "IP.order": ["-"],
            "BOLD.species": ["-"],
            "BOLD.genus": ["-"],
            "BOLD.family": ["-"],
            "BOLD.order": ["-"],
        }
    )
    assert filter_by_occurrence(df, source="bold").empty


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
            "IP.species": ["7"],
            "IP.genus": ["-"],
            "IP.family": ["-"],
            "IP.order": ["-"],
            "BOLD.species": ["-"],
            "BOLD.genus": ["-"],
            "BOLD.family": ["-"],
            "BOLD.order": ["-"],
        }
    )
    result = filter_by_occurrence(df, source="ncbi")
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
