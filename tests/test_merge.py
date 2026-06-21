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
