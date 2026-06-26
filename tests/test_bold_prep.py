import pandas as pd
import pytest

from metatax.bold_prep import reshape


def _boldigger_frame():
    """A frame shaped like a BOLDigger identification result."""
    return pd.DataFrame(
        {
            "id": ["OTU0001", "OTU0002", "OTU0003"],
            "Phylum": ["Arthropoda", "no-match", "Ascomycota"],
            "Class": ["Insecta", "no-match", "Leotiomycetes"],
            "Order": ["Diptera", "no-match", "Helotiales"],
            "Family": ["Chironomidae", "no-match", None],
            "Genus": [None, "no-match", None],
            "Species": [None, "no-match", None],
            # extra BOLDigger columns that must be dropped
            "pct_identity": [99.0, 0.0, 95.0],
            "status": ["public", None, "public"],
            "records": [12, None, 3],
        }
    )


def test_keeps_only_id_and_taxonomy_columns():
    out = reshape(_boldigger_frame())
    assert list(out.columns) == [
        "id",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]


def test_no_match_and_blank_cells_become_na():
    out = reshape(_boldigger_frame()).set_index("id")
    # a resolved row keeps its values, with blank sub-ranks as NA
    assert out.loc["OTU0001", "Family"] == "Chironomidae"
    assert out.loc["OTU0001", "Genus"] == "NA"
    # a fully unmatched row is all NA
    assert (out.loc["OTU0002"] == "NA").all()
    # a partially resolved row: blank below the resolved rank becomes NA
    assert out.loc["OTU0003", "Order"] == "Helotiales"
    assert out.loc["OTU0003", "Family"] == "NA"


def test_missing_expected_column_raises():
    df = _boldigger_frame().drop(columns=["Species"])
    with pytest.raises(SystemExit, match="Species"):
        reshape(df)
