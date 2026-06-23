"""
exporters.py
------------
Functions for exporting DataFrames to CSV and Excel (XLSX) formats.
"""

import io
import pandas as pd


def df_to_csv_bytes(df: pd.DataFrame) -> bytes:
    """Serialize a DataFrame to UTF-8 CSV bytes."""
    return df.to_csv(index=False).encode("utf-8")


def df_to_xlsx_bytes(df: pd.DataFrame, sheet_name: str = "Sheet1") -> bytes:
    """
    Serialize a DataFrame to XLSX bytes using openpyxl.

    Booleans are written as TRUE/FALSE strings so that Excel displays them
    correctly without converting to 1/0.
    """
    # Convert bool columns to strings for friendlier Excel display
    df_export = df.copy()
    for col in df_export.columns:
        if df_export[col].dtype == bool:
            df_export[col] = df_export[col].map({True: "TRUE", False: "FALSE"})

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df_export.to_excel(writer, index=False, sheet_name=sheet_name)
    return buffer.getvalue()


def analyses_export_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a clean export-ready copy of the Analyses_Sequences DataFrame.
    Column order: original CSV columns first, then computed columns.
    """
    # Identify computed/added columns
    computed_cols = [
        "AS_Sequence_valid",
        "Off-Targets_valid",
        "First_base_AU_valid",
        "AU_percentage",
        "AU_content_valid",
        "AU_1_14",
        "AU_1_14_percentage",
        "AU_1_14_percentage_valid",
        "AU_15_21",
        "AU_ratio",
        "AU_ratio_valid",
        "Position_13_14_valid",
        "No_A_position_20_valid",
        "No_forbidden_motifs_valid",
        "Final_valid",
        "Warnings",
    ]
    original_cols = [c for c in df.columns if c not in computed_cols]
    ordered_cols = original_cols + [c for c in computed_cols if c in df.columns]
    return df[ordered_cols]


def constructs_export_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a clean export-ready copy of the constructs DataFrame
    with the canonical column order.
    """
    desired = [
        "shRNA_Name",
        "Format",
        "Rank",
        "AS Sequence",
        "Antisense_Sequence",
        "Sense_Sequence",
        "LinkerA",
        "LinkerB",
        "Final_Construct",
        "Off-Targets",
        "Corrected_Score",
        "Position",
        "Warnings",
    ]
    available = [c for c in desired if c in df.columns]
    return df[available]
