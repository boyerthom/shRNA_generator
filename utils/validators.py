"""
validators.py
-------------
Validates CSV structure and individual AS Sequence entries.
"""

import pandas as pd

REQUIRED_CSV_COLUMNS = ["AS Sequence", "Off-Targets", "Corrected_Score", "Position"]
VALID_SEQ_BASES = set("AUTCG")


def validate_csv_columns(df: pd.DataFrame) -> list:
    """
    Check that all required columns exist in the uploaded CSV.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    list of error strings (empty if valid).
    """
    missing = [col for col in REQUIRED_CSV_COLUMNS if col not in df.columns]
    if missing:
        return [f"CSV is missing required column(s): {', '.join(missing)}"]
    return []


def validate_as_sequence(seq) -> tuple[bool, str]:
    """
    Validate a single AS Sequence entry.

    Rules:
    - Must be exactly 21 nucleotides
    - Must contain only A, U, T, C, G (case-insensitive)

    Parameters
    ----------
    seq : any
        Raw value from the AS Sequence column.

    Returns
    -------
    (is_valid: bool, warning_message: str)
    """
    if not isinstance(seq, str):
        return False, f"AS Sequence is not a string (got {type(seq).__name__})."

    seq_upper = seq.upper().strip()

    if len(seq_upper) != 21:
        return False, f"AS Sequence '{seq_upper}' has length {len(seq_upper)} (expected 21)."

    invalid = set(seq_upper) - VALID_SEQ_BASES
    if invalid:
        return False, (
            f"AS Sequence '{seq_upper}' contains invalid characters: "
            f"{', '.join(sorted(invalid))}."
        )

    return True, ""
