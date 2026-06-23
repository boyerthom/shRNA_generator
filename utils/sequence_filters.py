"""
sequence_filters.py
-------------------
Computes all per-row analysis columns and validation flags for the AS Sequence table.
"""

import math
import pandas as pd
from utils.validators import validate_as_sequence


# Nucleotides treated as A/U (AU-like)
AU_BASES = {"A", "U", "T"}

FORBIDDEN_MOTIFS = ["AAAAAA", "TTTTT", "UUUUU", "CCCC", "GGGG"]


def _normalize(seq: str) -> str:
    """Convert sequence to uppercase. U and T are treated equivalently downstream."""
    return seq.upper().strip()


def _au_count(seq: str) -> int:
    """Count A + U/T bases in a sequence."""
    return sum(1 for b in seq if b in AU_BASES)


def analyze_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add all analysis and validation columns to the DataFrame.

    This function operates row-by-row on the 'AS Sequence' column.
    All original columns are preserved; computed columns are appended.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'AS Sequence', 'Off-Targets', 'Corrected_Score', 'Position'.

    Returns
    -------
    pd.DataFrame with additional columns.
    """
    results = []

    for _, row in df.iterrows():
        r = row.to_dict()
        warnings_list = []

        raw_seq = r.get("AS Sequence", "")

        # --- AS Sequence validation ---
        is_seq_valid, seq_warning = validate_as_sequence(raw_seq)
        if not is_seq_valid:
            warnings_list.append(seq_warning)
            # Fill all computed columns with None/False and skip analysis
            r.update(_empty_analysis(warnings_list))
            results.append(r)
            continue

        seq = _normalize(raw_seq)  # uppercase

        # --- 1. Off-Targets validation ---
        try:
            off_targets = float(r.get("Off-Targets", math.inf))
            off_targets_valid = off_targets <= 1
        except (ValueError, TypeError):
            off_targets_valid = False
            warnings_list.append("Off-Targets value could not be parsed as a number.")

        # --- 2. First base A/U validation ---
        first_base = seq[0]
        first_base_au_valid = first_base in AU_BASES

        # --- 3. AU content ---
        au_count_total = _au_count(seq)
        au_percentage = au_count_total / len(seq)
        au_content_valid = 0.40 <= au_percentage <= 0.80

        # --- 4 & 5. AU counts by region ---
        au_1_14 = _au_count(seq[0:14])   # positions 1–14 (0-indexed 0:14)
        au_15_21 = _au_count(seq[14:21]) # positions 15–21 (0-indexed 14:21)

        # --- 4b. AU percentage in positions 1–14 ---
        au_1_14_pct = au_1_14 / 14
        au_1_14_percentage_valid = au_1_14_pct > 0.50

        # --- 6. AU ratio ---
        if au_15_21 == 0:
            au_ratio = math.inf
            au_ratio_valid = True  # infinity > 1
            warnings_list.append(
                f"AU_15_21 = 0 for sequence '{seq}'; AU_ratio set to infinity."
            )
        else:
            au_ratio = (au_1_14 / 14) / (au_15_21 / 7)
            au_ratio_valid = au_ratio > 1

        # --- 7. Position 13/14 rule ---
        pos13 = seq[12]  # 0-indexed position 12 → nucleotide 13
        pos14 = seq[13]  # 0-indexed position 13 → nucleotide 14
        pos_13_14_valid = (pos13 in AU_BASES) or (pos14 in {"U", "T"})

        # --- 8. Position 20 rule ---
        pos20 = seq[19]  # 0-indexed position 19 → nucleotide 20
        no_a_pos20_valid = pos20 != "A"

        # --- 9. Forbidden motifs ---
        seq_for_motif = seq.replace("U", "T")  # normalize U→T before motif check
        has_forbidden = any(motif.replace("U", "T") in seq_for_motif for motif in FORBIDDEN_MOTIFS)
        no_forbidden_motifs_valid = not has_forbidden

        # --- 10. Final validity ---
        final_valid = all([
            is_seq_valid,
            off_targets_valid,
            first_base_au_valid,
            au_content_valid,
            au_1_14_percentage_valid,
            au_ratio_valid,
            pos_13_14_valid,
            no_a_pos20_valid,
            no_forbidden_motifs_valid,
        ])

        r.update({
            "AS_Sequence_valid": is_seq_valid,
            "Off-Targets_valid": off_targets_valid,
            "First_base_AU_valid": first_base_au_valid,
            "AU_percentage": round(au_percentage * 100, 2),
            "AU_content_valid": au_content_valid,
            "AU_1_14": au_1_14,
            "AU_1_14_percentage": round(au_1_14_pct * 100, 2),
            "AU_1_14_percentage_valid": au_1_14_percentage_valid,
            "AU_15_21": au_15_21,
            "AU_ratio": au_ratio if not math.isinf(au_ratio) else "∞",
            "AU_ratio_valid": au_ratio_valid,
            "Position_13_14_valid": pos_13_14_valid,
            "No_A_position_20_valid": no_a_pos20_valid,
            "No_forbidden_motifs_valid": no_forbidden_motifs_valid,
            "Final_valid": final_valid,
            "Warnings": "; ".join(warnings_list) if warnings_list else "",
        })

        results.append(r)

    return pd.DataFrame(results)


def _empty_analysis(warnings_list: list) -> dict:
    """Return a dict of None/False values for all computed columns when AS Sequence is invalid."""
    return {
        "AS_Sequence_valid": False,
        "Off-Targets_valid": False,
        "First_base_AU_valid": False,
        "AU_percentage": None,
        "AU_content_valid": False,
        "AU_1_14": None,
        "AU_1_14_percentage": None,
        "AU_1_14_percentage_valid": False,
        "AU_15_21": None,
        "AU_ratio": None,
        "AU_ratio_valid": False,
        "Position_13_14_valid": False,
        "No_A_position_20_valid": False,
        "No_forbidden_motifs_valid": False,
        "Final_valid": False,
        "Warnings": "; ".join(warnings_list),
    }
