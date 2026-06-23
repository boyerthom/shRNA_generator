"""
sequence_builder.py
-------------------
Constructs 90mer or 97mer shRNA sequences for all Final_valid candidates.
"""

import math
import pandas as pd

# --- Fixed structural components ---
LOOP = "TAGTGAAGCCACAGATGTA"

# 90mer flanking sequences
FIVE_PRIME_END_90  = "GTTGACAGTGAGCG"
THREE_PRIME_END_90 = "TGCCTACTGCCTC"

# 97mer flanking sequences (original miR30 design, Dow et al.)
FIVE_PRIME_END_97  = "TGCTGTTGACAGTGAGCG"
THREE_PRIME_END_97 = "TGCCTACTGCCTCGGA"

FORMAT_COMPONENTS = {
    "90mer": (FIVE_PRIME_END_90, THREE_PRIME_END_90),
    "97mer": (FIVE_PRIME_END_97, THREE_PRIME_END_97),
}

COMPLEMENT = str.maketrans("ATCG", "TAGC")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (A↔T, C↔G)."""
    return seq.translate(COMPLEMENT)[::-1]


def get_linker_b(sense_seq: str, clean_fasta: str) -> tuple[str, str]:
    """
    Find the base immediately upstream of the sense sequence in the FASTA.

    Parameters
    ----------
    sense_seq   : str – the 21-nt sense sequence (DNA, uppercase)
    clean_fasta : str – the cleaned FASTA string (DNA, uppercase)

    Returns
    -------
    (linker_b: str, warning: str)
        linker_b is "" if it cannot be determined.
    """
    matches = []
    start = 0
    while True:
        idx = clean_fasta.find(sense_seq, start)
        if idx == -1:
            break
        matches.append(idx)
        start = idx + 1

    if len(matches) == 0:
        return "", "Sense sequence not found in FASTA."

    if len(matches) > 1:
        return "", "Multiple FASTA matches found."

    pos = matches[0]
    if pos == 0:
        return "", "No upstream nucleotide available."

    return clean_fasta[pos - 1], ""


def get_linker_a(linker_b: str) -> tuple[str, str]:
    """
    Determine LinkerA based on LinkerB.

    Rules:
        LinkerB ∈ {A, T, U} → LinkerA = C
        LinkerB ∈ {C, G}    → LinkerA = A
        Otherwise            → warning

    Returns
    -------
    (linker_a: str, warning: str)
    """
    if not linker_b:
        return "", "LinkerA cannot be determined: LinkerB is missing."

    lb = linker_b.upper()
    if lb in {"A", "T", "U"}:
        return "C", ""
    elif lb in {"C", "G"}:
        return "A", ""
    else:
        return "", f"LinkerA cannot be determined: unexpected LinkerB value '{linker_b}'."


def build_constructs(
    df: pd.DataFrame,
    clean_fasta: str,
    project_name: str = "",
    construct_format: str = "97mer",
) -> pd.DataFrame:
    """
    Generate shRNA constructs (90mer or 97mer) for all rows where Final_valid is True.

    Parameters
    ----------
    df               : pd.DataFrame – output of analyze_sequences(), must have Final_valid column.
    clean_fasta      : str          – cleaned FASTA string (DNA).
    project_name     : str          – optional project name for shRNA naming.
    construct_format : str          – "90mer" or "97mer" (default: "97mer").

    Returns
    -------
    pd.DataFrame with one row per valid candidate, containing all construct components.
    """
    if construct_format not in FORMAT_COMPONENTS:
        raise ValueError(f"Unknown construct_format '{construct_format}'. Choose '90mer' or '97mer'.")

    five_prime_end, three_prime_end = FORMAT_COMPONENTS[construct_format]

    valid_df = df[df["Final_valid"] == True].copy().reset_index(drop=True)

    if valid_df.empty:
        return pd.DataFrame()

    # --- Sort by Off-Targets (asc) then Corrected_Score (desc) ---
    valid_df = valid_df.sort_values(
        by=["Off-Targets", "Corrected_Score"],
        ascending=[True, False]
    ).reset_index(drop=True)

    records = []

    for rank_0, (_, row) in enumerate(valid_df.iterrows()):
        rank = rank_0 + 1
        warnings_list = []

        # Build shRNA name
        prefix = f"shRNA_{project_name}_" if project_name.strip() else "shRNA_"
        shrna_name = f"{prefix}{rank}"

        # Antisense: U → T
        as_raw = str(row["AS Sequence"]).upper().strip()
        antisense_seq = as_raw.replace("U", "T")

        # Sense: reverse complement of antisense
        sense_seq = reverse_complement(antisense_seq)

        # LinkerB
        linker_b, lb_warning = get_linker_b(sense_seq, clean_fasta)
        if lb_warning:
            warnings_list.append(lb_warning)

        # LinkerA
        linker_a, la_warning = get_linker_a(linker_b)
        if la_warning:
            warnings_list.append(la_warning)

        # Final construct (even if linkers are missing, construct what we can)
        final_construct = (
            five_prime_end
            + linker_a
            + sense_seq
            + LOOP
            + antisense_seq
            + linker_b
            + three_prime_end
        )

        # Carry over existing warnings from analysis step
        existing_warnings = str(row.get("Warnings", "")).strip()
        if existing_warnings:
            warnings_list.insert(0, existing_warnings)

        records.append({
            "shRNA_Name":          shrna_name,
            "Format":              construct_format,
            "Rank":                rank,
            "AS Sequence":         row["AS Sequence"],
            "Antisense_Sequence":  antisense_seq,
            "Sense_Sequence":      sense_seq,
            "LinkerA":             linker_a,
            "LinkerB":             linker_b,
            "5prime_end":          five_prime_end,
            "Loop":                LOOP,
            "3prime_end":          three_prime_end,
            "Final_Construct":     final_construct,
            "Off-Targets":         row.get("Off-Targets"),
            "Corrected_Score":     row.get("Corrected_Score"),
            "Position":            row.get("Position"),
            "Warnings":            "; ".join(warnings_list) if warnings_list else "",
        })

    return pd.DataFrame(records)


def add_proximity_warnings(df_90mers: pd.DataFrame) -> pd.DataFrame:
    """
    Add proximity warnings when two candidates are within 15 nt of each other.

    The warning is appended to the Warnings column of both affected rows.
    Ranking is NOT changed.

    Parameters
    ----------
    df_90mers : pd.DataFrame – output of build_90mers()

    Returns
    -------
    pd.DataFrame with updated Warnings column.
    """
    if df_90mers.empty:
        return df_90mers

    df = df_90mers.copy()

    try:
        positions = df["Position"].astype(float).tolist()
    except (ValueError, TypeError):
        return df  # cannot compute proximity if positions aren't numeric

    n = len(positions)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(positions[i] - positions[j]) <= 15:
                msg = (
                    f"Candidate positions are within 15 nucleotides "
                    f"(#{df.at[i, 'Rank']} and #{df.at[j, 'Rank']})."
                )
                for idx in (i, j):
                    existing = df.at[idx, "Warnings"]
                    if existing:
                        df.at[idx, "Warnings"] = existing + "; " + msg
                    else:
                        df.at[idx, "Warnings"] = msg

    return df
