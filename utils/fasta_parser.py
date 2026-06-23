"""
fasta_parser.py
---------------
Handles parsing, validation, and cleaning of FASTA sequences.
"""

import re


VALID_BASES = set("ATUCGN")


def parse_fasta(text: str) -> dict:
    """
    Parse and validate a FASTA-format string.

    Parameters
    ----------
    text : str
        Raw FASTA text (either from file upload or text box).

    Returns
    -------
    dict with keys:
        - "header"     : str  – the '>' header line (without '>')
        - "Clean_FASTA": str  – cleaned, uppercase, DNA sequence (U→T)
        - "errors"     : list – validation error messages
        - "warnings"   : list – non-fatal warning messages
    """
    errors = []
    warnings = []
    header = ""
    clean_seq = ""

    if not text or not text.strip():
        errors.append("FASTA input is empty.")
        return {"header": header, "Clean_FASTA": clean_seq, "errors": errors, "warnings": warnings}

    lines = text.strip().splitlines()

    # --- Validate header line ---
    if not lines[0].startswith(">"):
        errors.append("FASTA format error: first line must start with '>'.")
        return {"header": header, "Clean_FASTA": clean_seq, "errors": errors, "warnings": warnings}

    header = lines[0][1:].strip()

    # --- Extract sequence lines ---
    seq_lines = lines[1:]
    if not seq_lines or all(l.strip() == "" for l in seq_lines):
        errors.append("FASTA format error: no sequence found after header.")
        return {"header": header, "Clean_FASTA": clean_seq, "errors": errors, "warnings": warnings}

    # --- Clean sequence: remove whitespace, join, uppercase ---
    raw_seq = "".join(l.strip() for l in seq_lines).upper()
    raw_seq = raw_seq.replace(" ", "").replace("\t", "")

    # --- Validate characters ---
    invalid_chars = set(raw_seq) - VALID_BASES
    if invalid_chars:
        errors.append(
            f"FASTA sequence contains invalid characters: {', '.join(sorted(invalid_chars))}. "
            "Only A, T, U, C, G, N are allowed."
        )
        return {"header": header, "Clean_FASTA": clean_seq, "errors": errors, "warnings": warnings}

    # --- Convert RNA → DNA (U → T) ---
    if "U" in raw_seq:
        warnings.append("RNA sequence detected (U found). Converting U → T for downstream analysis.")
        raw_seq = raw_seq.replace("U", "T")

    clean_seq = raw_seq

    return {
        "header": header,
        "Clean_FASTA": clean_seq,
        "errors": errors,
        "warnings": warnings,
    }
