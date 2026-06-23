"""
app.py
------
shRNA 90mer / 97mer Generator – Streamlit main application.

Run with:
    streamlit run app.py
"""

import io
import math
import pandas as pd
import streamlit as st

from utils.fasta_parser import parse_fasta
from utils.validators import validate_csv_columns
from utils.sequence_filters import analyze_sequences
from utils.sequence_builder import build_constructs, add_proximity_warnings
from utils.exporters import (
    df_to_csv_bytes,
    df_to_xlsx_bytes,
    analyses_export_columns,
    constructs_export_columns,
)

# ──────────────────────────────────────────────
# Page configuration
# ──────────────────────────────────────────────
st.set_page_config(
    page_title="shRNA Construct Generator",
    page_icon="🧬",
    layout="wide",
)

# ──────────────────────────────────────────────
# Color palette for construct segments
# ──────────────────────────────────────────────
SEGMENT_COLORS = {
    "5prime_end":         "#4A90D9",   # blue
    "LinkerA":            "#E67E22",   # orange
    "Sense_Sequence":     "#27AE60",   # green
    "Loop":               "#8E44AD",   # purple
    "Antisense_Sequence": "#E74C3C",   # red
    "LinkerB":            "#F39C12",   # amber
    "3prime_end":         "#16A085",   # teal
}

SEGMENT_LABELS = {
    "5prime_end":         "5′ End",
    "LinkerA":            "LinkerA",
    "Sense_Sequence":     "Sense Sequence",
    "Loop":               "Loop",
    "Antisense_Sequence": "Antisense Sequence",
    "LinkerB":            "LinkerB",
    "3prime_end":         "3′ End",
}

# ──────────────────────────────────────────────
# Tooltip text constants
# ──────────────────────────────────────────────
HOW_TO_USE = """**How to Use This Application**

1. Enter a project name (optional).
2. Upload or paste the target gene FASTA sequence.
3. Upload the CSV file containing candidate shRNA sequences generated on \
[DSIR](http://biodev.extra.cea.fr/DSIR/DSIR.html).
4. Select the desired construct format (90mer or 97mer).
5. Review the automatically generated sequence analysis table, which displays \
all filtering criteria and candidate rankings.
6. Inspect the generated shRNA constructs and associated warnings.
7. Export the analysis table and final shRNA sequences in CSV or Excel format.
"""

FORMAT_INFO = """**90mer and 97mer Formats**

Two construct formats are available:

- **97mer**: original miR30-based shRNA design described by Dow et al.
- **90mer**: shortened version containing reduced 5′ and 3′ flanking sequences.

The 90mer format has been experimentally validated in our laboratory and provides \
comparable functionality in our cloning workflow. Both formats are therefore provided \
to allow users to choose the construct design that best fits their experimental requirements.
"""

# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────

def color_bool(val):
    """Pandas Styler function: green for True, red for False."""
    if val is True or val == "TRUE":
        return "background-color: #059669; color: #ffffff; font-weight: 600;"
    if val is False or val == "FALSE":
        return "background-color: #dc2626; color: #ffffff; font-weight: 600;"
    return ""


def color_warnings(val):
    """Pandas Styler: amber background for non-empty warning strings."""
    if isinstance(val, str) and val.strip():
        return "background-color: #d97706; color: #ffffff; font-weight: 600;"
    return ""


def render_colored_construct(row: pd.Series):
    """Render one construct row as HTML with per-segment color spans."""
    segments = [
        ("5prime_end",         row.get("5prime_end", "")),
        ("LinkerA",            row.get("LinkerA", "")),
        ("Sense_Sequence",     row.get("Sense_Sequence", "")),
        ("Loop",               row.get("Loop", "")),
        ("Antisense_Sequence", row.get("Antisense_Sequence", "")),
        ("LinkerB",            row.get("LinkerB", "")),
        ("3prime_end",         row.get("3prime_end", "")),
    ]
    html_parts = []
    for key, seq in segments:
        color = SEGMENT_COLORS.get(key, "#000")
        html_parts.append(
            f'<span style="color:{color}; font-weight:bold; '
            f'font-family:monospace; font-size:1rem;">{seq}</span>'
        )
    return "".join(html_parts)


def render_legend():
    """Render the shared color legend as an HTML row."""
    parts = []
    for key, label in SEGMENT_LABELS.items():
        color = SEGMENT_COLORS[key]
        parts.append(
            f'<span style="display:inline-block; margin-right:12px;">'
            f'<span style="background:{color}; color:white; padding:2px 8px; '
            f'border-radius:4px; font-size:0.8rem;">{label}</span></span>'
        )
    st.markdown("**Legend:** " + "".join(parts), unsafe_allow_html=True)


# ──────────────────────────────────────────────
# CSV reader — handles metadata headers & semicolon delimiters
# ──────────────────────────────────────────────
def _read_candidate_csv(csv_file) -> pd.DataFrame:
    """
    Read the candidate CSV produced by third-party shRNA tools.

    Handles:
    - Semicolon-delimited files (common export format)
    - Metadata/comment rows above the actual column header
    - Quoted column names (e.g. "AS Sequence")
    - Standard comma-delimited files

    The function scans lines from the top until it finds a row that
    contains all required column names, then uses that row as the header.
    """
    from utils.validators import REQUIRED_CSV_COLUMNS

    raw = csv_file.read().decode("utf-8", errors="replace")
    csv_file.seek(0)

    lines = raw.splitlines()

    # Auto-detect delimiter: count semicolons vs commas in first 10 non-empty lines
    sample = [l for l in lines if l.strip()][:10]
    n_semi  = sum(l.count(";") for l in sample)
    n_comma = sum(l.count(",") for l in sample)
    sep = ";" if n_semi > n_comma else ","

    # Find the header row: first line whose unquoted tokens contain all required columns
    header_line_idx = None
    for i, line in enumerate(lines):
        tokens = {t.strip().strip('"').strip("'") for t in line.split(sep)}
        if all(col in tokens for col in REQUIRED_CSV_COLUMNS):
            header_line_idx = i
            break

    if header_line_idx is None:
        return pd.read_csv(io.StringIO(raw), sep=sep)

    data_text = "\n".join(lines[header_line_idx:])
    df = pd.read_csv(io.StringIO(data_text), sep=sep)
    df.columns = [c.strip().strip('"').strip("'") for c in df.columns]
    return df


# ──────────────────────────────────────────────
# Session state initialisation
# ──────────────────────────────────────────────
def init_state():
    defaults = {
        "clean_fasta":      None,
        "fasta_header":     None,
        "fasta_warnings":   [],
        "df_analyses":      None,
        "df_constructs":    None,
        "construct_format": "97mer",
        "top_n_set":        set(),
        "other_set":        set(),
        "n_top":            5,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


init_state()

# ──────────────────────────────────────────────
# Sidebar – Inputs
# ──────────────────────────────────────────────
with st.sidebar:
    # ── Header row: title + info icon ────────
    title_col, info_col = st.columns([5, 1])
    with title_col:
        st.title("🧬 shRNA Generator")
    with info_col:
        st.markdown("<div style='padding-top:28px;'></div>", unsafe_allow_html=True)
        with st.popover("ℹ️"):
            st.markdown(HOW_TO_USE)

    st.markdown("---")

    # Project name
    project_name = st.text_input(
        "Project Name (optional)",
        placeholder="e.g. ABDC1",
        help="Used to name generated shRNAs: shRNA_ABDC1_1, shRNA_ABDC1_2, …",
    )

    st.markdown("---")

    # ── FASTA input ──────────────────────────
    st.subheader("1. FASTA Sequence")
    fasta_mode = st.radio("Input method", ["Upload file", "Paste sequence"], horizontal=True)

    fasta_text = None
    if fasta_mode == "Upload file":
        fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
        if fasta_file:
            fasta_text = fasta_file.read().decode("utf-8", errors="replace")
    else:
        fasta_text = st.text_area(
            "Paste FASTA sequence",
            height=150,
            placeholder=">NM_001234\nATGCATGCATGC...",
        )

    # ── CSV input ────────────────────────────
    st.markdown("---")
    st.subheader("2. shRNA Candidates CSV")
    csv_file = st.file_uploader("Upload CSV file", type=["csv"])

    # ── Construct format ─────────────────────
    st.markdown("---")
    fmt_label_col, fmt_info_col = st.columns([5, 1])
    with fmt_label_col:
        st.subheader("3. Construct Format")
    with fmt_info_col:
        st.markdown("<div style='padding-top:28px;'></div>", unsafe_allow_html=True)
        with st.popover("ℹ️"):
            st.markdown(FORMAT_INFO)

    construct_format = st.radio(
        "Select format",
        ["97mer", "90mer"],
        index=0,
        horizontal=True,
        label_visibility="collapsed",
    )

    st.markdown("---")

    run_btn = st.button("▶ Run Analysis", type="primary", use_container_width=True)

# ──────────────────────────────────────────────
# Main area – Introduction
# ──────────────────────────────────────────────
st.title("shRNA Construct Generator")

st.markdown("""
This application is designed to facilitate the generation of miR30-based shRNA constructs
for gene silencing experiments. Starting from a target transcript sequence and a list of
candidate antisense sequences, the software automatically applies established shRNA design
criteria, identifies high-confidence candidates, and generates ready-to-order shRNA 90mer
or 97mer constructs suitable for downstream PCR cloning workflows.

The filtering strategy implemented in this application is based on the Sensor-derived design
criteria described by Dow et al. and was developed to enrich for highly potent shRNAs while
minimizing ineffective candidates. These criteria were originally established as part of a
pipeline for the generation of inducible miR30-based shRNA transgenic mouse models.

**References:**

Dow LE, Premsrirut PK, Zuber J, Fellmann C, McJunkin K, Miething C, Park Y, Dickins RA,
Hannon GJ, Lowe SW. A pipeline for the generation of shRNA transgenic mice.
*Nature Protocols.* 2012;7(2):374–393. doi:10.1038/nprot.2011.446.

---
*Disclaimer: The ranking and filtering criteria implemented in this software are intended to
prioritize shRNA candidates with a higher predicted probability of efficacy. Final shRNA
performance should always be experimentally validated in the relevant biological system.*
""")

st.markdown("---")

# ──────────────────────────────────────────────
# Processing
# ──────────────────────────────────────────────
if run_btn:
    # Reset state
    for k in ("clean_fasta", "fasta_header", "fasta_warnings", "df_analyses", "df_constructs"):
        st.session_state[k] = None
    st.session_state["fasta_warnings"] = []
    st.session_state["top_n_set"] = set()
    st.session_state["other_set"] = set()
    st.session_state["construct_format"] = construct_format

    all_ok = True

    # ── Parse FASTA ──────────────────────────
    if not fasta_text or not fasta_text.strip():
        st.error("Please provide a FASTA sequence (upload or paste).")
        all_ok = False
    else:
        fasta_result = parse_fasta(fasta_text)
        if fasta_result["errors"]:
            for e in fasta_result["errors"]:
                st.error(f"FASTA Error: {e}")
            all_ok = False
        else:
            st.session_state["clean_fasta"]    = fasta_result["Clean_FASTA"]
            st.session_state["fasta_header"]   = fasta_result["header"]
            st.session_state["fasta_warnings"] = fasta_result["warnings"]

    # ── Parse CSV ────────────────────────────
    if csv_file is None:
        st.error("Please upload a CSV file.")
        all_ok = False
    else:
        try:
            df_csv = _read_candidate_csv(csv_file)
        except Exception as exc:
            st.error(f"Could not read CSV: {exc}")
            all_ok = False
            df_csv = None

        if df_csv is not None:
            col_errors = validate_csv_columns(df_csv)
            if col_errors:
                for e in col_errors:
                    st.error(f"CSV Error: {e}")
                all_ok = False

    if all_ok:
        fmt = st.session_state["construct_format"]

        with st.spinner("Analysing sequences…"):
            df_analyses = analyze_sequences(df_csv)
            st.session_state["df_analyses"] = df_analyses

        with st.spinner(f"Generating {fmt} constructs…"):
            df_constructs = build_constructs(
                df_analyses,
                st.session_state["clean_fasta"],
                project_name=project_name,
                construct_format=fmt,
            )
            if not df_constructs.empty:
                df_constructs = add_proximity_warnings(df_constructs)
            st.session_state["df_constructs"] = df_constructs

        # Initialise Top N set
        n = st.session_state["n_top"]
        if not df_constructs.empty:
            st.session_state["top_n_set"] = set(df_constructs.index[:n])
            st.session_state["other_set"] = set(df_constructs.index[n:])

        st.success("Analysis complete!")


# ──────────────────────────────────────────────
# Display FASTA info
# ──────────────────────────────────────────────
if st.session_state["clean_fasta"]:
    with st.expander("FASTA Info", expanded=False):
        st.markdown(f"**Header:** `{st.session_state['fasta_header']}`")
        seq = st.session_state["clean_fasta"]
        st.markdown(f"**Length:** {len(seq):,} nt")
        st.code(seq[:200] + ("…" if len(seq) > 200 else ""), language=None)
        for w in st.session_state["fasta_warnings"]:
            st.warning(w)


# ──────────────────────────────────────────────
# Analyses Sequences Table
# ──────────────────────────────────────────────
if st.session_state["df_analyses"] is not None:
    st.header("Analyses Sequences")

    df_a = st.session_state["df_analyses"]

    # Search / filter
    search_term = st.text_input("Search table (filters all text columns)", "")
    if search_term:
        mask = df_a.apply(
            lambda col: col.astype(str).str.contains(search_term, case=False, na=False)
        ).any(axis=1)
        df_display = df_a[mask]
    else:
        df_display = df_a

    # Styled display
    bool_cols = [c for c in df_display.columns if df_display[c].dtype == bool]
    styler = df_display.style
    for col in bool_cols:
        styler = styler.map(color_bool, subset=[col])
    if "Warnings" in df_display.columns:
        styler = styler.map(color_warnings, subset=["Warnings"])

    st.dataframe(styler, use_container_width=True, height=400)

    # Export buttons
    col1, col2 = st.columns(2)
    export_df = analyses_export_columns(df_a)
    with col1:
        st.download_button(
            "⬇ Download CSV",
            data=df_to_csv_bytes(export_df),
            file_name="Analyses_Sequences.csv",
            mime="text/csv",
        )
    with col2:
        st.download_button(
            "⬇ Download Excel",
            data=df_to_xlsx_bytes(export_df, sheet_name="Analyses_Sequences"),
            file_name="Analyses_Sequences.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )


# ──────────────────────────────────────────────
# Construct Results
# ──────────────────────────────────────────────
if st.session_state["df_constructs"] is not None:
    df_c = st.session_state["df_constructs"]
    fmt  = st.session_state["construct_format"]

    st.markdown("---")
    st.header(f"Generated {fmt} Constructs")

    if df_c.empty:
        st.warning("No valid candidates passed all filters. No constructs generated.")
    else:
        render_legend()
        st.markdown("---")

        # ── Controls ─────────────────────────
        col_n, col_all = st.columns([3, 1])
        with col_n:
            n_top = st.number_input(
                "Top N to display",
                min_value=1,
                max_value=len(df_c),
                value=min(st.session_state["n_top"], len(df_c)),
                step=1,
            )
            if n_top != st.session_state["n_top"]:
                st.session_state["n_top"] = n_top
                st.session_state["top_n_set"] = set(df_c.index[:n_top])
                st.session_state["other_set"]  = set(df_c.index[n_top:])

        with col_all:
            show_all = st.checkbox(f"Show All {fmt}s", value=False)

        # ── Top N section ─────────────────────
        st.subheader(f"Top {n_top} {fmt} Constructs")

        def render_construct_card(idx, row, section="top"):
            """Render a single construct as an expander card."""
            name    = row.get("shRNA_Name", f"shRNA_{idx+1}")
            rank    = row.get("Rank", "—")
            score   = row.get("Corrected_Score", "—")
            off_tgt = row.get("Off-Targets", "—")
            pos     = row.get("Position", "—")
            warn    = row.get("Warnings", "")
            final   = row.get("Final_Construct", "")
            cfmt    = row.get("Format", fmt)

            title = f"#{rank}  {name}  |  Score: {score}  |  Off-Targets: {off_tgt}  |  Position: {pos}"

            with st.expander(title, expanded=(rank == 1)):
                # Colored sequence
                colored_html = render_colored_construct(row)
                st.markdown(
                    f'<div style="word-break:break-all; font-family:monospace; '
                    f'font-size:1rem; padding:8px; background:#f8f8f8; '
                    f'border-radius:6px;">{colored_html}</div>',
                    unsafe_allow_html=True,
                )

                # Plain sequence for copying
                st.text_area(
                    f"Full {cfmt} construct (copy from here)",
                    value=final,
                    height=80,
                    key=f"seq_{section}_{idx}",
                )

                if warn:
                    st.warning(f"⚠ {warn}")

                # Move button
                if section == "top":
                    if st.button("Move to Other ↓", key=f"move_down_{idx}"):
                        st.session_state["top_n_set"].discard(idx)
                        st.session_state["other_set"].add(idx)
                        st.rerun()
                else:
                    if st.button("Move to Top N ↑", key=f"move_up_{idx}"):
                        st.session_state["other_set"].discard(idx)
                        st.session_state["top_n_set"].add(idx)
                        st.rerun()

        # Render top section
        current_top = sorted(st.session_state["top_n_set"])
        if current_top:
            for idx in current_top:
                if idx in df_c.index:
                    render_construct_card(idx, df_c.loc[idx], section="top")
        else:
            st.info("No candidates in Top N section.")

        # ── Other section ─────────────────────
        current_other = sorted(st.session_state["other_set"])
        if current_other or show_all:
            st.markdown("---")
            st.subheader(f"Other {fmt} Constructs")
            if current_other:
                for idx in current_other:
                    if idx in df_c.index:
                        render_construct_card(idx, df_c.loc[idx], section="other")
            else:
                st.info("No candidates in Other section.")

        # ── Export constructs ─────────────────
        st.markdown("---")
        st.subheader(f"Export {fmt} Constructs")
        export_c = constructs_export_columns(df_c)
        ecol1, ecol2 = st.columns(2)
        with ecol1:
            st.download_button(
                f"⬇ Download {fmt} CSV",
                data=df_to_csv_bytes(export_c),
                file_name=f"Generated_{fmt}_constructs.csv",
                mime="text/csv",
            )
        with ecol2:
            st.download_button(
                f"⬇ Download {fmt} Excel",
                data=df_to_xlsx_bytes(export_c, sheet_name=f"{fmt}_Constructs"),
                file_name=f"Generated_{fmt}_constructs.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )
