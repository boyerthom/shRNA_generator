# shRNA Construct Generator

A local web application that automates the generation of ready-to-order miR30-based shRNA constructs (90mer or 97mer) from a target FASTA sequence and a list of candidate antisense sequences.

The filtering strategy is based on the Sensor-derived design criteria described by **Dow et al. (2012)** and was developed to enrich for highly potent shRNAs while minimizing ineffective candidates.

> **Disclaimer:** The ranking and filtering criteria implemented in this software are intended to prioritize shRNA candidates with a higher predicted probability of efficacy. Final shRNA performance should always be experimentally validated in the relevant biological system.

---

## Features

- Upload or paste a target FASTA sequence (DNA or RNA)
- Upload a candidate CSV exported from [DSIR](http://biodev.extra.cea.fr/DSIR/DSIR.html)
- Automatic filtering based on established shRNA design rules
- Generation of 90mer or 97mer constructs with color-coded segment visualization
- Ranked candidate output with proximity warnings
- Export analysis table and final constructs to CSV or Excel

---

## Requirements

- Python 3.10+
- pip

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/shRNA_Generator.git
cd shRNA_Generator
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv venv

# Windows
venv\Scripts\activate

# macOS / Linux
source venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

---

## Running the Application

```bash
streamlit run app.py
```

The app will open automatically in your browser at `http://localhost:8501`.

---

## Usage

### Step 1 – Project Name (optional)
Enter a project name (e.g. `ABDC1`). Generated shRNAs will be named:
- `shRNA_ABDC1_1`, `shRNA_ABDC1_2`, …

If left blank: `shRNA_1`, `shRNA_2`, …

### Step 2 – FASTA Sequence
Upload a `.fa` / `.fasta` / `.txt` file **or** paste a FASTA sequence directly.

- First line must start with `>`
- Sequence may contain only: `A T U C G N` (case-insensitive)
- RNA sequences (`U`) are automatically converted to DNA (`T`)

### Step 3 – Candidates CSV
Upload the CSV exported from [DSIR](http://biodev.extra.cea.fr/DSIR/DSIR.html).

Required columns:

| Column | Description |
|---|---|
| `AS Sequence` | 21-nt antisense RNA sequence |
| `Off-Targets` | Number of off-target sites |
| `Corrected_Score` | Prediction score |
| `Position` | Position in the target sequence |

> The application automatically handles the semicolon-delimited format and metadata header rows produced by DSIR.

### Step 4 – Construct Format
Choose between **97mer** (default, original Dow et al. design) or **90mer** (shortened format).

| Format | 5′ End | 3′ End |
|---|---|---|
| 97mer | `TGCTGTTGACAGTGAGCG` | `TGCCTACTGCCTCGGA` |
| 90mer | `GTTGACAGTGAGCG` | `TGCCTACTGCCTC` |

### Step 5 – Run Analysis
Click **▶ Run Analysis**.

---

## Filtering Criteria

Each candidate is evaluated against the following rules:

| Column | Criterion |
|---|---|
| `Off-Targets_valid` | Off-Targets ≤ 1 |
| `First_base_AU_valid` | Position 1 is A, U, or T |
| `AU_content_valid` | 40–80% A+U/T content across full sequence |
| `AU_1_14_percentage_valid` | >50% A+U/T in positions 1–14 |
| `AU_ratio_valid` | (AU₁₋₁₄/14) / (AU₁₅₋₂₁/7) > 1 |
| `Position_13_14_valid` | Position 13 = A/U/T **or** position 14 = U/T |
| `No_A_position_20_valid` | Position 20 ≠ A |
| `No_forbidden_motifs_valid` | No AAAAAA, TTTTT, CCCC, or GGGG motifs |
| `Final_valid` | All criteria above are TRUE |

Candidates passing all criteria are used to generate constructs, ranked by Off-Targets (ascending) then Corrected_Score (descending).

---

## Construct Structure

```
5′ End  +  LinkerA  +  Sense  +  Loop  +  Antisense  +  LinkerB  +  3′ End
```

- **Loop:** `TAGTGAAGCCACAGATGTA`
- **LinkerB:** base immediately upstream of the sense sequence in the FASTA
- **LinkerA:** `C` if LinkerB ∈ {A, T, U} · `A` if LinkerB ∈ {C, G}

---

## Project Structure

```
shRNA_Generator/
│
├── app.py                      # Streamlit main application
├── requirements.txt
├── README.md
│
├── utils/
│   ├── __init__.py
│   ├── fasta_parser.py         # FASTA parsing and validation
│   ├── validators.py           # CSV and AS Sequence validation
│   ├── sequence_filters.py     # Per-row sequence analysis and filtering
│   ├── sequence_builder.py     # 90mer / 97mer construct generation
│   └── exporters.py            # CSV and XLSX export helpers
│
└── test_data/
    ├── test_fasta.fa            # Example FASTA file
    └── test_candidates.csv      # Example candidates CSV
```

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `ModuleNotFoundError` | Run `pip install -r requirements.txt` |
| FASTA validation error | Ensure sequence contains only A/T/U/C/G/N |
| CSV column missing | Ensure CSV has: `AS Sequence`, `Off-Targets`, `Corrected_Score`, `Position` |
| No constructs generated | No candidates passed all filters — check the Analyses table for which criteria failed |
| App does not open | Make sure the virtual environment is activated before running `streamlit run app.py` |

---

## Reference

Dow LE, Premsrirut PK, Zuber J, Fellmann C, McJunkin K, Miething C, Park Y, Dickins RA, Hannon GJ, Lowe SW.
**A pipeline for the generation of shRNA transgenic mice.**
*Nature Protocols.* 2012;7(2):374–393.
doi:[10.1038/nprot.2011.446](https://doi.org/10.1038/nprot.2011.446)
