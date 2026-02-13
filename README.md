# NCBI Sequence QC & Annotation Pipeline (Google Colab Free Course)

Created for **eng. Taha bilel chalbi**.

Welcome! This repository contains a **Python 3 Colab-ready mini pipeline** that teaches and executes an end-to-end bioinformatics workflow for FASTA sequences:

- Accept pasted or uploaded FASTA.
- Detect if sequence is **DNA / RNA / Protein**.
- Perform **IUPAC sanity checks** and quick QC.
- Compute:
  - Sequence complexity (entropy + low-complexity indicator)
  - GC content (DNA/RNA)
  - Molecular weight (Protein)
- Ask for your NCBI credentials (email + optional API key).
- Query NCBI using alignment-based matching (BLAST) and return annotations.
- Decide whether the sequence appears **already represented** or **likely novel/variant**.

---

## 1) Free Course Roadmap

Think of this as a compact, practical course in 7 lessons.

## Lesson 1 — Setup in Google Colab

Open a new Colab notebook and install dependencies:

```python
!pip install biopython requests
```

Then upload or clone this repository and run:

```python
!python colab_ncbi_sequence_pipeline.py
```

---

## Lesson 2 — FASTA Input Basics

The script supports two modes:

1. **Paste FASTA** directly in terminal input
2. **Upload file** (Colab mode)

Example FASTA:

```fasta
>example_seq
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

When pasting manually, finish with:

```text
END
```

---

## Lesson 3 — Sequence Type Detection

The pipeline inspects symbols and classifies each record as:

- `DNA`
- `RNA`
- `PROTEIN`
- `NUCLEIC_ACID_AMBIGUOUS`
- `UNKNOWN`

It uses IUPAC-compatible character sets (including ambiguous symbols) to make this decision.

---

## Lesson 4 — QC and Sanity Checks

For each input sequence, the pipeline reports:

- Length
- Invalid characters (if any)
- Ambiguity fraction
- Shannon entropy
- Complexity index (0..1)
- Low-complexity warning

Extra calculations by type:

- **DNA/RNA**: GC content (%)
- **Protein**: molecular weight (Da)

---

## Lesson 5 — NCBI Access Requirements

Before NCBI lookup, the script asks for:

- Your **email** (required by NCBI Entrez usage policy)
- Your **API key** (optional but recommended for higher request limits)

> Tip: create an NCBI account and API key for better throughput.

---

## Lesson 6 — Database Search with Alignment Evidence

Because users may provide only a raw sequence (without accession ID), the pipeline uses **BLAST alignment-based evidence**:

- `blastn` against `nt` for DNA/RNA
- `blastp` against `nr` for proteins

It reports top hits with:

- Accession
- Title
- Identity %
- Query coverage %
- E-value
- Alignment snippets

Then it pulls NCBI summary/GenBank details (name, organism, authors when available).

---

## Lesson 7 — Novelty Interpretation

Current rule of thumb in the script:

- Likely known if best hit has high identity and high coverage.
- Candidate novel/variant if no near-exact full-length match.
- Novel if no significant alignment hits are returned.

> Important: “novel” in this script is a computational indicator, not a publication-grade novelty claim.

---

## 2) What This Pipeline Teaches You

- Practical FASTA handling
- Sequence QC fundamentals
- NCBI Entrez + BLAST integration
- Annotation retrieval workflow
- Reproducible analysis structure for Colab

---

## 3) How to Run (Quick Start)

```bash
python colab_ncbi_sequence_pipeline.py
```

Flow:

1. Choose input mode.
2. Provide FASTA.
3. Review QC.
4. Choose NCBI lookup (`y`/`n`).
5. Enter email (+ optional API key).
6. Review BLAST evidence and annotations.

---

## 4) Suggested Practice Exercises (Course Style)

1. Run one clean DNA FASTA and inspect GC%.
2. Run RNA with ambiguous symbols and inspect ambiguity fraction.
3. Run a protein sequence and verify molecular weight output.
4. Compare one known sequence vs one synthetic/random sequence to observe novelty interpretation differences.

---

## 5) Notes & Limitations

- BLAST/Entrez depend on NCBI network availability.
- Runtime depends on query size and NCBI server load.
- Classification and novelty decisions are heuristic and should be reviewed by a domain expert for critical applications.

---

## 6) File in This Branch

- `colab_ncbi_sequence_pipeline.py` — main interactive pipeline script
- `README.md` — this free-course style guide

---

## 7) Acknowledgment

Prepared with your requested framing and naming: **eng. Taha bilel chalbi**.
