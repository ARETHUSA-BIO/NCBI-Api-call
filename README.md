# NCBI Sequence QC & Annotation Pipeline

Prepared by Taha bilel chalbi**.

A short, Colab-ready pipeline to:
- read FASTA,
- detect DNA/RNA/protein,
- perform IUPAC/QC checks,
- compute sequence metrics,
- query NCBI with BLAST evidence,
- report annotations or likely novelty.

---

## Features

For each sequence record, the script returns:
- type and length,
- invalid characters,
- ambiguity fraction,
- Shannon entropy,
- complexity index + low-complexity flag,
- GC% (DNA/RNA),
- molecular weight (protein).

For NCBI search, it asks for:
- email (required),
- API key (optional).

Then it:
- runs BLAST (`blastn` for DNA/RNA, `blastp` for protein),
- shows top hits (identity, coverage, e-value),
- fetches summary/annotation details (authors when available).

---

## Requirements

```bash
pip install biopython requests
```

Main script:
- `colab_ncbi_sequence_pipeline.py`

---

## Run

Local:

```bash
python colab_ncbi_sequence_pipeline.py
```

Colab:

```python
!pip install biopython requests
!python colab_ncbi_sequence_pipeline.py
```

---

## FASTA input

Modes:
1. Paste FASTA manually (finish with `END`)
2. Upload FASTA file (Colab)

Example:

```fasta
>example_seq
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

---

## Presence / novelty interpretation

The script uses alignment evidence (not accession-only exact lookup):
- high identity + high coverage: likely already represented,
- weaker near-hit: candidate novel/variant,
- no significant hit: likely novel.

> This is computational guidance, not a formal novelty claim for publication.

---

## Notes

- Requires internet and NCBI availability.
- Runtime depends on sequence size and NCBI load.
- Expert review is recommended for critical decisions.

---

## Files

- `colab_ncbi_sequence_pipeline.py`
- `README.md`

---
