"""
Colab-ready Python 3 pipeline for sequence intake, QC, characterization, and NCBI lookup.

Workflow:
1) Accept FASTA input (paste or upload).
2) Detect whether the sequence is DNA, RNA, or protein.
3) Perform IUPAC sanity checks and QC summary.
4) Compute sequence metrics:
   - Complexity (Shannon entropy + low-complexity estimate)
   - GC content for DNA/RNA
   - Molecular weight for proteins
5) Query NCBI for exact-sequence matches and return annotations.

Usage in Colab:
    !python colab_ncbi_sequence_pipeline.py
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from getpass import getpass
from typing import Dict, List, Optional, Sequence, Tuple

# Runtime dependency notes:
# pip install biopython requests
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# IUPAC alphabets (including ambiguous symbols)
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")
IUPAC_PROTEIN = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
GAP_CHARS = set("-.")


@dataclass
class SequenceQC:
    record_id: str
    seq_type: str
    length: int
    invalid_characters: List[str]
    ambiguity_fraction: float
    gc_content: Optional[float]
    shannon_entropy_bits: float
    complexity_index: float
    low_complexity_warning: bool
    molecular_weight_da: Optional[float]


def normalize_sequence(raw: str) -> str:
    """Remove whitespace and FASTA header lines, return uppercase sequence string."""
    lines = [ln.strip() for ln in raw.strip().splitlines() if ln.strip()]
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    seq = "".join(seq_lines).upper().replace(" ", "")
    return seq


def parse_fasta(raw_fasta: str) -> List[Tuple[str, str]]:
    """Very lightweight FASTA parser that supports multi-record input."""
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []

    for line in raw_fasta.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, normalize_sequence("\n".join(seq_parts))))
            header = line[1:].strip() or "unnamed_record"
            seq_parts = []
        else:
            seq_parts.append(line)

    if header is not None:
        records.append((header, normalize_sequence("\n".join(seq_parts))))
    elif raw_fasta.strip():
        # fallback: no FASTA header provided; treat as a single sequence
        records.append(("unnamed_record", normalize_sequence(raw_fasta)))

    return records


def detect_sequence_type(seq: str) -> str:
    """Detect sequence nature: DNA, RNA, PROTEIN, or UNKNOWN."""
    cleaned = [c for c in seq if c not in GAP_CHARS]
    if not cleaned:
        return "UNKNOWN"

    chars = set(cleaned)
    dna_ok = chars.issubset(IUPAC_DNA)
    rna_ok = chars.issubset(IUPAC_RNA)
    prot_ok = chars.issubset(IUPAC_PROTEIN)

    has_u = "U" in chars
    has_t = "T" in chars

    if dna_ok and not has_u:
        return "DNA"
    if rna_ok and not has_t:
        return "RNA"
    if dna_ok and rna_ok:
        return "NUCLEIC_ACID_AMBIGUOUS"
    if prot_ok:
        return "PROTEIN"
    return "UNKNOWN"


def find_invalid_characters(seq: str, seq_type: str) -> List[str]:
    chars = set(c for c in seq if c not in GAP_CHARS)
    if seq_type in {"DNA", "NUCLEIC_ACID_AMBIGUOUS"}:
        valid = IUPAC_DNA
    elif seq_type == "RNA":
        valid = IUPAC_RNA
    elif seq_type == "PROTEIN":
        valid = IUPAC_PROTEIN
    else:
        valid = IUPAC_DNA | IUPAC_RNA | IUPAC_PROTEIN

    invalid = sorted(ch for ch in chars if ch not in valid)
    return invalid


def ambiguity_fraction(seq: str, seq_type: str) -> float:
    if not seq:
        return 0.0
    if seq_type in {"DNA", "RNA", "NUCLEIC_ACID_AMBIGUOUS"}:
        ambiguous = set("NRYSWKMBDHV")
    elif seq_type == "PROTEIN":
        ambiguous = set("BXZJUO*")
    else:
        ambiguous = set("NRYSWKMBDHVBXZJUO*")

    clean = [c for c in seq if c not in GAP_CHARS]
    if not clean:
        return 0.0
    amb_count = sum(1 for c in clean if c in ambiguous)
    return amb_count / len(clean)


def gc_content(seq: str) -> Optional[float]:
    clean = [c for c in seq if c not in GAP_CHARS]
    if not clean:
        return None
    nuc = [c for c in clean if c in set("ACGTUN")]
    if not nuc:
        return None
    gc = sum(1 for c in nuc if c in {"G", "C"})
    return 100.0 * gc / len(nuc)


def shannon_entropy(seq: str) -> float:
    clean = [c for c in seq if c not in GAP_CHARS]
    if not clean:
        return 0.0
    counts: Dict[str, int] = {}
    for c in clean:
        counts[c] = counts.get(c, 0) + 1
    entropy = 0.0
    n = len(clean)
    for k in counts.values():
        p = k / n
        entropy -= p * math.log2(p)
    return entropy


def complexity_index(seq: str, window: int = 12) -> float:
    """
    Heuristic complexity score in [0,1].
    Combines normalized entropy and repetitive-window uniqueness.
    """
    clean = "".join(c for c in seq if c not in GAP_CHARS)
    if not clean:
        return 0.0

    ent = shannon_entropy(clean)
    alphabet = len(set(clean))
    max_ent = math.log2(alphabet) if alphabet > 1 else 1.0
    norm_ent = min(1.0, ent / max_ent) if max_ent > 0 else 0.0

    if len(clean) < window:
        return norm_ent

    windows = [clean[i : i + window] for i in range(0, len(clean) - window + 1)]
    unique_ratio = len(set(windows)) / len(windows)
    return max(0.0, min(1.0, 0.6 * norm_ent + 0.4 * unique_ratio))


def molecular_weight(seq: str) -> Optional[float]:
    prot = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    if not prot:
        return None
    try:
        return ProteinAnalysis(prot).molecular_weight()
    except Exception:
        return None


def qc_sequence(record_id: str, seq: str) -> SequenceQC:
    seq_type = detect_sequence_type(seq)
    invalid = find_invalid_characters(seq, seq_type)
    amb_frac = ambiguity_fraction(seq, seq_type)
    gc = gc_content(seq) if seq_type in {"DNA", "RNA", "NUCLEIC_ACID_AMBIGUOUS"} else None
    ent = shannon_entropy(seq)
    cpx = complexity_index(seq)
    mw = molecular_weight(seq) if seq_type == "PROTEIN" else None

    return SequenceQC(
        record_id=record_id,
        seq_type=seq_type,
        length=len([c for c in seq if c not in GAP_CHARS]),
        invalid_characters=invalid,
        ambiguity_fraction=amb_frac,
        gc_content=gc,
        shannon_entropy_bits=ent,
        complexity_index=cpx,
        low_complexity_warning=cpx < 0.35,
        molecular_weight_da=mw,
    )


def configure_ncbi(email: str, api_key: str = "") -> None:
    Entrez.email = email.strip()
    if api_key.strip():
        Entrez.api_key = api_key.strip()


def search_exact_sequence_in_ncbi(seq: str, seq_type: str, retmax: int = 5) -> Dict:
    """
    Try exact sequence search in NCBI Entrez.
    Note: exact matching behavior can vary by sequence length and indexing.
    """
    if seq_type in {"DNA", "RNA", "NUCLEIC_ACID_AMBIGUOUS"}:
        db = "nuccore"
    elif seq_type == "PROTEIN":
        db = "protein"
    else:
        return {"status": "skipped", "reason": "Unknown sequence type"}

    term = f'"{seq}"[Sequence]'
    with Entrez.esearch(db=db, term=term, retmax=retmax) as handle:
        res = Entrez.read(handle)

    ids = res.get("IdList", [])
    if not ids:
        return {
            "status": "novel",
            "db": db,
            "message": "No exact indexed match found in NCBI for this query.",
            "ids": [],
        }

    with Entrez.esummary(db=db, id=",".join(ids), retmode="xml") as handle:
        summaries = Entrez.read(handle)

    annotations = []
    for uid in ids:
        doc = None
        for d in summaries:
            if str(d.get("Id", "")) == str(uid):
                doc = d
                break
        if doc is None:
            continue
        annotations.append(
            {
                "id": str(uid),
                "title": str(doc.get("Title", "")),
                "accession": str(doc.get("Caption", "")),
                "length": int(doc.get("Length", 0)) if str(doc.get("Length", "")).isdigit() else doc.get("Length", ""),
                "organism": str(doc.get("Organism", "")),
            }
        )

    return {
        "status": "found",
        "db": db,
        "ids": ids,
        "annotations": annotations,
    }


def blast_similarity_search(seq: str, seq_type: str, hitlist_size: int = 10) -> Dict:
    """
    Run NCBI BLAST similarity search and return alignment-derived hit metrics.

    Why this is needed:
    Users may submit only a raw sequence (no accession), so sequence-similarity
    alignment is required before deciding whether the sequence is already known.
    """
    if seq_type in {"DNA", "RNA", "NUCLEIC_ACID_AMBIGUOUS"}:
        program = "blastn"
        database = "nt"
    elif seq_type == "PROTEIN":
        program = "blastp"
        database = "nr"
    else:
        return {"status": "skipped", "reason": "Unknown sequence type"}

    handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=seq,
        hitlist_size=hitlist_size,
        format_type="XML",
    )
    blast_record = NCBIXML.read(handle)

    hits = []
    for aln in blast_record.alignments[:hitlist_size]:
        if not aln.hsps:
            continue
        hsp = aln.hsps[0]
        identity_pct = (100.0 * hsp.identities / hsp.align_length) if hsp.align_length else 0.0
        coverage_pct = (100.0 * hsp.align_length / max(1, len(seq)))
        hits.append(
            {
                "accession": aln.accession,
                "title": aln.title,
                "length": aln.length,
                "score": hsp.score,
                "evalue": hsp.expect,
                "identity_pct": round(identity_pct, 3),
                "query_coverage_pct": round(min(coverage_pct, 100.0), 3),
                "query_alignment": hsp.query,
                "match_alignment": hsp.match,
                "subject_alignment": hsp.sbjct,
            }
        )

    if not hits:
        return {
            "status": "novel",
            "program": program,
            "database": database,
            "message": "No significant BLAST hits found.",
            "hits": [],
        }

    best = hits[0]
    known = best["identity_pct"] >= 99.0 and best["query_coverage_pct"] >= 95.0
    return {
        "status": "found" if known else "candidate_novel",
        "program": program,
        "database": database,
        "decision": (
            "Sequence likely already represented in NCBI (high identity + coverage)."
            if known
            else "No near-exact full-length hit; sequence may be novel/variant."
        ),
        "best_hit": best,
        "hits": hits,
    }


def fetch_accession_summaries(accessions: Sequence[str], seq_type: str) -> List[Dict]:
    """Fetch Entrez summaries for BLAST hit accessions."""
    if not accessions:
        return []
    db = "protein" if seq_type == "PROTEIN" else "nuccore"
    with Entrez.esearch(db=db, term=" OR ".join(f"{acc}[Accession]" for acc in accessions), retmax=len(accessions)) as h:
        res = Entrez.read(h)
    ids = res.get("IdList", [])
    if not ids:
        return []
    with Entrez.esummary(db=db, id=",".join(ids), retmode="xml") as h:
        docs = Entrez.read(h)
    out = []
    for doc in docs:
        out.append(
            {
                "id": str(doc.get("Id", "")),
                "accession": str(doc.get("Caption", "")),
                "title": str(doc.get("Title", "")),
                "organism": str(doc.get("Organism", "")),
                "length": str(doc.get("Length", "")),
            }
        )
    return out


def fetch_genbank_details(db: str, ids: Sequence[str], limit: int = 3) -> List[Dict]:
    """Fetch richer annotations including authors from GenBank records."""
    out: List[Dict] = []
    use_ids = list(ids)[:limit]
    if not use_ids:
        return out

    with Entrez.efetch(db=db, id=",".join(use_ids), rettype="gb", retmode="text") as handle:
        records = list(SeqIO.parse(handle, "genbank"))

    for rec in records:
        refs = rec.annotations.get("references", [])
        authors = []
        if refs:
            for r in refs[:2]:
                if getattr(r, "authors", None):
                    authors.append(r.authors)
        out.append(
            {
                "id": rec.id,
                "name": rec.name,
                "description": rec.description,
                "organism": rec.annotations.get("organism", ""),
                "taxonomy": rec.annotations.get("taxonomy", []),
                "authors": authors,
                "keywords": rec.annotations.get("keywords", []),
            }
        )

    return out


def get_fasta_from_user() -> str:
    print("Choose FASTA input mode:")
    print("  [1] Paste FASTA manually")
    print("  [2] Upload file (Colab only)")
    mode = input("Enter 1 or 2 [default 1]: ").strip() or "1"

    if mode == "2":
        try:
            from google.colab import files  # type: ignore
        except Exception:
            print("Colab file-upload not available; falling back to paste mode.")
            mode = "1"

    if mode == "2":
        uploaded = files.upload()  # type: ignore[name-defined]
        if not uploaded:
            raise RuntimeError("No file uploaded.")
        first_name = next(iter(uploaded.keys()))
        return uploaded[first_name].decode("utf-8")

    print("Paste FASTA content. End with a line containing only: END")
    lines: List[str] = []
    while True:
        line = input()
        if line.strip() == "END":
            break
        lines.append(line)
    return "\n".join(lines)


def run_pipeline() -> None:
    print("=== NCBI Sequence QC & Annotation Pipeline (Colab-ready) ===")
    raw_fasta = get_fasta_from_user()
    records = parse_fasta(raw_fasta)

    if not records:
        raise ValueError("No sequence records detected.")

    print(f"Detected {len(records)} FASTA record(s).")

    qc_results: List[SequenceQC] = []
    for rid, seq in records:
        qc = qc_sequence(rid, seq)
        qc_results.append(qc)

    print("\n=== QC SUMMARY ===")
    for qc in qc_results:
        print(f"\nRecord: {qc.record_id}")
        print(f"  Type: {qc.seq_type}")
        print(f"  Length: {qc.length}")
        print(f"  Invalid chars: {qc.invalid_characters if qc.invalid_characters else 'None'}")
        print(f"  Ambiguity fraction: {qc.ambiguity_fraction:.4f}")
        print(f"  Shannon entropy (bits): {qc.shannon_entropy_bits:.4f}")
        print(f"  Complexity index [0-1]: {qc.complexity_index:.4f}")
        print(f"  Low complexity warning: {qc.low_complexity_warning}")
        if qc.gc_content is not None:
            print(f"  GC content (%): {qc.gc_content:.2f}")
        if qc.molecular_weight_da is not None:
            print(f"  Molecular weight (Da): {qc.molecular_weight_da:.2f}")

    use_ncbi = input("\nRun NCBI lookup now? [y/N]: ").strip().lower() == "y"
    if not use_ncbi:
        print("NCBI lookup skipped.")
        return

    email = input("Enter your email for NCBI Entrez: ").strip()
    api_key = getpass("Enter your NCBI API key (optional, hidden input): ")
    configure_ncbi(email=email, api_key=api_key)

    print("\n=== NCBI LOOKUP ===")
    for (rid, seq), qc in zip(records, qc_results):
        print(f"\nRecord: {rid}")
        print("Running BLAST alignment search (MSA-style evidence from top hits)...")
        blast_result = blast_similarity_search(seq=seq, seq_type=qc.seq_type, hitlist_size=5)
        print(json.dumps(blast_result, indent=2))

        if blast_result.get("hits"):
            top_accessions = [h["accession"] for h in blast_result["hits"][:3] if h.get("accession")]
            summaries = fetch_accession_summaries(top_accessions, qc.seq_type)
            print("NCBI summaries for top alignment hits:")
            print(json.dumps(summaries, indent=2))

            db = "protein" if qc.seq_type == "PROTEIN" else "nuccore"
            ids = [s["id"] for s in summaries if s.get("id")]
            details = fetch_genbank_details(db=db, ids=ids, limit=3)
            print("Detailed annotations (subset):")
            print(json.dumps(details, indent=2, default=str))
        else:
            print("Interpretation: no significant alignment hit; sequence is likely novel.")


if __name__ == "__main__":
    run_pipeline()
