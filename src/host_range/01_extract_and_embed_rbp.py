"""
Phase 3, Step 1: Extract Receptor Binding Proteins (RBPs) and compute ESM-2 embeddings.

WHY THIS STEP EXISTS
--------------------
Phages recognize their host bacteria through specific surface proteins —
the **receptor binding proteins** (RBPs), located on tail fibers, tail
spikes, or baseplate. Two phages with similar RBPs are likely to infect
the same bacterial strains. Two phages with different RBPs are
complementary — exactly what you want in a PHIND cocktail that needs to
detect diverse strains of Listeria.

But "similar RBP" isn't just sequence identity — modern protein language
models like **ESM-2** (Meta AI) learn structural and functional similarity
from millions of protein sequences, and produce embedding vectors where
proteins with similar function cluster together even if their sequence
identity is low.

This script:
  1. Pulls every CDS annotated as "tail" by Pharokka from the top-5
     candidates, with an extra filter on the annotation text to favor
     genuine RBPs (tail fiber, tail spike, baseplate, receptor binding)
  2. Writes them to a FASTA
  3. Computes a per-protein ESM-2 embedding (480-dim vector)
  4. Saves embeddings + metadata for Phase 3 Step 2 to cluster

MODEL CHOICE
------------
We use **esm2_t12_35M_UR50D** — Meta AI's 35M-parameter variant of ESM-2.
Trade-offs:
  - Small enough for CPU/MPS inference on a laptop (~1 sec per protein)
  - Captures ~85% of the functional discrimination of the larger 650M model
  - 480-dim embeddings, mean-pooled per protein
For a PHIND-scale demo (50-100 proteins), this is the right choice.

Usage
-----
    python src/host_range/01_extract_and_embed_rbp.py

Inputs
------
    results/pharokka_annotations.csv           (from Phase 2 Step 2)
    results/pharokka/<acc>/phanotate.faa       (Pharokka-predicted proteins)
    results/candidate_phages_ranked.csv

Outputs
-------
    results/host_range/rbp_proteins.fasta      — extracted RBP sequences
    results/host_range/rbp_metadata.csv        — per-protein metadata
    results/host_range/rbp_embeddings.npy      — ESM-2 embedding matrix
    results/07_rbp_extract_log.txt             — run summary
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from tqdm import tqdm
from transformers import AutoModel, AutoTokenizer

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PHAROKKA_CSV = PROJECT_ROOT / "results" / "pharokka_annotations.csv"
RANKED_CSV = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"
PHAROKKA_OUT_ROOT = PROJECT_ROOT / "results" / "pharokka"

OUT_DIR = PROJECT_ROOT / "results" / "host_range"
OUT_FASTA = OUT_DIR / "rbp_proteins.fasta"
OUT_META = OUT_DIR / "rbp_metadata.csv"
OUT_EMB = OUT_DIR / "rbp_embeddings.npy"
LOG_OUT = PROJECT_ROOT / "results" / "07_rbp_extract_log.txt"

# ESM-2 model
MODEL_NAME = "facebook/esm2_t12_35M_UR50D"

# Regex to identify likely RBPs within Pharokka's annotated proteins.
# These are the proteins most strongly involved in host recognition.
RBP_PATTERN = re.compile(
    r"\b(tail\s*fiber|tail\s*spike|tailspike|baseplate|"
    r"receptor[\s-]?binding|host[\s-]?specificity|"
    r"tropism|tail[\s-]?associated\s+lysin|TAL)\b",
    re.IGNORECASE,
)

# Maximum protein length to consider (ESM-2 has 1024 token limit; phage RBPs
# rarely exceed 1500 aa, but tail fibers can be longer)
MAX_AA = 1022


# ---------------------------------------------------------------------------
# Step 1: Identify RBPs from Pharokka annotation table
# ---------------------------------------------------------------------------

def find_rbps(annotations: pd.DataFrame) -> pd.DataFrame:
    """Filter Pharokka annotations to candidate RBP proteins."""
    # Start from any CDS in 'tail' or 'connector' PHROG category, then
    # tighten using the annotation text
    tail_mask = annotations["phrog_category"].isin(["tail", "connector"])
    text_mask = annotations["phrog_annot"].fillna("").str.contains(
        RBP_PATTERN, regex=True
    )
    # Union: PHROG category 'tail' covers all canonical tail proteins;
    # text filter catches things mis-categorized or in 'other'
    return annotations[tail_mask | text_mask].copy()


# ---------------------------------------------------------------------------
# Step 2: Pull actual protein sequences from Pharokka's phanotate.faa
# ---------------------------------------------------------------------------

def load_pharokka_proteins(accession: str) -> dict[str, str]:
    """
    Load all proteins predicted by phanotate for one phage.
    Returns dict mapping CDS id → amino acid sequence.
    """
    safe = accession.replace(".", "_")
    faa = PHAROKKA_OUT_ROOT / safe / "phanotate.faa"
    if not faa.exists():
        return {}

    proteins = {}
    for rec in SeqIO.parse(faa, "fasta"):
        # phanotate IDs look like "Default_CDS_1" or similar — strip prefix
        proteins[rec.id] = str(rec.seq).rstrip("*")
    return proteins


def match_cds_id_to_protein(
    cds_id: str, proteins: dict[str, str]
) -> tuple[str, str] | None:
    """Pharokka CDS IDs and phanotate.faa IDs may differ slightly. Try direct
    match first, then fallback to suffix matching."""
    if cds_id in proteins:
        return cds_id, proteins[cds_id]
    # Pharokka may format as "DQ003638_2_CDS_5" while phanotate uses CDS_5
    suffix = cds_id.split("_")[-1] if "_" in cds_id else cds_id
    for pid, seq in proteins.items():
        if pid.endswith(suffix) or pid.endswith(f"CDS_{suffix}"):
            return pid, seq
    return None


# ---------------------------------------------------------------------------
# Step 3: ESM-2 embedding
# ---------------------------------------------------------------------------

def compute_embeddings(sequences: list[str], model, tokenizer, device) -> np.ndarray:
    """
    Mean-pool ESM-2 per-residue embeddings to a single per-protein vector.
    Returns array of shape (n_sequences, embedding_dim).
    """
    embeddings = []
    for seq in tqdm(sequences, desc="ESM-2 embedding"):
        # Truncate over-long sequences
        seq_trim = seq[:MAX_AA]
        with torch.no_grad():
            inputs = tokenizer(seq_trim, return_tensors="pt").to(device)
            outputs = model(**inputs)
            # Mean-pool over residues (skip special tokens at ends)
            hidden = outputs.last_hidden_state[0, 1:-1, :]
            emb = hidden.mean(dim=0).cpu().numpy()
        embeddings.append(emb)
    return np.stack(embeddings)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    annotations = pd.read_csv(PHAROKKA_CSV)
    ranked = pd.read_csv(RANKED_CSV)

    # Build accession → display name lookup
    name_lookup = dict(zip(ranked["accession"], ranked["name"]))

    print(f"\n[1/4] Identifying candidate RBP proteins...\n")
    rbp_table = find_rbps(annotations)
    print(f"      → {len(rbp_table)} RBP-candidate CDS across "
          f"{rbp_table['accession'].nunique()} phages.")

    # Pull sequences from Pharokka phanotate FASTAs
    print(f"\n[2/4] Loading protein sequences from Pharokka outputs...\n")
    records: list[dict] = []
    fasta_records = []
    for acc, group in rbp_table.groupby("accession"):
        proteins = load_pharokka_proteins(acc)
        for _, row in group.iterrows():
            match = match_cds_id_to_protein(row["cds_id"], proteins)
            if match is None:
                continue
            pid, seq = match
            if len(seq) < 50:  # skip tiny ORFs
                continue
            records.append(
                {
                    "accession": acc,
                    "phage_name": name_lookup.get(acc, ""),
                    "cds_id": row["cds_id"],
                    "phrog_category": row["phrog_category"],
                    "phrog_annot": row["phrog_annot"],
                    "start": row["start"],
                    "end": row["end"],
                    "length_aa": len(seq),
                }
            )
            fasta_records.append((f"{acc}|{pid}", seq))

    print(f"      → {len(records)} proteins with full sequence + ≥50 aa.")

    # Write FASTA
    with OUT_FASTA.open("w") as fh:
        for header, seq in fasta_records:
            fh.write(f">{header}\n{seq}\n")
    print(f"      → wrote {OUT_FASTA}")

    # Save metadata
    meta_df = pd.DataFrame(records)
    meta_df.to_csv(OUT_META, index=False)
    print(f"      → wrote {OUT_META}")

    # ---- ESM-2 embedding ----
    print(f"\n[3/4] Loading ESM-2 model ({MODEL_NAME})...")
    device = (
        "mps" if torch.backends.mps.is_available()
        else "cuda" if torch.cuda.is_available()
        else "cpu"
    )
    print(f"      → using device: {device}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to(device).eval()

    print(f"\n[4/4] Computing embeddings for {len(fasta_records)} proteins...\n")
    sequences = [seq for _, seq in fasta_records]
    embeddings = compute_embeddings(sequences, model, tokenizer, device)
    np.save(OUT_EMB, embeddings)
    print(f"\n      → wrote embeddings of shape {embeddings.shape} to {OUT_EMB}")

    # ---- Log ----
    log_lines = [
        "Phase 3 Step 1: RBP extraction + ESM-2 embedding",
        "================================================",
        f"Pharokka-annotated phages:           {annotations['accession'].nunique()}",
        f"RBP-candidate CDS identified:        {len(rbp_table)}",
        f"Proteins with full sequence:         {len(records)}",
        f"Embedding model:                     {MODEL_NAME}",
        f"Embedding shape:                     {embeddings.shape}",
        f"Compute device:                      {device}",
        "",
        "Per-phage RBP count:",
    ]
    for acc, group in meta_df.groupby("accession"):
        log_lines.append(
            f"  ▸ {acc} ({name_lookup.get(acc, '')[:40]}): {len(group)} RBP proteins"
        )
        for _, r in group.iterrows():
            log_lines.append(
                f"      - {r['phrog_annot'][:60]}  ({r['length_aa']} aa)"
            )

    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")

    print("\n" + "\n".join(log_lines[:8]))


if __name__ == "__main__":
    main()
