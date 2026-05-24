"""
Phase 2, Step 2: Re-annotate top candidates with Pharokka.

WHY THIS STEP EXISTS
--------------------
A511 and P100 — the two phages our Phase 1 ranking identified as top
reporter-engineering backbones — have minimal NCBI annotation. Every CDS
is labeled "gp1", "gp2", etc., with no functional descriptions. As a
result, Phase 2 Step 1's keyword classifier could not categorize their
flanking genes, capping their insertion-site scores at 55/100. That's a
data-annotation artifact, not a biological deficiency.

Pharokka resolves this by re-predicting all CDS with PHANOTATE and
annotating against the PHROG (Phage Family Orthologous Groups) database —
the standard functional annotation for phage genomes. Each CDS gets
assigned to one of ~3,500 PHROG families with a functional category
(e.g., "lysis", "tail", "DNA, RNA and nucleotide metabolism").

This script:
  1. Reads the top N candidates from the Phase 1 ranking
  2. Extracts each genome to FASTA
  3. Runs Pharokka against the PHROG database
  4. Collects Pharokka's per-CDS functional categories
  5. Writes a unified per-CDS table for Phase 2.3 re-scoring

Prerequisites
-------------
- Pharokka installed in conda env `pharokka_env`
- PHROG database downloaded to ~/pharokka_db (via install_databases.py)

Usage
-----
    python src/insertion/02_pharokka_annotate.py

Inputs
------
    results/candidate_phages_ranked.csv  (top candidates from Phase 1 Step 3)
    data/raw/*.gb                        (downloaded GenBank files)

Outputs
-------
    results/pharokka/<accession>/        (Pharokka output per candidate)
    results/pharokka_annotations.csv     (unified CDS table with PHROG categories)
    results/05_pharokka_log.txt          (run summary)
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
RANKED_CSV = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"
PHAROKKA_OUT_ROOT = PROJECT_ROOT / "results" / "pharokka"
ANNOTATIONS_CSV = PROJECT_ROOT / "results" / "pharokka_annotations.csv"
LOG_OUT = PROJECT_ROOT / "results" / "05_pharokka_log.txt"

# Pharokka environment and database
PHAROKKA_BIN = Path.home() / "miniconda3/envs/pharokka_env/bin/pharokka.py"
PHAROKKA_DB = Path.home() / "pharokka_db"

# Number of top candidates to re-annotate (start small to validate the pipeline)
TOP_N = 5

# How many CPU threads Pharokka can use
THREADS = 4


# ---------------------------------------------------------------------------
# Pipeline functions
# ---------------------------------------------------------------------------

def find_gb_path(accession: str) -> Path | None:
    """GenBank files are stored by NCBI UID, not accession — find by parsing."""
    for f in RAW_DIR.glob("*.gb"):
        try:
            rec = SeqIO.read(f, "genbank")
            if rec.id == accession:
                return f
        except Exception:  # noqa: BLE001
            continue
    return None


def gb_to_fasta(gb_path: Path, fasta_out: Path) -> None:
    """Convert a GenBank file to FASTA (Pharokka takes FASTA input)."""
    record = SeqIO.read(gb_path, "genbank")
    # Make ID safe (no spaces, no colons) for Pharokka
    record.id = record.id.replace(".", "_").replace(":", "_")
    record.description = ""
    SeqIO.write([record], fasta_out, "fasta")


def run_pharokka(fasta: Path, out_dir: Path) -> int:
    """
    Invoke pharokka on a single FASTA. We use `conda run` to activate the
    pharokka_env so pharokka's subprocess calls (phanotate, mmseqs, etc.)
    find their binaries in PATH.

    IMPORTANT: `conda run` re-parses arguments through a shell, so paths
    containing spaces break it. As a workaround we stage the input FASTA
    and output directory in `/tmp` (no spaces guaranteed), run pharokka
    there, then copy the results back to the project's results/ folder.
    """
    import shutil
    import tempfile

    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="pharokka_") as tmp:
        tmp_path = Path(tmp)
        tmp_fasta = tmp_path / fasta.name
        tmp_out = tmp_path / "out"
        shutil.copy(fasta, tmp_fasta)

        cmd = [
            str(Path.home() / "miniconda3/bin/conda"),
            "run", "-n", "pharokka_env", "--no-capture-output",
            "pharokka.py",
            "-i", str(tmp_fasta),
            "-o", str(tmp_out),
            "-d", str(PHAROKKA_DB),
            "-t", str(THREADS),
            "--fast",          # mmseqs2-only mode (sufficient for our purpose)
            "-f",              # force overwrite
            "-p", fasta.stem,  # output prefix
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"\n  stderr (last 600 chars): {result.stderr[-600:]}")
            return result.returncode

        # Copy outputs back
        for f in tmp_out.iterdir():
            target = out_dir / f.name
            if f.is_dir():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(f, target)
            else:
                shutil.copy(f, target)
    return 0


def parse_pharokka_output(out_dir: Path, accession: str, prefix: str) -> list[dict]:
    """
    Parse Pharokka's `<prefix>_cds_final_merged_output.tsv` into a list of dicts.
    Each row = one CDS with its PHROG annotation.
    """
    tsv = out_dir / f"{prefix}_cds_final_merged_output.tsv"
    if not tsv.exists():
        return []

    df = pd.read_csv(tsv, sep="\t")
    rows = []
    for _, r in df.iterrows():
        rows.append(
            {
                "accession": accession,
                "cds_id": r.get("ID", ""),
                "start": int(r.get("start", 0)),
                "end": int(r.get("stop", 0)),
                "strand": r.get("frame", ""),
                "phrog_id": r.get("phrog", ""),
                "phrog_category": r.get("category", ""),  # PHROG functional category
                "phrog_annot": r.get("annot", ""),         # PHROG product name
                "method": r.get("Method", ""),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # Sanity checks
    if not PHAROKKA_BIN.exists():
        sys.exit(f"❌ Pharokka not found at {PHAROKKA_BIN}. "
                 f"Install with: CONDA_SUBDIR=osx-64 conda create -n pharokka_env "
                 f"-c bioconda -c conda-forge -y pharokka python=3.10")
    if not PHAROKKA_DB.exists():
        sys.exit(f"❌ Pharokka database not found at {PHAROKKA_DB}. "
                 f"Download with: install_databases.py -o {PHAROKKA_DB}")

    PHAROKKA_OUT_ROOT.mkdir(parents=True, exist_ok=True)

    ranked = pd.read_csv(RANKED_CSV).head(TOP_N)
    print(f"\n[1/3] Re-annotating top {len(ranked)} candidates with Pharokka.\n")

    all_annotations: list[dict] = []
    log_lines = ["Phase 2 Step 2: Pharokka re-annotation summary",
                 "============================================="]

    for _, row in tqdm(ranked.iterrows(), total=len(ranked), desc="phages"):
        acc = row["accession"]
        safe_id = acc.replace(".", "_")
        gb_path = find_gb_path(acc)
        if gb_path is None:
            log_lines.append(f"  ! {acc}: GenBank file not found, skipping")
            continue

        # Step a: convert to FASTA
        out_dir = PHAROKKA_OUT_ROOT / safe_id
        out_dir.mkdir(parents=True, exist_ok=True)
        fasta = out_dir / f"{safe_id}.fasta"
        gb_to_fasta(gb_path, fasta)

        # Step b: run Pharokka
        rc = run_pharokka(fasta, out_dir)
        if rc != 0:
            log_lines.append(f"  ! {acc}: Pharokka exited with code {rc}")
            continue

        # Step c: parse output
        rows = parse_pharokka_output(out_dir, acc, safe_id)
        all_annotations.extend(rows)

        # Summary stats for this phage
        if rows:
            df_p = pd.DataFrame(rows)
            cat_counts = df_p["phrog_category"].value_counts().to_dict()
            log_lines.append(f"  ▸ {acc} ({row['name'][:50]})")
            log_lines.append(f"     CDS predicted: {len(rows)}")
            log_lines.append(f"     PHROG categories:")
            for cat, n in sorted(cat_counts.items(), key=lambda x: -x[1])[:8]:
                log_lines.append(f"       {cat:50s}  {n}")
        else:
            log_lines.append(f"  ! {acc}: no CDS parsed from Pharokka output")

    # Write unified table
    with ANNOTATIONS_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "accession", "cds_id", "start", "end", "strand",
                "phrog_id", "phrog_category", "phrog_annot", "method",
            ],
        )
        writer.writeheader()
        writer.writerows(all_annotations)

    print(f"\n[2/3] Wrote {len(all_annotations)} annotated CDS to {ANNOTATIONS_CSV}\n")
    print("[3/3] Per-phage summary:\n")
    for line in log_lines:
        print("    " + line)
    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
