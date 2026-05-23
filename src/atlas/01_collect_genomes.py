"""
Phase 1, Step 1: Collect publicly available Listeria phage genomes from NCBI.

This script queries NCBI's Nucleotide database for complete genomes of phages
that infect Listeria species, downloads them as GenBank files, and writes a
summary CSV of basic metadata (accession, length, host, isolation source).

Usage
-----
    python src/atlas/01_collect_genomes.py

Outputs
-------
    data/raw/*.gb              — one GenBank file per phage genome
    data/listeria_phages.csv   — metadata summary table
    results/01_collection_log.txt — run log

Notes
-----
- Uses anonymous NCBI Entrez access. NCBI requests up to 3 queries/sec without
  an API key; we throttle accordingly.
- Re-running the script skips genomes already downloaded (idempotent).
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Tell NCBI who is querying (required courtesy). Replace with your own email.
Entrez.email = "yushanz5@illinois.edu"

# Search query: complete phage genomes infecting Listeria, in the size range
# typical of tailed phages (15kb – 250kb). Excludes prophage fragments and
# very short partial sequences.
SEARCH_QUERY = (
    '("Listeria"[Organism] OR Listeria[All Fields]) '
    'AND ("complete genome"[Title] OR "complete sequence"[Title]) '
    'AND (phage[Title] OR bacteriophage[Title]) '
    'AND 15000:250000[SLEN]'
)

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
RESULTS_DIR = PROJECT_ROOT / "results"
CSV_OUT = DATA_DIR / "listeria_phages.csv"
LOG_OUT = RESULTS_DIR / "01_collection_log.txt"

# Throttle between NCBI requests (seconds). 0.34s = ~3 req/sec, NCBI safe.
REQUEST_DELAY = 0.34


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def search_ncbi(query: str) -> list[str]:
    """Return a list of NCBI nucleotide accession IDs matching the query."""
    print(f"\n[1/3] Searching NCBI Nucleotide:\n      {query}\n")

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=2000)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    print(f"      → Found {len(id_list)} candidate records.")
    time.sleep(REQUEST_DELAY)
    return id_list


def fetch_genome(uid: str, out_path: Path) -> dict | None:
    """Download one GenBank record. Returns a metadata dict, or None on error."""
    if out_path.exists():
        # Already downloaded — just parse for metadata
        record = SeqIO.read(out_path, "genbank")
    else:
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=uid, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()
            SeqIO.write([record], out_path, "genbank")
        except Exception as exc:  # noqa: BLE001
            print(f"      ! Failed to fetch UID {uid}: {exc}")
            return None
        time.sleep(REQUEST_DELAY)

    # Extract metadata from the GenBank record
    src_feature = next(
        (f for f in record.features if f.type == "source"), None
    )
    host = ""
    isolation_source = ""
    country = ""
    if src_feature:
        host = "; ".join(src_feature.qualifiers.get("host", []))
        isolation_source = "; ".join(
            src_feature.qualifiers.get("isolation_source", [])
        )
        country = "; ".join(src_feature.qualifiers.get("country", []))

    return {
        "accession": record.id,
        "name": record.description,
        "length_bp": len(record.seq),
        "gc_percent": round(
            100 * (record.seq.count("G") + record.seq.count("C")) / len(record.seq),
            2,
        ),
        "host": host,
        "isolation_source": isolation_source,
        "country": country,
        "organism": record.annotations.get("organism", ""),
        "taxonomy": "; ".join(record.annotations.get("taxonomy", [])),
    }


def main() -> None:
    # Set up output folders
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: search
    uids = search_ncbi(SEARCH_QUERY)
    if not uids:
        print("No records found. Check the search query.")
        return

    # Step 2: fetch each genome (skip already-downloaded ones)
    print(f"\n[2/3] Downloading GenBank files to {RAW_DIR}/\n")
    rows: list[dict] = []
    for uid in tqdm(uids, desc="phages"):
        out_path = RAW_DIR / f"{uid}.gb"
        meta = fetch_genome(uid, out_path)
        if meta is not None:
            rows.append(meta)

    # Step 3: write summary CSV
    print(f"\n[3/3] Writing metadata table to {CSV_OUT}\n")
    fieldnames = [
        "accession",
        "name",
        "length_bp",
        "gc_percent",
        "host",
        "isolation_source",
        "country",
        "organism",
        "taxonomy",
    ]
    with CSV_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    # Write a small log file
    with LOG_OUT.open("w") as fh:
        fh.write(f"NCBI query: {SEARCH_QUERY}\n")
        fh.write(f"Records found: {len(uids)}\n")
        fh.write(f"Records downloaded successfully: {len(rows)}\n")
        fh.write(f"Genomes stored under: {RAW_DIR}\n")
        fh.write(f"Metadata CSV: {CSV_OUT}\n")

    print(f"✅ Done. {len(rows)} Listeria phage genomes collected.")
    print(f"   CSV: {CSV_OUT}")
    print(f"   Log: {LOG_OUT}")


if __name__ == "__main__":
    main()
