"""
Phase 1, Step 2: Filter to Listeria hosts and classify lytic vs temperate lifestyle.

WHY THIS STEP EXISTS
--------------------
PHIND's reporter phage must be **lytic** (virulent). Lytic phages always
kill the host cell after infection, producing enough luciferase signal to
detect. Temperate phages integrate their DNA into the bacterial chromosome
(lysogeny) and stay dormant — useless as a reporter.

We also discovered in Step 1 that the NCBI keyword search caught some
non-Listeria phages (Streptococcus, Bacillus, etc.) that happened to mention
Listeria somewhere in their metadata. We filter those out here.

CLASSIFICATION HEURISTIC
------------------------
A phage is flagged as **putatively temperate** if its annotated CDS list
contains any of the canonical lysogeny machinery genes:

    integrase, excisionase, recombinase, CI/CII-like repressor,
    antirepressor, transposase, attP/attB sites

If none of these are present, we call the phage **putatively lytic**.

This is the same gene-content signal used by BACPHLIP (the standard ML
classifier for phage lifestyle), applied transparently from existing GenBank
annotations — no external dependencies required, and every prediction can be
audited by reading the evidence list.

Usage
-----
    python src/atlas/02_classify_lifestyle.py

Inputs
------
    data/raw/*.gb        (downloaded in Step 1)

Outputs
-------
    data/listeria_phages_classified.csv   — Listeria-host phages with lifestyle
    results/02_classification_log.txt     — run summary
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
CSV_OUT = PROJECT_ROOT / "data" / "listeria_phages_classified.csv"
LOG_OUT = PROJECT_ROOT / "results" / "02_classification_log.txt"

# Keywords that indicate a phage is putatively temperate (lysogenic).
# Each is a regex pattern matched against CDS product names (case-insensitive).
TEMPERATE_KEYWORDS = [
    r"\bintegrase\b",
    r"\bexcisionase\b",
    r"\brecombinase\b",
    r"\bCI repressor\b",
    r"\bCII repressor\b",
    r"\bCro-like\b",
    r"\bantirepressor\b",
    r"\btranspos",          # transposase, transposon
    r"\battP\b",
    r"\battB\b",
    r"\blysogen",           # lysogeny module / lysogenic
]
TEMPERATE_RE = re.compile("|".join(TEMPERATE_KEYWORDS), re.IGNORECASE)


# ---------------------------------------------------------------------------
# Per-phage analyses
# ---------------------------------------------------------------------------

def is_listeria_phage(record) -> bool:
    """True if the phage's host (or fallback organism name) is Listeria."""
    src = next((f for f in record.features if f.type == "source"), None)
    if src:
        host_values = src.qualifiers.get("host", [])
        if host_values:
            return any("listeria" in h.lower() for h in host_values)
    # Fallback: many records have no host field but say "Listeria phage X"
    organism = record.annotations.get("organism", "").lower()
    return "listeria" in organism


def classify_lifestyle(record) -> tuple[str, list[str]]:
    """
    Inspect CDS product names for lysogeny machinery.
    Returns (lifestyle, list_of_evidence_genes).
    """
    products = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        products.extend(feat.qualifiers.get("product", []))

    evidence = [prod for prod in products if TEMPERATE_RE.search(prod)]
    lifestyle = "temperate" if evidence else "lytic"
    return lifestyle, evidence


def extract_metadata(record) -> dict:
    """Pull the same metadata fields as Step 1, plus a few we want here."""
    src = next((f for f in record.features if f.type == "source"), None)
    host = isolation_source = country = ""
    if src:
        host = "; ".join(src.qualifiers.get("host", []))
        isolation_source = "; ".join(src.qualifiers.get("isolation_source", []))
        country = "; ".join(src.qualifiers.get("country", []))

    seq = record.seq
    gc = round(100 * (seq.count("G") + seq.count("C")) / len(seq), 2)

    return {
        "accession": record.id,
        "name": record.description,
        "length_bp": len(seq),
        "gc_percent": gc,
        "host": host,
        "isolation_source": isolation_source,
        "country": country,
        "organism": record.annotations.get("organism", ""),
        "taxonomy": "; ".join(record.annotations.get("taxonomy", [])),
    }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    gb_files = sorted(RAW_DIR.glob("*.gb"))
    print(f"\n[1/3] Found {len(gb_files)} GenBank files in {RAW_DIR}\n")

    rows: list[dict] = []
    print("[2/3] Filtering Listeria hosts + classifying lifestyle...\n")

    for gb_path in tqdm(gb_files, desc="phages"):
        try:
            record = SeqIO.read(gb_path, "genbank")
        except Exception as exc:  # noqa: BLE001
            print(f"  ! Failed to parse {gb_path.name}: {exc}")
            continue

        meta = extract_metadata(record)
        meta["is_listeria_host"] = is_listeria_phage(record)
        meta["cds_count"] = sum(1 for f in record.features if f.type == "CDS")

        lifestyle, evidence = classify_lifestyle(record)
        meta["lifestyle"] = lifestyle
        meta["n_lysogeny_genes"] = len(evidence)
        meta["lifestyle_evidence"] = "; ".join(evidence[:5])  # top 5 hits

        rows.append(meta)

    # Write CSV (Listeria-host phages only, with full classification)
    listeria_rows = [r for r in rows if r["is_listeria_host"]]
    fieldnames = list(listeria_rows[0].keys()) if listeria_rows else []
    with CSV_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(listeria_rows)

    # ---- Stats ----
    n_total = len(rows)
    n_listeria = sum(1 for r in rows if r["is_listeria_host"])
    n_listeria_lytic = sum(
        1 for r in rows if r["is_listeria_host"] and r["lifestyle"] == "lytic"
    )
    n_listeria_temperate = sum(
        1 for r in rows if r["is_listeria_host"] and r["lifestyle"] == "temperate"
    )

    log_lines = [
        "Step 2 classification summary",
        "=============================",
        f"Total .gb files processed:               {n_total}",
        f"  → with Listeria host:                  {n_listeria}",
        f"      → lytic (REPORTER CANDIDATES):    {n_listeria_lytic}",
        f"      → temperate (excluded):           {n_listeria_temperate}",
        "",
        f"Output CSV: {CSV_OUT}",
    ]

    print("\n[3/3]")
    for line in log_lines:
        print("      " + line)

    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")

    # Preview top 5 lytic Listeria candidates by genome size
    listeria_lytic = [
        r for r in rows
        if r["is_listeria_host"] and r["lifestyle"] == "lytic"
    ]
    listeria_lytic.sort(key=lambda r: -r["length_bp"])
    print("\n      Top 5 lytic Listeria-host candidates (by genome size):")
    for r in listeria_lytic[:5]:
        print(f"        • {r['accession']:14s} {r['length_bp']:>7,} bp   {r['name'][:50]}")


if __name__ == "__main__":
    main()
