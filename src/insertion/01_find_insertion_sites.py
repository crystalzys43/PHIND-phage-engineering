"""
Phase 2, Step 1: Predict safe luciferase insertion sites in top phage candidates.

WHY THIS STEP EXISTS
--------------------
Phase 1 told us WHICH phages to engineer. Phase 2 tells us WHERE in each
phage genome to drop the luciferase (luxAB / NanoLuc) reporter cassette.

The wrong choice destroys the phage's lytic cycle and wastes 3+ months of
wet-lab time. The right choice gives synchronous, strong luminescence
during phage replication — exactly what PHIND needs.

BIOLOGICAL CRITERIA FOR A SAFE INSERTION SITE
---------------------------------------------
A. INTERGENIC: must be located in a non-coding region between two CDS,
   never inside a CDS or in a promoter/RBS zone (we use a ≥50 bp gap as
   the working definition of "intergenic with room").

B. NOT ESSENTIAL FLANKED: the flanking genes should not be both
   essential (capsid, portal, terminase, etc.). One essential + one
   non-essential is acceptable; two non-essentials is best.

C. LATE-GENE PROXIMITY: lysis cluster proximity is strongly preferred,
   because (1) the lysis genes share a late promoter that times
   luciferase peak to host-cell lysis, and (2) Loessner's seminal 1996
   A511-luxAB construct was placed near this region — established
   biological precedent.

D. SIZE: gap must accommodate the cassette + a small buffer. Standard
   luxAB cassette with promoter = ~2.5 kb. NanoLuc cassette = ~1 kb.
   We score gaps ≥1000 bp as "comfortable", ≥500 bp as "tight",
   ≥50 bp as "minimal."

GENE CATEGORIES (keyword-based, transparent)
--------------------------------------------
ESSENTIAL (don't disrupt):
  - structural:  capsid, portal, tail, baseplate, head, scaffold
  - packaging:   terminase, packaging
  - replication: polymerase, helicase, primase, ligase, DNA repair

LYSIS (preferred flanks):
  - endolysin, holin, lysin, spanin, hydrolase, peptidase, amidase

PERMISSIVE (good to flank with):
  - hypothetical, putative, unknown — usually accessory/dispensable
  - HNH endonuclease, homing endonuclease — often dispensable

Usage
-----
    python src/insertion/01_find_insertion_sites.py

Inputs
------
    data/raw/*.gb                          (GenBank files from Phase 1 Step 1)
    results/candidate_phages_ranked.csv    (from Phase 1 Step 3)

Outputs
-------
    results/insertion_sites_top10.csv      — ranked sites for top 10 phages
    results/04_insertion_log.txt           — per-phage site summary
"""

from __future__ import annotations

import csv
import re
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
CSV_OUT = PROJECT_ROOT / "results" / "insertion_sites_top10.csv"
LOG_OUT = PROJECT_ROOT / "results" / "04_insertion_log.txt"

# How many top phages to analyze (start with top 10; can re-run with more)
TOP_N = 10

# Minimum intergenic gap to qualify as a candidate insertion site
MIN_GAP_BP = 50

# Cassette size thresholds (used for scoring gap size)
COMFORTABLE_BP = 1000
TIGHT_BP = 500

# Gene category keywords (case-insensitive regex)
ESSENTIAL_KEYWORDS = [
    r"\bcapsid\b", r"\bportal\b", r"\btail\b", r"\bbaseplate\b",
    r"\bhead\b", r"\bscaffold\b", r"\bterminase\b", r"\bpackaging\b",
    r"\bpolymerase\b", r"\bhelicase\b", r"\bprimase\b", r"\bligase\b",
    r"\bsingle.?strand DNA binding\b", r"\bDNA\s+repair\b",
]
LYSIS_KEYWORDS = [
    r"\bendolysin\b", r"\bholin\b", r"\blysin\b", r"\bspanin\b",
    r"\bhydrolase\b", r"\bpeptidase\b", r"\bamidase\b",
    r"\bN-acetylmuramoyl\b", r"\bcell\s+wall\b",
]
PERMISSIVE_KEYWORDS = [
    r"\bhypothetical\b", r"\bputative\b", r"\bunknown\b",
    r"\bHNH endonuclease\b", r"\bhoming endonuclease\b",
]

ESSENTIAL_RE = re.compile("|".join(ESSENTIAL_KEYWORDS), re.IGNORECASE)
LYSIS_RE = re.compile("|".join(LYSIS_KEYWORDS), re.IGNORECASE)
PERMISSIVE_RE = re.compile("|".join(PERMISSIVE_KEYWORDS), re.IGNORECASE)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def categorize(product: str) -> str:
    """Bucket a CDS product name into one of: essential, lysis, permissive, other."""
    if not product:
        return "other"
    if LYSIS_RE.search(product):
        # Check lysis first — lysin/hydrolase win over essential-ish keywords
        return "lysis"
    if ESSENTIAL_RE.search(product):
        return "essential"
    if PERMISSIVE_RE.search(product):
        return "permissive"
    return "other"


def find_gb_path(accession: str) -> Path | None:
    """The .gb files are named by NCBI UID, not accession. Find by parsing."""
    for f in RAW_DIR.glob("*.gb"):
        try:
            rec = SeqIO.read(f, "genbank")
            if rec.id == accession:
                return f
        except Exception:  # noqa: BLE001
            continue
    return None


def extract_cds_table(record) -> pd.DataFrame:
    """Build a table of CDS with positions and categorized function."""
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        try:
            start = int(feat.location.start)
            end = int(feat.location.end)
        except Exception:  # noqa: BLE001
            continue
        product = feat.qualifiers.get("product", ["unknown"])[0]
        rows.append(
            {
                "start": start,
                "end": end,
                "strand": feat.location.strand,
                "product": product,
                "category": categorize(product),
            }
        )
    df = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
    return df


def score_site(
    gap_bp: int, flank_left_cat: str, flank_right_cat: str,
    distance_to_lysis_bp: int, genome_size: int,
) -> tuple[float, list[str]]:
    """
    Score an intergenic insertion site (0–100) with reasoning.
    Returns (score, list_of_reasons).
    """
    score = 0.0
    reasons = []

    # 1. Gap size (max 25)
    if gap_bp >= COMFORTABLE_BP:
        score += 25; reasons.append(f"comfortable gap ({gap_bp} bp ≥ 1000)")
    elif gap_bp >= TIGHT_BP:
        score += 15; reasons.append(f"tight gap ({gap_bp} bp)")
    else:
        score += 5; reasons.append(f"minimal gap ({gap_bp} bp)")

    # 2. Flank essentiality (max 30)
    cats = (flank_left_cat, flank_right_cat)
    if "essential" not in cats:
        score += 30; reasons.append("neither flank essential")
    elif cats.count("essential") == 1:
        score += 15; reasons.append("one flank essential (acceptable)")
    else:
        score += 0; reasons.append("⚠ both flanks essential")

    # 3. Lysis proximity (max 30) — closer = better
    if "lysis" in cats:
        score += 30; reasons.append("⭐ adjacent to lysis gene")
    elif distance_to_lysis_bp < 2000:
        score += 20; reasons.append(f"near lysis cluster ({distance_to_lysis_bp} bp)")
    elif distance_to_lysis_bp < 5000:
        score += 10; reasons.append(f"within 5 kb of lysis cluster")
    else:
        score += 0; reasons.append("far from lysis cluster")

    # 4. Permissive flank bonus (max 15)
    n_permissive = cats.count("permissive")
    if n_permissive == 2:
        score += 15; reasons.append("both flanks permissive (hypothetical)")
    elif n_permissive == 1:
        score += 8; reasons.append("one permissive flank")

    return round(score, 1), reasons


def find_lysis_positions(cds: pd.DataFrame) -> list[int]:
    """Return list of midpoints of lysis-category CDS."""
    lysis_rows = cds[cds["category"] == "lysis"]
    return [int((r["start"] + r["end"]) / 2) for _, r in lysis_rows.iterrows()]


def analyze_phage(accession: str, name: str) -> list[dict]:
    """For one phage, return a list of scored insertion-site dicts."""
    gb_path = find_gb_path(accession)
    if gb_path is None:
        return []

    record = SeqIO.read(gb_path, "genbank")
    cds = extract_cds_table(record)
    if len(cds) < 2:
        return []

    lysis_positions = find_lysis_positions(cds)

    sites = []
    for i in range(len(cds) - 1):
        left = cds.iloc[i]
        right = cds.iloc[i + 1]
        gap = int(right["start"]) - int(left["end"])
        if gap < MIN_GAP_BP:
            continue

        gap_midpoint = int(left["end"]) + gap // 2
        if lysis_positions:
            dist_to_lysis = min(abs(gap_midpoint - lp) for lp in lysis_positions)
        else:
            dist_to_lysis = 10_000_000  # essentially infinite if no lysis genes

        score, reasons = score_site(
            gap_bp=gap,
            flank_left_cat=left["category"],
            flank_right_cat=right["category"],
            distance_to_lysis_bp=dist_to_lysis,
            genome_size=len(record.seq),
        )

        sites.append(
            {
                "accession": accession,
                "name": name,
                "site_start": int(left["end"]),
                "site_end": int(right["start"]),
                "gap_bp": gap,
                "left_gene": left["product"][:60],
                "left_category": left["category"],
                "right_gene": right["product"][:60],
                "right_category": right["category"],
                "dist_to_lysis_bp": dist_to_lysis if dist_to_lysis < 1e6 else None,
                "insertion_score": score,
                "rationale": " · ".join(reasons),
            }
        )

    sites.sort(key=lambda s: -s["insertion_score"])
    return sites


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    ranked = pd.read_csv(RANKED_CSV).head(TOP_N)
    print(f"\n[1/3] Analyzing top {len(ranked)} phages from the ranked list.\n")

    all_sites: list[dict] = []
    for _, row in tqdm(ranked.iterrows(), total=len(ranked), desc="phages"):
        sites = analyze_phage(row["accession"], row["name"])
        all_sites.extend(sites)

    out_df = pd.DataFrame(all_sites)
    out_df.to_csv(CSV_OUT, index=False)
    print(f"\n[2/3] Wrote {len(out_df)} candidate insertion sites to {CSV_OUT}")

    # ---- Per-phage summary ----
    log_lines = [
        "Phase 2 Step 1: Insertion site analysis",
        "=========================================",
        f"Phages analyzed: {len(ranked)}",
        f"Total candidate sites found: {len(out_df)}",
        "",
        "Per-phage top-3 sites:",
    ]
    for _, row in ranked.iterrows():
        acc = row["accession"]
        phage_sites = out_df[out_df["accession"] == acc].head(3)
        log_lines.append("")
        log_lines.append(f"▸ {acc} — {row['name'][:60]} (rank #{int(row['rank'])})")
        if phage_sites.empty:
            log_lines.append("    no qualifying sites found")
            continue
        for _, s in phage_sites.iterrows():
            log_lines.append(
                f"    score={s['insertion_score']:5.1f}  "
                f"pos {s['site_start']:>7,}..{s['site_end']:<7,}  "
                f"gap={s['gap_bp']:>5} bp"
            )
            log_lines.append(f"        left  ({s['left_category']:10s}): {s['left_gene']}")
            log_lines.append(f"        right ({s['right_category']:10s}): {s['right_gene']}")
            log_lines.append(f"        rationale: {s['rationale']}")

    print("\n[3/3] Per-phage summary:\n")
    for line in log_lines:
        print("    " + line)
    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
