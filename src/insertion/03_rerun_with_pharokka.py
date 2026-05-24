"""
Phase 2, Step 3: Re-score insertion sites using Pharokka PHROG categories.

WHY THIS STEP EXISTS
--------------------
Step 2 (02_pharokka_annotate.py) re-annotated top candidates with functional
PHROG categories. Now we redo the insertion-site scoring using those
categories instead of NCBI's variable-quality CDS product names. This is
what fixes the A511 / P100 underscoring problem from Phase 2 Step 1.

PHROG CATEGORY MAPPING
----------------------
PHROG assigns each phage CDS to one of these functional categories:

  - head and packaging              → ESSENTIAL  (capsid, terminase, portal)
  - tail                            → ESSENTIAL  (tail fiber, tail tube, etc.)
  - connector                       → ESSENTIAL  (head-tail joining)
  - DNA, RNA and nucleotide metabolism → ESSENTIAL  (polymerase, helicase)
  - integration and excision        → LYSOGENY   (should not appear in lytic)
  - transcription regulation        → REGULATORY (acceptable flanks)
  - moron, auxiliary metabolic gene and host takeover → ACCESSORY (great flanks)
  - lysis                           → LYSIS     (PREFERRED flanks for late genes)
  - other                           → PERMISSIVE
  - unknown function                → PERMISSIVE (likely hypothetical/accessory)

This is the **gold-standard** functional categorization for phage genomes,
replacing the keyword-matching heuristic from Phase 2 Step 1.

Usage
-----
    python src/insertion/03_rerun_with_pharokka.py

Inputs
------
    results/pharokka_annotations.csv      (from Phase 2 Step 2)
    results/candidate_phages_ranked.csv   (from Phase 1 Step 3)

Outputs
-------
    results/insertion_sites_pharokka.csv  — re-scored sites for top candidates
    results/06_pharokka_insertion_log.txt — summary + before/after comparison
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PHAROKKA_CSV = PROJECT_ROOT / "results" / "pharokka_annotations.csv"
RANKED_CSV = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"
OLD_INSERTION_CSV = PROJECT_ROOT / "results" / "insertion_sites_top10.csv"
CSV_OUT = PROJECT_ROOT / "results" / "insertion_sites_pharokka.csv"
LOG_OUT = PROJECT_ROOT / "results" / "06_pharokka_insertion_log.txt"

# Scoring parameters (mirror Phase 2 Step 1 for an apples-to-apples comparison)
MIN_GAP_BP = 50
COMFORTABLE_BP = 1000
TIGHT_BP = 500

# PHROG category → our four-class scheme
PHROG_TO_CLASS = {
    "head and packaging": "essential",
    "tail": "essential",
    "connector": "essential",
    "DNA, RNA and nucleotide metabolism": "essential",
    "integration and excision": "essential",   # should be absent in lytic; flag as essential not to disrupt
    "transcription regulation": "permissive",
    "moron, auxiliary metabolic gene and host takeover": "permissive",
    "lysis": "lysis",
    "other": "permissive",
    "unknown function": "permissive",
}


def phrog_to_class(category: str) -> str:
    if pd.isna(category) or not category:
        return "permissive"
    return PHROG_TO_CLASS.get(category.strip(), "permissive")


def score_site(
    gap_bp: int, flank_left_cat: str, flank_right_cat: str,
    distance_to_lysis_bp: int,
) -> tuple[float, list[str]]:
    """Same scoring scheme as Phase 2 Step 1 — for direct comparability."""
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

    # 3. Lysis proximity (max 30)
    if "lysis" in cats:
        score += 30; reasons.append("⭐ adjacent to lysis gene")
    elif distance_to_lysis_bp < 2000:
        score += 20; reasons.append(f"near lysis cluster ({distance_to_lysis_bp} bp)")
    elif distance_to_lysis_bp < 5000:
        score += 10; reasons.append("within 5 kb of lysis cluster")
    else:
        score += 0; reasons.append("far from lysis cluster")

    # 4. Permissive flank bonus (max 15)
    n_permissive = cats.count("permissive")
    if n_permissive == 2:
        score += 15; reasons.append("both flanks permissive")
    elif n_permissive == 1:
        score += 8; reasons.append("one permissive flank")

    return round(score, 1), reasons


def analyze_phage(accession: str, name: str, cds_df: pd.DataFrame) -> list[dict]:
    """Score insertion sites for one phage using Pharokka annotations."""
    if cds_df.empty:
        return []

    cds_df = cds_df.sort_values("start").reset_index(drop=True)
    cds_df["category"] = cds_df["phrog_category"].apply(phrog_to_class)

    # Lysis gene midpoints (for distance calculation)
    lysis_positions = [
        int((r["start"] + r["end"]) / 2)
        for _, r in cds_df[cds_df["category"] == "lysis"].iterrows()
    ]

    sites = []
    for i in range(len(cds_df) - 1):
        left = cds_df.iloc[i]
        right = cds_df.iloc[i + 1]
        gap = int(right["start"]) - int(left["end"])
        if gap < MIN_GAP_BP:
            continue

        mid = int(left["end"]) + gap // 2
        dist_lysis = min((abs(mid - p) for p in lysis_positions), default=10_000_000)

        score, reasons = score_site(
            gap_bp=gap,
            flank_left_cat=left["category"],
            flank_right_cat=right["category"],
            distance_to_lysis_bp=dist_lysis,
        )

        sites.append(
            {
                "accession": accession,
                "name": name,
                "site_start": int(left["end"]),
                "site_end": int(right["start"]),
                "gap_bp": gap,
                "left_gene": str(left.get("phrog_annot", ""))[:80],
                "left_category": left["category"],
                "left_phrog_category": left["phrog_category"],
                "right_gene": str(right.get("phrog_annot", ""))[:80],
                "right_category": right["category"],
                "right_phrog_category": right["phrog_category"],
                "dist_to_lysis_bp": dist_lysis if dist_lysis < 1e6 else None,
                "insertion_score": score,
                "rationale": " · ".join(reasons),
            }
        )

    sites.sort(key=lambda s: -s["insertion_score"])
    return sites


def main() -> None:
    if not PHAROKKA_CSV.exists():
        raise SystemExit(f"❌ Pharokka annotations not found at {PHAROKKA_CSV}. "
                         "Run 02_pharokka_annotate.py first.")

    annotations = pd.read_csv(PHAROKKA_CSV)
    ranked = pd.read_csv(RANKED_CSV)
    old_sites = pd.read_csv(OLD_INSERTION_CSV) if OLD_INSERTION_CSV.exists() else pd.DataFrame()

    accessions = annotations["accession"].unique().tolist()
    print(f"\n[1/2] Re-scoring insertion sites for {len(accessions)} phages.\n")

    all_sites: list[dict] = []
    log_lines = ["Phase 2 Step 3: Pharokka-based insertion scoring",
                 "================================================"]
    log_lines.append("Before/after comparison (top site only):")
    log_lines.append("")

    for acc in accessions:
        row = ranked[ranked["accession"] == acc].iloc[0]
        cds_df = annotations[annotations["accession"] == acc]
        sites = analyze_phage(acc, row["name"], cds_df)
        all_sites.extend(sites)

        # Log: best site, with before/after comparison
        new_best = sites[0] if sites else None
        old_best = None
        if not old_sites.empty:
            old_match = old_sites[old_sites["accession"] == acc]
            if not old_match.empty:
                old_best = old_match.sort_values("insertion_score", ascending=False).iloc[0]

        log_lines.append(f"▸ {acc} — {row['name'][:50]}")
        if old_best is not None:
            log_lines.append(
                f"     BEFORE (Phase 2.1):  score={old_best['insertion_score']:.1f}   "
                f"flanks: {old_best['left_category']} | {old_best['right_category']}"
            )
        if new_best is not None:
            log_lines.append(
                f"     AFTER  (Pharokka):   score={new_best['insertion_score']:.1f}   "
                f"flanks: {new_best['left_category']} | {new_best['right_category']}"
            )
            log_lines.append(f"        left  ({new_best['left_phrog_category']}): {new_best['left_gene']}")
            log_lines.append(f"        right ({new_best['right_phrog_category']}): {new_best['right_gene']}")
            log_lines.append(f"        rationale: {new_best['rationale']}")
        log_lines.append("")

    out_df = pd.DataFrame(all_sites)
    out_df.to_csv(CSV_OUT, index=False)
    print(f"[2/2] Wrote {len(out_df)} re-scored sites to {CSV_OUT}\n")

    for line in log_lines:
        print("    " + line)
    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
