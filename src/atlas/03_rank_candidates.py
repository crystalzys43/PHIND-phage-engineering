"""
Phase 1, Step 3: Rank lytic Listeria phages by engineering readiness.

WHY THIS STEP EXISTS
--------------------
After Step 2 we have 32 lytic Listeria phages — far too many to engineer
all at once. PHIND wet lab needs a prioritized list: "order phage #1 first,
then #2, etc." This step assigns each candidate a transparent score (0-100)
across five biological criteria, then sorts the candidates.

SCORING (composite 0-100)
-------------------------
1. annotation_quality (0-25): can we even SEE the genome's genes?
   - cds_count / expected_cds_for_genome_size, capped at 25
   - records with cds_count == 0 (e.g., UNVERIFIED) score 0

2. size_class_score (0-25): is the genome in the engineering sweet spot?
   - 60–200 kb → 25  (room for reporter insertion, established backbones)
   - 30–60 kb  → 12  (compact, fewer non-essential genes, riskier)
   - other     → 5

3. lytic_confidence (0-20): how cleanly lytic is this phage?
   - 20 if n_lysogeny_genes == 0
   - decreases by 5 per residual signal (defensive, since Step 2 already filtered)

4. taxonomy_score (0-20): is this in a reporter-phage-friendly clade?
   - Herelleviridae or Pecentumvirus       → 20  (A511, P100 family)
   - Other Caudoviricetes (tailed phages)  → 15
   - Anything else                          → 5

5. literature_bonus (0-10): published reporter-phage precedent?
   - A511 (Loessner 1996), P100 (LISTEX), P70 → 10
   - Otherwise 0

DEDUPLICATION
-------------
RefSeq (NC_xxxx) records are exact copies of GenBank originals (e.g.
NC_018831 == JX442241 for Listeria phage P70). We keep one per (length_bp,
name) cluster, preferring the original GenBank accession.

UNVERIFIED RECORDS
------------------
NCBI marks unannotated submissions with "UNVERIFIED:" in the title and
gives them 0 CDS features. These auto-score very low and sink to the
bottom, but we also flag them explicitly.

Usage
-----
    python src/atlas/03_rank_candidates.py

Inputs
------
    data/listeria_phages_classified.csv   (from Step 2)

Outputs
-------
    results/candidate_phages_ranked.csv   — ranked, deduplicated CSV
    results/03_ranking_log.txt            — run summary + top-10 preview
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CSV_IN = PROJECT_ROOT / "data" / "listeria_phages_classified.csv"
CSV_OUT = PROJECT_ROOT / "results" / "candidate_phages_ranked.csv"
LOG_OUT = PROJECT_ROOT / "results" / "03_ranking_log.txt"

# Known reporter-phage precedents in the published literature.
# Add to this list as Crystal does deeper literature review.
LITERATURE_PRECEDENTS = {
    "DQ003638": "A511 (Loessner et al. 1996, Mol Microbiol — original luxAB reporter)",
    "DQ004855": "P100 (FDA-approved LISTEX P100 biocontrol; Carlton et al. 2005)",
    "JX442241": "P70 (Klumpp et al. 2014 — group B Listeria phage characterization)",
}


# ---------------------------------------------------------------------------
# Scoring functions — each returns a float
# ---------------------------------------------------------------------------

def expected_cds(length_bp: int) -> float:
    """Phage genomes have roughly 1 CDS per ~1 kb. Use 1100 bp/CDS as a
    conservative expectation — better-annotated phages should still hit
    ratio close to 1 with this denominator."""
    return length_bp / 1100


def score_annotation(cds_count: int, length_bp: int) -> float:
    if cds_count == 0:
        return 0.0
    ratio = cds_count / expected_cds(length_bp)
    return min(1.0, ratio) * 25


def score_size(length_bp: int) -> float:
    if 60_000 <= length_bp <= 200_000:
        return 25.0
    if 30_000 <= length_bp < 60_000:
        return 12.0
    return 5.0


def score_lytic(n_lysogeny_genes: int) -> float:
    # Step 2 already filtered to lytic, so n_lysogeny_genes should be 0.
    # The decay protects against borderline cases if we ever loosen Step 2.
    return max(0.0, 20.0 - 5.0 * n_lysogeny_genes)


def score_taxonomy(taxonomy: str) -> float:
    tax_lower = taxonomy.lower()
    if "herelleviridae" in tax_lower or "pecentumvirus" in tax_lower:
        return 20.0
    if "caudoviricetes" in tax_lower:
        return 15.0
    return 5.0


def score_literature(accession: str) -> tuple[float, str]:
    """Return (bonus_points, citation_string)."""
    base = accession.split(".")[0]
    if base in LITERATURE_PRECEDENTS:
        return 10.0, LITERATURE_PRECEDENTS[base]
    return 0.0, ""


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop RefSeq/GenBank duplicates. We cluster by (length_bp, phage_name_root)
    and keep the row with the shortest accession prefix that is NOT a
    RefSeq (NC_*) accession, when available.
    """
    # Phage "name root" = the description up to the first comma
    df = df.copy()
    df["name_root"] = df["name"].str.split(",").str[0].str.strip()
    df["is_refseq"] = df["accession"].str.startswith("NC_")

    # Sort so non-RefSeq comes first (False sorts before True)
    df = df.sort_values(["length_bp", "name_root", "is_refseq"])
    deduped = df.drop_duplicates(subset=["length_bp", "name_root"], keep="first")
    return deduped.drop(columns=["name_root", "is_refseq"])


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    df = pd.read_csv(CSV_IN)
    print(f"\n[1/4] Loaded {len(df)} classified Listeria phages.")

    # Keep only lytic candidates
    lytic = df[df["lifestyle"] == "lytic"].copy()
    print(f"      → {len(lytic)} lytic candidates pre-dedup.")

    # De-duplicate (RefSeq vs GenBank)
    deduped = deduplicate(lytic)
    print(f"[2/4] After deduplication: {len(deduped)} unique candidates.")

    # Compute each scoring component
    print(f"[3/4] Scoring on 5 dimensions...")
    deduped["annotation_quality"] = deduped.apply(
        lambda r: score_annotation(r["cds_count"], r["length_bp"]), axis=1
    )
    deduped["size_class_score"] = deduped["length_bp"].apply(score_size)
    deduped["lytic_confidence"] = deduped["n_lysogeny_genes"].apply(score_lytic)
    deduped["taxonomy_score"] = deduped["taxonomy"].apply(score_taxonomy)
    deduped[["literature_bonus", "literature_citation"]] = deduped["accession"].apply(
        lambda a: pd.Series(score_literature(a))
    )

    # Total
    deduped["engineering_readiness_score"] = (
        deduped["annotation_quality"]
        + deduped["size_class_score"]
        + deduped["lytic_confidence"]
        + deduped["taxonomy_score"]
        + deduped["literature_bonus"]
    ).round(1)

    # Sort descending by score
    ranked = deduped.sort_values(
        "engineering_readiness_score", ascending=False
    ).reset_index(drop=True)
    ranked["rank"] = ranked.index + 1

    # Reorder columns: rank, score, identity, score components, metadata
    front = [
        "rank",
        "engineering_readiness_score",
        "accession",
        "name",
        "length_bp",
        "gc_percent",
        "cds_count",
        "annotation_quality",
        "size_class_score",
        "lytic_confidence",
        "taxonomy_score",
        "literature_bonus",
        "literature_citation",
    ]
    other = [c for c in ranked.columns if c not in front]
    ranked = ranked[front + other]

    ranked.to_csv(CSV_OUT, index=False)
    print(f"[4/4] Wrote {CSV_OUT}")

    # ---- Log + top-10 preview ----
    top10 = ranked.head(10)[
        ["rank", "engineering_readiness_score", "accession", "length_bp", "name", "literature_citation"]
    ]

    log = ["Step 3 ranking summary", "======================"]
    log.append(f"Input lytic candidates:       {len(lytic)}")
    log.append(f"After deduplication:          {len(deduped)}")
    log.append(f"Output: {CSV_OUT}")
    log.append("")
    log.append("Top 10 engineering-ready candidates:")
    log.append("-" * 80)
    for _, r in top10.iterrows():
        lit = f"   ✦ {r['literature_citation']}" if r["literature_citation"] else ""
        log.append(
            f"  #{r['rank']:2d}  score={r['engineering_readiness_score']:5.1f}   "
            f"{r['accession']:14s} {r['length_bp']:>7,} bp   {r['name'][:55]}"
        )
        if lit:
            log.append("       " + lit)

    print("\n" + "\n".join(log))
    with LOG_OUT.open("w") as fh:
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
