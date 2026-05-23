# PHIND Phage Engineering Pipeline

> An open-source pipeline that helps phage engineers decide **which phage to engineer**, **where to insert reporter genes**, **which bacterial strains it will detect**, and **how to design optimal phage cocktails** — built to support [PHIND](https://example.com)'s reporter phage development.

PHIND is a phage-based bacterial contamination screening platform for food manufacturing QA. Building it requires engineering **reporter phages** — phages modified to produce a measurable signal (luciferase) when they infect target food pathogens. Before any wet-lab work begins, four decisions have to be made:

1. **Which phage** should we engineer as the reporter backbone?
2. **Where** in its genome should we insert the luciferase reporter?
3. **Which bacterial strains** will the engineered phage actually detect?
4. **What cocktail** of phages should the PHIND cartridge contain?

This pipeline answers each of those questions using publicly available phage genomic data.

---

## Project Roadmap

| Phase | PHIND Question | Module | Status |
|---|---|---|---|
| **1. Atlas** | Which phage to engineer? | `src/atlas/` | ✅ Complete |
| **2. Insertion Site Finder** | Where to insert luciferase? | `src/insertion/` | 🚧 v1 complete; Pharokka re-annotation planned |
| **3. Host Range Predictor** | Which strains can it detect? | `src/host_range/` | 📋 Planned |
| **4. Cocktail Designer** | What cocktail composition? | `src/cocktail/` | 📋 Planned |

### Interactive dashboard

```bash
streamlit run dashboard/app.py
```

Opens a local browser dashboard with five tabs: ranked candidates, score breakdown, genome landscape, phage detail card, and Phase 2 insertion sites with linear genome maps.

---

## Phase 1: Listeria Phage Atlas

The first deliverable is a ranked, searchable atlas of every publicly sequenced *Listeria* phage, annotated with the features that matter for reporter-phage engineering:

- Genome size and GC content
- Lytic vs. temperate lifestyle (transparent gene-content heuristic)
- Functional gene annotation (planned: Pharokka)
- Reported host range (planned: INPHARED metadata)
- Genome similarity clusters (planned: mash distance)
- Engineering-readiness score (planned: composite ranking)

**Output:** a CSV of ranked candidates and an interactive Streamlit dashboard.

### Progress so far

| Step | Script | Records | Outcome |
|---|---|---|---|
| 1. Collect genomes | `src/atlas/01_collect_genomes.py` | 230 candidate phage genomes from NCBI | `data/listeria_phages.csv` |
| 2. Filter + classify lifestyle | `src/atlas/02_classify_lifestyle.py` | 230 → 128 Listeria hosts → **32 lytic reporter candidates** | `data/listeria_phages_classified.csv` |
| 3. Rank by engineering readiness | `src/atlas/03_rank_candidates.py` | 32 lytic → 25 unique → **ranked 1–25** | `results/candidate_phages_ranked.csv` |

### Phase 1 final output: top 5 PHIND reporter-phage candidates

| Rank | Score | Accession | Phage | Why this matters |
|---:|---:|---|---|---|
| 1 | 100/100 | DQ003638.2 | **A511** | Loessner 1996 — the original luxAB reporter phage backbone |
| 2 | 100/100 | DQ004855.1 | **P100** | FDA-approved LISTEX biocontrol; large existing safety dossier |
| 3 | 95/100 | JX442241.1 | P70 | Group B Listeria phage; Klumpp 2014 characterization |
| 4 | 90/100 | MN939539.1 | vB_Lino_VEfB7 | Recently characterized giant; similar architecture to A511 |
| 5 | 85/100 | JX126919.1 | LP-110 | Mid-size, broadly representative of LP-series |

PHIND wet-lab implication: start with **A511** as primary backbone (deepest reporter-engineering literature), **P100** as secondary (regulatory advantages), and the **LP-series** as cocktail-diversity options.

---

## Phase 2: Insertion Site Finder (v1)

For the top-10 ranked phages, Phase 2 identifies intergenic regions where a luciferase cassette could be inserted without disrupting the lytic cycle.

### Scoring biology

Each candidate site scores 0–100 across four criteria:

- **Gap size** (≥1 kb comfortable; ≥500 bp tight; ≥50 bp minimal)
- **Flank essentiality** (neither essential = best; both essential = excluded)
- **Lysis-cluster proximity** (adjacent to endolysin/holin = preferred; this is the published Loessner 1996 A511-*luxAB* insertion locale)
- **Permissive-flank bonus** (hypothetical-protein flanks are safest)

Gene categories (essential / lysis / permissive / other) are assigned by keyword matching on NCBI CDS product names.

### Findings & known limitation

- For phages with rich functional annotation (P70, LP-series, vB_Lino_VEfB7), Phase 2 successfully identifies sites adjacent to endolysin / lysis cluster — the textbook Loessner insertion locale.
- For **A511 and P100**, the NCBI records use minimal gene-product names (`gp1`, `gp2`, …) without functional descriptions, so our keyword classifier cannot categorize flanks. Their top sites cap at 55/100 — a *data-annotation* artifact, not a biological deficiency.
- **Phase 2 next step:** re-annotate the top-5 candidates with [Pharokka](https://github.com/gbouras13/pharokka) to recover functional categories and produce calibrated insertion scores across all candidates.

### Step 2 sanity check

The Step 2 classifier rediscovered **A511** and **P100** as top candidates from the raw public data — these are the same two phages the Listeria reporter-phage literature has converged on since 1996. This validates that the lifestyle-classification heuristic captures the right biology.

**Top 5 lytic Listeria phages by genome size:**

| Accession | Size (bp) | Phage | Note |
|---|---:|---|---|
| OQ999172.1 | 181,606 | LIS04 | Giant Herelleviridae |
| DQ003638.2 | 137,619 | **A511** | Loessner 1996 reporter backbone |
| MN939539.1 | 135,461 | vB_Lino_VEfB7 | Recent giant |
| DQ004855.1 | 131,384 | **P100** | FDA-approved LISTEX biocontrol |
| OK283618.1 | 87,038 | LPML1 | Mid-size |

---

## Why this exists

PHIND won the Reimagine Our Future Grand Prize (Dec 2025) and the Cozad New Venture Challenge Agriculture Innovation Prize (Spring 2026). The wet-lab phase begins Summer 2026, and we want every reagent we order to be backed by evidence — not chosen by convenience. This pipeline is that evidence layer.

---

## Quick start

```bash
# clone
git clone https://github.com/crystalzys43/PHIND-phage-engineering.git
cd PHIND-phage-engineering

# install
pip install -r requirements.txt

# run Phase 1 step 1: collect Listeria phage genomes
python src/atlas/01_collect_genomes.py
```

---

## Author

**Crystal (Yushan) Zhao** — Founder & Team Lead, PHIND
Undergraduate researcher in phage biology, University of Illinois Urbana-Champaign
yushanz5@illinois.edu

I designed this pipeline, made all scientific decisions (target pathogen, ranking criteria, classification heuristic, deduplication strategy), and own the interpretation of every result. I am happy to walk through any step in detail.

## How this was built

This project was developed using **Claude (Anthropic) as a pair-programming partner**. I drove the scientific design and the engineering decisions; Claude helped translate those decisions into Python and explained library choices as we went. Individual commits from Step 2 onward use the `Co-authored-by` trailer; the Step 1 setup commit predates that practice but the same pairing applied.

This is, deliberately, an honest portfolio of how a modern undergraduate scientist works in 2026: bringing biological judgment to the front, using AI tooling to amplify implementation throughput, and being transparent about both.

## License

MIT
