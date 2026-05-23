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
| **1. Atlas** | Which phage to engineer? | `src/atlas/` | 🚧 In progress |
| **2. Insertion Site Finder** | Where to insert luciferase? | `src/insertion/` | 📋 Planned |
| **3. Host Range Predictor** | Which strains can it detect? | `src/host_range/` | 📋 Planned |
| **4. Cocktail Designer** | What cocktail composition? | `src/cocktail/` | 📋 Planned |

---

## Phase 1: Listeria Phage Atlas

The first deliverable is a ranked, searchable atlas of every publicly sequenced *Listeria* phage, annotated with the features that matter for reporter-phage engineering:

- Genome size and GC content
- Lytic vs. temperate lifestyle (BACPHLIP classifier)
- Functional gene annotation (Pharokka)
- Reported host range (INPHARED metadata)
- Genome similarity clusters (mash distance)
- Engineering-readiness score (composite ranking)

**Output:** a CSV of ranked candidates and an interactive Streamlit dashboard.

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

Crystal (Yushan) Zhao — Founder & Team Lead, PHIND
Undergraduate researcher in phage biology, University of Illinois Urbana-Champaign
yushanz5@illinois.edu

## License

MIT
