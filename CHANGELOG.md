# Changelog

All notable changes to `tale_linker_design` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- PyPI publication after bioRxiv goes live
- Extended Monte Carlo ensemble for supplementary data (HPC phase)
- AlphaFold3-based benchmark of published fusion systems

---

## [0.1.0] — 2026-04-19 — Initial preprint release

### Added
- `tale_linker_design` Python package (6 modules):
  - `structures.py` — PDB download, cleaning, chain extraction, annotation
  - `frames.py` — canonical reference frame, scissile-phosphate coordinate table
  - `linkers.py` — analytical WLC + helical-rod linker library (45 ensembles)
  - `reachability.py` — Monte Carlo reachability maps, P(reach) computation
  - `design.py` — inverse design (target position → ranked linker recommendations)
  - `visualize.py` — 6-panel publication figure generator
- Full analysis pipeline (`scripts/run_full_analysis.py`) with `--no-download` flag
- PDB download and cleaning scripts (`scripts/download_pdbs.py`, `scripts/clean_pdbs.py`)
- 24 unit tests across 3 test files; all pass on Python 3.12
- 6 Jupyter notebooks (Steps 01–07)
- Complete LaTeX manuscript (`manuscript/main.tex`, `supplementary.tex`, `references.bib`)
- Literature survey CSV with 50 TALE-fusion entries (`data/literature_survey.csv`)
- GENESIS linker specification (`data/genesis_linker_recommendations.json`)
- GitHub Actions CI (`python 3.11`, `pytest`)

### Scientific results
- Primary GENESIS target: bp +4 top strand, 16.2 Å from TALE C-terminus
- Design 1: (EAAAK)₂ 10-res helical, P(reach, 12 Å) = 19.1%
- Design 2: (GGS)₃GG 10-res flexible, P(reach, 12 Å) = 11.5%
- Design 3: GGS-EAAAK-GGS 10-res mixed, P(reach, 12 Å) = 9.7%
- New finding: length-matching principle (optimal n ≈ d_target / rise_per_residue)

### Known limitations
- WLC model uses analytical approximation; explicit-solvent MD recommended pre-synthesis
- Steric exclusion uses cylindrical analytical model, not all-atom clash detection
- AlphaFold3 benchmark deferred to HPC phase

---

[Unreleased]: https://github.com/anees-ahmed/tale-linker-design/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/anees-ahmed/tale-linker-design/releases/tag/v0.1.0
