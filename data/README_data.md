# Data Directory — C5 Analysis Outputs

This directory contains all data files produced by the C5 analysis pipeline.

## Directory Structure

```
data/
├── pdb_raw/                        # Raw downloaded PDB/mmCIF files (8 structures)
├── pdb_cleaned/                    # Biopython-cleaned PDB files (solvent/altconf removed)
├── pdb_annotations/                # JSON structural annotations per structure
├── reachability_maps/              # Per-linker reachability statistics CSV files
├── linker_summary_statistics.csv   # Summary of all 45 linker ensembles
├── literature_survey.csv           # Literature survey (16 TALE-fusion papers)
├── scissile_phosphate_coordinates.csv  # Canonical-frame phosphate positions
└── genesis_linker_recommendations.json # Frozen GENESIS Design 1/2/3 specification
```

## Key Files

### `genesis_linker_recommendations.json`
The **frozen GENESIS linker specification** — the primary deliverable of this paper.
Contains Design 1 (H-class helical EAAAK), Design 2 (GGS flexible), and Design 3 (Mixed).
This file is the authoritative recommendation for the GENESIS research programme.

**Design 1 (Primary Recommendation):** 15-residue `(EAAAK)₃` helical linker
- P(reach) = 9.9% at the GENESIS primary target (top strand, bp +4)
- Nearest approach: 0.4 Å  
- Entropic cost: 8.7 kcal/mol (helical ordering penalty compared to GGS)

### `scissile_phosphate_coordinates.csv`
Canonical-frame (Ångström) coordinates of scissile phosphate positions derived from
the 3V6T crystal structure (1.85 Å resolution). 22 positions: 11 bp on each strand.

### `linker_summary_statistics.csv`
WLC and helical-rod model statistics for all 45 linker ensembles:
5 classes (F/N/P/H/M) × 9 lengths (4–20 residues)

### `reachability_maps/`
Per-linker reachability statistics (50,000 WLC samples each):
- `genesis_target_reachability.csv` — P(reach) matrix for all GENESIS target positions

### `literature_survey.csv`
Survey of 16 published TALE-fusion and CRISPR-fusion papers with experimental data
on linker sequences, lengths, efficiencies, and off-target rates.

## Reproducing the Data

```bash
# From the repository root:
python scripts/download_pdbs.py          # Download PDB files (~3 min, network)
python scripts/clean_pdbs.py             # Clean structures (~10 sec)
python scripts/run_full_analysis.py --no-download  # Full analysis (~4 min on i3)
```

All scripts are deterministic with rng seed=42. Results should be numerically identical
across runs on any platform.

## Frozen Outputs

The following files are **version-controlled** and must not be overwritten without
updating the manuscript accordingly:
- `genesis_linker_recommendations.json`
- `scissile_phosphate_coordinates.csv`
- `linker_summary_statistics.csv`
