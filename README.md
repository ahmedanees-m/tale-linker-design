# tale_linker_design

[![CI](https://github.com/ahmedanees-m/tale-linker-design/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ahmedanees-m/tale-linker-design/actions)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Geometric constraints on catalytic domain fusion to TALE arrays: a design framework for programmable genome engineering.**

This package provides the first systematic computational framework for designing peptide linkers that connect Transcription Activator-Like Effector (TALE) DNA-binding arrays to catalytic effector domains. It maps the three-dimensional "reachability envelope" of TALE C-terminal fusions as a function of linker length, composition, and flexibility—enabling rational design of TALE-nuclease, TALE-recombinase, and TALE-editor architectures.

---

## Why This Matters

Programmable genome engineering relies on fusing DNA-binding modules to catalytic effectors. While TALE arrays provide programmable recognition of any 20-bp sequence, the geometric constraints of these fusions have remained empirically determined—linker lengths are chosen by trial and error rather than physical principles.

This framework replaces guesswork with quantitative prediction. By modeling linkers as worm-like chains (flexible), helical rods (rigid), or composite sequences, we compute the probability that a catalytic domain placed at the TALE C-terminus can reach a target DNA position. The method enables:

- **Rational linker design** for new TALE-fusion architectures (nucleases, integrases, base editors)
- **Conditional assembly engineering** (ensuring catalytic domains are positioned only at target-bound configurations)
- **Avoidance of steric clashes** between TALE superhelices and effector domains

The framework was developed to support [GENESIS](link-to-paper-when-live), a program to design de novo DNA-modifying enzymes, but applies broadly to any TALE-fusion system.

---

## Installation

```bash
# Clone repository
git clone https://github.com/ahmedanees-m/tale-linker-design.git
cd tale-linker-design

# Create environment
conda env create -f environment.yml
conda activate tale-linker-design

# Install package
pip install -e .
```

---

## Quick Start

```python
import tale_linker_design as tld

# Load reference TALE-DNA structure
tale = tld.load_reference("3V6T")

# Define target: bp +4 on top strand (16.2 Å from TALE C-terminus)
target = tld.Target(strand="top", bp_offset=4)

# Get linker recommendations
candidates = tld.recommend_linkers(tale, target, top_k=3)

# Results: ranked by probability of reaching target within 5 Å
for link in candidates:
    print(f"{link['class']}-{link['length']}: P(reach) = {link['probability']:.1%}")
```

**Output:**
```
H-15: P(reach) = 9.9%   # Helical (EAAAK)×3
H-18: P(reach) = 5.9%   # Helical (EAAAK)×3.6
P-8:  P(reach) = 4.0%   # Proline-rich
```

---

## Key Features

- **Five linker classes**: Flexible (GGS), Helical (EAAAK), Mixed, Natural (TALEN variants), and Proline-rich
- **45 pre-computed ensembles**: 5–30 residue lengths, 50,000 conformations each
- **Canonical coordinate frame**: DNA-centric reference system for comparing across TALE structures
- **Inverse design**: Specify target position → get optimal linker sequence
- **Validation**: Benchmarked against published TALEN-FokI, TALE-PvuII, and cTALEN geometries

---

## Reproducing Results

To regenerate all reachability maps and linker statistics (requires ~250 CPU-hours):

```bash
python scripts/run_full_analysis.py --no-download
```

Large data files (PDB structures, linker ensembles, reachability maps) are deposited on Zenodo (DOI: [pending]). Small reference tables are included in `data/`.

---

## License

MIT License — see [LICENSE](LICENSE) file.

Copyright (c) 2026 Anees Ahmed Mahaboob Ali
