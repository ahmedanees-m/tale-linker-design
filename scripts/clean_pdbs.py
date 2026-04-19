#!/usr/bin/env python3
"""
Clean and standardize downloaded TALE-DNA PDB structures.

Usage:
    python scripts/clean_pdbs.py [--raw-dir data/pdb_raw] [--clean-dir data/pdb_cleaned]

For each PDB file in raw-dir:
  - Removes waters, buffer ligands, alternate conformations
  - Separates TALE protein chain(s) and DNA strand(s)
  - Applies consistent residue numbering
  - Saves cleaned PDB to clean-dir
  - Saves annotation JSON to data/pdb_annotations/

Run after: python scripts/download_pdbs.py
"""

import argparse
import json
import sys
from pathlib import Path

# Add src to path if running directly
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import (
    PRIORITY_STRUCTURES,
    clean_structure,
    load_reference,
)


def parse_args():
    p = argparse.ArgumentParser(description="Clean priority TALE-DNA PDB structures.")
    p.add_argument("--raw-dir", default="data/pdb_raw")
    p.add_argument("--clean-dir", default="data/pdb_cleaned")
    p.add_argument("--annotations-dir", default="data/pdb_annotations")
    return p.parse_args()


def main():
    args = parse_args()
    base_dir = Path(__file__).parent.parent
    raw_dir = base_dir / args.raw_dir
    clean_dir = base_dir / args.clean_dir
    ann_dir = base_dir / args.annotations_dir

    clean_dir.mkdir(parents=True, exist_ok=True)
    ann_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("C5 PDB Cleaning: Structural Preparation Pipeline")
    print("=" * 60)
    print(f"\nRaw directory:     {raw_dir.resolve()}")
    print(f"Cleaned directory: {clean_dir.resolve()}")
    print(f"Annotations dir:   {ann_dir.resolve()}\n")

    success, failure = [], []

    for pdb_id in PRIORITY_STRUCTURES:
        # Check for raw file (Biopython saves as lowercase)
        raw_candidates = [
            raw_dir / f"{pdb_id.lower()}.pdb",
            raw_dir / f"{pdb_id.upper()}.pdb",
            raw_dir / f"pdb{pdb_id.lower()}.ent",
        ]
        raw_file = next((f for f in raw_candidates if f.exists()), None)

        if raw_file is None:
            print(f"  [SKIP] {pdb_id}: raw file not found in {raw_dir} -- skipping")
            failure.append(pdb_id)
            continue

        try:
            clean_pdb = clean_structure(pdb_id, raw_dir, clean_dir)
            # Load cleaned structure to generate annotation JSON
            tale = load_reference(pdb_id, cleaned_dir=clean_dir)
            ann_path = ann_dir / f"{pdb_id}.json"
            ann_path.write_text(json.dumps(tale.to_annotation_json(), indent=2))
            print(f"  [OK] {pdb_id}: cleaned -> {clean_pdb.name} | annotation -> {ann_path.name}")
            success.append(pdb_id)
        except Exception as exc:
            print(f"  [FAIL] {pdb_id}: cleaning failed -- {exc}")
            failure.append(pdb_id)

    print(f"\nDone. {len(success)}/{len(PRIORITY_STRUCTURES)} structures cleaned successfully.")
    if failure:
        print(f"Failed: {', '.join(failure)}")
    print("\nNext step: python scripts/run_full_analysis.py --no-download")


if __name__ == "__main__":
    main()
