#!/usr/bin/env python3
"""
Download priority TALE-DNA crystal structures from the RCSB PDB.

Usage:
    python scripts/download_pdbs.py [--output-dir data/pdb_raw]

Downloads 7 priority structures listed in C5 Execution Plan Appendix B.
Files are saved as <PDB_ID>.pdb in the output directory.
"""

import argparse
import sys
from pathlib import Path

# Add src to path if running directly
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import PRIORITY_STRUCTURES, download_priority_structures


def parse_args():
    p = argparse.ArgumentParser(description="Download priority TALE-DNA PDB structures.")
    p.add_argument(
        "--output-dir",
        default="data/pdb_raw",
        help="Directory to save PDB files (default: data/pdb_raw)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    base_dir = Path(__file__).parent.parent
    out_dir = base_dir / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("C5 PDB Download: Priority TALE-DNA Crystal Structures")
    print("=" * 60)
    print(f"\nTarget directory: {out_dir.resolve()}")
    print(f"Structures to download: {', '.join(PRIORITY_STRUCTURES)}\n")

    downloaded = download_priority_structures(out_dir)

    print(f"\nDone. Downloaded {len(downloaded)}/{len(PRIORITY_STRUCTURES)} structures.")
    for pdb_id in PRIORITY_STRUCTURES:
        pdb_file = out_dir / f"{pdb_id.lower()}.pdb"
        if pdb_file.exists():
            size_kb = pdb_file.stat().st_size // 1024
            print(f"  [OK] {pdb_id:6s}  {size_kb} KB  {pdb_file}")
        else:
            print(f"  ✗ {pdb_id:6s} — NOT FOUND (download may have failed)")

    print("\nNext step: python scripts/clean_pdbs.py")


if __name__ == "__main__":
    main()
