#!/usr/bin/env python3
"""
Master analysis script for C5: runs Steps 2–7 end-to-end.

Usage:
    python scripts/run_full_analysis.py [--pdb-dir data/pdb_cleaned] [--no-download]

Runtime on Intel i3 / 8 GB RAM: ~5–15 minutes (analytical models, no PyRosetta).
"""

import argparse
import json
import time
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# Add src to path if running directly
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import (
    load_reference, PRIORITY_STRUCTURES, download_priority_structures, clean_structure
)
from tale_linker_design.frames import (
    ReferenceFrame, build_scissile_phosphate_table, genesis_target_positions
)
from tale_linker_design.linkers import build_linker_library, linker_summary_dataframe
from tale_linker_design.reachability import (
    compute_all_reachability_maps, genesis_target_reachability
)
from tale_linker_design.design import (
    recommend_linkers, genesis_linker_specification, PUBLISHED_FUSIONS
)
from tale_linker_design.visualize import (
    figure1_reference_frame, figure2_linker_library, figure3_reachability_maps,
    figure4_published_benchmark, figure5_genesis_recommendations,
    figure6_decision_flowchart,
)

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb-dir", default=str(DATA_DIR / "pdb_cleaned"))
    p.add_argument("--no-download", action="store_true",
                   help="Skip PDB download and use literature parameters instead")
    p.add_argument("--n-samples", type=int, default=50_000,
                   help="MC samples per reachability map (default 50 000)")
    p.add_argument("--primary-structure", default="3V6T",
                   help="PDB ID to use as primary reference (default 3V6T)")
    p.add_argument("--save-ensembles", action="store_true",
                   help="Save the 3D numpy arrays to data/linker_library/ for Zenodo deposition")
    return p.parse_args()


def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)
    t0 = time.time()

    print("=" * 65)
    print("C5 ANALYSIS: Geometric Constraints on TALE-Fusion Linker Design")
    print("=" * 65)

    # ── Step 2: Load / prepare structures ──────────────────────────────
    print("\n[Step 2] Loading TALE-DNA reference structures...")
    if not args.no_download:
        try:
            raw_dir = DATA_DIR / "pdb_raw"
            downloaded = download_priority_structures(raw_dir)
            print(f"  Downloaded {len(downloaded)} PDB files")
            for pdb_id in downloaded:
                try:
                    clean_structure(pdb_id, raw_dir, pdb_dir)
                    print(f"  Cleaned {pdb_id}")
                except Exception as exc:
                    print(f"  Warning: cleaning {pdb_id} failed: {exc}")
        except Exception as exc:
            print(f"  PDB download failed ({exc}); using literature parameters.")

    primary_id = args.primary_structure
    tale = load_reference(primary_id, cleaned_dir=pdb_dir)
    print(f"  Loaded primary structure: {tale.pdb_id} "
          f"({tale.n_repeats:.1f} repeats, {tale.resolution:.2f} Å)")

    # Save annotation JSON
    ann_dir = DATA_DIR / "pdb_annotations"
    ann_dir.mkdir(parents=True, exist_ok=True)
    (ann_dir / f"{primary_id}.json").write_text(
        json.dumps(tale.to_annotation_json(), indent=2)
    )

    # ── Step 3: Reference frame and scissile phosphates ────────────────
    print("\n[Step 3] Computing reference frame and scissile phosphate table...")
    frame = ReferenceFrame.from_tale_structure(tale)
    scissile_table = build_scissile_phosphate_table(tale, frame, bp_range=11)
    print(f"  Frame origin: {frame.origin.round(2)}")
    print(f"  Scissile phosphates defined: {len(scissile_table)} positions")

    # Save CSV
    rows = []
    for (strand, bp), entry in scissile_table.items():
        c = entry["coords"]
        rows.append({
            "strand": strand, "bp_offset": bp,
            "x_A": round(c[0], 3), "y_A": round(c[1], 3), "z_A": round(c[2], 3),
            "distance_A": round(entry["distance_from_origin"], 3),
        })
    sc_df = pd.DataFrame(rows)
    sc_df.to_csv(DATA_DIR / "scissile_phosphate_coordinates.csv", index=False)
    print(f"  Saved: data/scissile_phosphate_coordinates.csv")

    # ── Step 4: Linker library ──────────────────────────────────────────
    print("\n[Step 4] Building analytical linker library (WLC + helical models)...")
    library = build_linker_library()
    summary_df = linker_summary_dataframe(library)
    summary_df.to_csv(DATA_DIR / "linker_summary_statistics.csv", index=False)
    print(f"  Built {len(library)} linker ensembles (5 classes × 9 lengths)")
    print(f"  Saved: data/linker_summary_statistics.csv")

    # ── Step 5: Reachability maps ──────────────────────────────────────
    print(f"\n[Step 5] Computing reachability maps ({args.n_samples:,} samples each)...")
    rng = np.random.default_rng(42)
    reach_maps = compute_all_reachability_maps(tale, library, n_samples=args.n_samples, rng=rng)
    print(f"  Computed {len(reach_maps)} reachability maps")

    if args.save_ensembles:
        library_dir = DATA_DIR / "linker_library"
        library_dir.mkdir(parents=True, exist_ok=True)
        for (cls, n), rm in reach_maps.items():
            if hasattr(rm, "samples"):
                np.save(library_dir / f"{cls.value}_{n}.npy", rm.samples)
        print(f"  Saved 3D coordinate ensembles to: data/linker_library/*.npy")

    # Report summary statistics
    reach_df = genesis_target_reachability(reach_maps, scissile_table, tolerance_A=5.0)
    reach_df.to_csv(DATA_DIR / "reachability_maps" / "genesis_target_reachability.csv", index=False)
    print(f"  Saved: data/reachability_maps/genesis_target_reachability.csv")

    top_results = reach_df[reach_df["bp_offset"] == 4].sort_values("p_reach_pct", ascending=False)
    print("\n  Top linker designs for primary GENESIS target (bp +4, top strand):")
    print(top_results[top_results["strand"] == "top"].head(5).to_string(index=False))

    # ── Step 7: GENESIS-specific recommendations ───────────────────────
    print("\n[Step 7] Generating GENESIS linker recommendations...")
    from tale_linker_design.frames import Target
    primary_target = Target(strand="top", bp_offset=4)
    designs = recommend_linkers(
        tale, primary_target, reach_maps, library, scissile_table, top_k=3
    )

    for i, d in enumerate(designs):
        print(f"\n  Design {i+1}: {d.name}")
        print(f"    Sequence motif: {d.sequence_motif}")
        print(f"    Sequence:       {d.recommended_sequence[:30]}...")
        print(f"    P(reach):       {d.p_reach_primary_pct:.1f}%")
        print(f"    Nearest dist:   {d.nearest_A:.1f} Å")
        print(f"    Entropic cost:  {d.entropic_cost_kcal_mol:.1f} kcal/mol")

    spec = genesis_linker_specification(designs, primary_target)
    spec_path = DATA_DIR / "genesis_linker_recommendations.json"
    spec_path.write_text(json.dumps(spec, indent=2))
    print(f"\n  Saved: data/genesis_linker_recommendations.json")

    # ── Step 9: Generate all figures ──────────────────────────────────
    print("\n[Step 9] Generating publication figures...")

    figure1_reference_frame(scissile_table, FIG_DIR / "fig1_reference_frame.png")
    print("  Fig 1: reference frame — done")

    figure2_linker_library(library, FIG_DIR / "fig2_linker_library.png")
    print("  Fig 2: linker library — done")

    figure3_reachability_maps(reach_maps, scissile_table, FIG_DIR / "fig3_reachability.png")
    print("  Fig 3: reachability maps — done")

    published_list = [
        {**v, "name": k, "cut_site_bp_offset": min(v["cut_site_bp_offset"], 10)}
        for k, v in PUBLISHED_FUSIONS.items()
    ]
    figure4_published_benchmark(published_list, reach_maps, scissile_table,
                                FIG_DIR / "fig4_benchmark.png")
    print("  Fig 4: published benchmark — done")

    figure5_genesis_recommendations(designs, reach_maps, scissile_table,
                                    FIG_DIR / "fig5_genesis_designs.png")
    print("  Fig 5: GENESIS designs — done")

    figure6_decision_flowchart(FIG_DIR / "fig6_flowchart.png")
    print("  Fig 6: decision flowchart — done")

    # ── Task 7.1: Provenance manifest ─────────────────────────────────
    print("\n[Step 10] Generating provenance manifest...")
    manifest = {
        "seed": 42,
        "n_samples": args.n_samples,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "primary_structure": args.primary_structure,
        "md5_hashes": {}
    }
    
    for f in list(DATA_DIR.rglob("*")) + list(FIG_DIR.rglob("*")):
        if f.is_file() and f.name != "provenance_manifest.json":
            hasher = hashlib.md5()
            with open(f, "rb") as fh:
                for chunk in iter(lambda: fh.read(4096), b""):
                    hasher.update(chunk)
            key = f.relative_to(BASE_DIR).as_posix()
            manifest["md5_hashes"][key] = hasher.hexdigest()

    manifest_path = DATA_DIR / "provenance_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print("  Saved: data/provenance_manifest.json")

    elapsed = time.time() - t0
    print(f"\n{'='*65}")
    print(f"Analysis complete in {elapsed:.1f} s")
    print(f"All outputs in: {DATA_DIR.resolve()}")
    print(f"All figures in: {FIG_DIR.resolve()}")
    print("=" * 65)


if __name__ == "__main__":
    main()
