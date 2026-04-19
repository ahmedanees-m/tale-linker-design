import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
from tale_linker_design.structures import load_reference
from tale_linker_design.frames import Target, ReferenceFrame, build_scissile_phosphate_table
from tale_linker_design.linkers import build_linker_library, LinkerClass
from tale_linker_design.reachability import compute_reachability

def main():
    print("Running baseline comparisons...")
    tale = load_reference("3V6T", cleaned_dir="data/pdb_cleaned")
    frame = ReferenceFrame.from_tale_structure(tale)
    scissile_table = build_scissile_phosphate_table(tale, frame, bp_range=11)
    target_coords = scissile_table[("top", 4)]["coords"]

    rng = np.random.default_rng(42)

    # 1. Baseline: TALEN standard linker (n=63, flexible)
    lib_talen = build_linker_library(classes=[LinkerClass.F], lengths=[63])
    ens_talen = lib_talen[(LinkerClass.F, 63)]
    rmap_talen = compute_reachability(tale, ens_talen, n_samples=50_000, rng=rng)
    p_talen = rmap_talen.bootstrap_p_reach(target_coords, tolerance_A=12.0)["estimate"]
    
    # 2. Baseline: GGS n=10 (Flexible matching the optimal length)
    lib_ggs = build_linker_library(classes=[LinkerClass.F], lengths=[10])
    ens_ggs = lib_ggs[(LinkerClass.F, 10)]
    rmap_ggs = compute_reachability(tale, ens_ggs, n_samples=50_000, rng=rng)
    p_ggs = rmap_ggs.bootstrap_p_reach(target_coords, tolerance_A=12.0)["estimate"]

    # 3. GENESIS Design 1: (EAAAK)2 n=10
    lib_h = build_linker_library(classes=[LinkerClass.H], lengths=[10])
    ens_h = lib_h[(LinkerClass.H, 10)]
    rmap_h = compute_reachability(tale, ens_h, n_samples=50_000, rng=rng)
    p_h = rmap_h.bootstrap_p_reach(target_coords, tolerance_A=12.0)["estimate"]
    
    print("\nTarget: top strand, bp +4, 16.2 \u00C5 separation")
    print(f"TALEN standard (GGS, n=63): P(reach, 12 \u00C5) = {p_talen:5.1%}")
    print(f"Flexible matched (GGS, n=10): P(reach, 12 \u00C5) = {p_ggs:5.1%}")
    print(f"GENESIS Design 1 (H, n=10): P(reach, 12 \u00C5)   = {p_h:5.1%}")

if __name__ == "__main__":
    main()
