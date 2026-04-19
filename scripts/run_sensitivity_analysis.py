import sys
import pandas as pd
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import load_reference
from tale_linker_design.frames import ReferenceFrame, build_scissile_phosphate_table
from tale_linker_design.linkers import build_linker_library, LinkerClass
from tale_linker_design.reachability import compute_reachability

def main():
    print("Running sensitivity analysis...")
    tale = load_reference("3V6T", cleaned_dir="data/pdb_cleaned")
    frame = ReferenceFrame.from_tale_structure(tale)
    scissile_table = build_scissile_phosphate_table(tale, frame, bp_range=11)
    target_coords = scissile_table[("top", 4)]["coords"]

    tolerances = [8, 10, 12, 15]          # \u00C5
    pers_lengths = [4.0, 5.0, 6.0]        # \u00C5
    linker_classes = [LinkerClass.H, LinkerClass.F, LinkerClass.M]
    n_residues = 10
    
    rng = np.random.default_rng(42)
    results = []

    for lp in pers_lengths:
        import tale_linker_design.linkers as linkers_mod
        linkers_mod.LP_FLEXIBLE = lp
        for lclass in linker_classes:
            lib = build_linker_library(classes=[lclass], lengths=[n_residues])
            ensemble = lib[(lclass, n_residues)]
            
            rmap = compute_reachability(tale, ensemble, n_samples=20_000, rng=rng)
            
            for tol in tolerances:
                boot_res = rmap.bootstrap_p_reach(target_coords, tolerance_A=tol)
                p = boot_res["estimate"]
                results.append({
                    "tolerance_A": tol, 
                    "lp": lp, 
                    "class": lclass.value, 
                    "p_reach": p
                })

    df = pd.DataFrame(results)
    out_dir = Path("data")
    out_dir.mkdir(exist_ok=True)
    out_csv = out_dir / "sensitivity_analysis.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

if __name__ == "__main__":
    main()
