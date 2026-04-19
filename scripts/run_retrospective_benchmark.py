import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

from tale_linker_design.structures import load_reference
from tale_linker_design.frames import ReferenceFrame, build_scissile_phosphate_table
from tale_linker_design.linkers import build_linker_library, LinkerClass
from tale_linker_design.reachability import compute_reachability

def main():
    print("Running retrospective benchmark...")
    LITERATURE_CONSTRUCTS = [
        # (paper, year, linker_class, n_residues, target_bp_offset, reported_efficiency_score, notes)
        ("Miller 2011",  2011, "N",  63, 4, 0.85, "NΔ152/C+63; ~25% human CCR5/NTF3 editing"),
        ("Cermak 2011",  2011, "N",  63, 4, 0.80, "C+63 variant"),
        ("Cermak 2011",  2011, "N",  40, 4, 0.55, "C+40 variant; narrower spacer tolerance"),
        ("Beurdeley 2013", 2013, "N", 59, 3, 0.75, "I-TevI native linker; >23 bp optimal"),
        ("Yanik 2013",   2013, "F",  12, 2, 0.70, "PvuII fusion; 34000x preference"),
        ("Sun 2014",     2014, "F",  95, 4, 0.65, "scTALEN 95-aa flexible bridge"),
        ("Maeder 2013",  2013, "F",   4, 1, 0.55, "TET1 4-aa GGGS; 21-30% demethylation"),
        ("Kim 2013",     2013, "N",  63, 4, 0.60, "Library scale; C+63 scaffold"),
        ("Guilinger 2014", 2014, "H", 16, 4, 0.30, "fCas9 XTEN; strict spacer"),
        ("Chavez 2015",  2015, "F",  20, 8, 0.70, "VPR; flexible; activation"),
        ("Yeo 2018",     2018, "F",  20, 8, 0.65, "KRAB-MeCP2; flexible; repression"),
    ]

    tale = load_reference("3V6T", cleaned_dir="data/pdb_cleaned")
    frame = ReferenceFrame.from_tale_structure(tale)
    scissile_table = build_scissile_phosphate_table(tale, frame, bp_range=11)

    rng = np.random.default_rng(42)
    results = []

    for name, year, lclass, n, bp_offset, obs_eff, notes in LITERATURE_CONSTRUCTS:
        target_coords = scissile_table[("top", bp_offset)]["coords"]
        
        lib_dict = build_linker_library(classes=[LinkerClass(lclass)], lengths=[n])
        ensemble = lib_dict[(LinkerClass(lclass), n)]
        
        rmap = compute_reachability(tale, ensemble, n_samples=50_000, rng=rng)
        
        boot_res = rmap.bootstrap_p_reach(target_coords, tolerance_A=12.0)
        p_reach = boot_res["estimate"]
        
        results.append({
            "paper": name, "year": year, "class": lclass, "n": n,
            "bp_offset": bp_offset, "observed_efficiency": obs_eff,
            "predicted_p_reach": p_reach, "notes": notes
        })

    df = pd.DataFrame(results)
    out_dir = Path("data")
    out_dir.mkdir(exist_ok=True)
    df.to_csv(out_dir / "retrospective_benchmark.csv", index=False)
    print(f"Saved: data/retrospective_benchmark.csv")

    rho, pval = spearmanr(df["observed_efficiency"], df["predicted_p_reach"])
    print(f"Spearman rho = {rho:.3f}, p = {pval:.4f}, n = {len(df)}")

    rho_baseline_length, _ = spearmanr(df["observed_efficiency"], df["n"])
    np.random.seed(42)
    rho_random = np.array([spearmanr(df["observed_efficiency"], np.random.permutation(df["predicted_p_reach"]))[0] for _ in range(10000)])
    print(f"Baseline (length only): rho = {rho_baseline_length:.3f}")
    print(f"Random baseline: mean rho = {rho_random.mean():.3f} \u00B1 {rho_random.std():.3f}")
    print(f"Model beats random: {np.mean(rho_random < rho):.1%} of permutations")

if __name__ == "__main__":
    main()
