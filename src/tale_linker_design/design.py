"""
Inverse linker design: from target position to linker specification.

Implements Step 7 of the C5 Execution Plan.

Given a target scissile phosphate position and reachability maps for all
linker classes and lengths, rank candidate linker designs by:
  1. P(reach) — probability of placing attachment point within tolerance of target
  2. Orientational concentration κ — higher = more directional control
  3. Length — shorter is preferred (lower entropic cost, simpler expression)

Then produce the top-3 GENESIS linker recommendations:
  Design 1 (Rigid):   α-helical (EAAAK-based)
  Design 2 (Flexible): GGS-based
  Design 3 (Mixed):   GGS–EAAAK–GGS
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .frames import Target, genesis_target_positions
from .linkers import LinkerClass, LinkerEnsemble


# Published TALEN-fusion benchmarks from the Step 1 literature review
# Format: {name: {linker_class, n_residues, reported_efficiency_pct, cut_site_bp_offset}}
PUBLISHED_FUSIONS = {
    "TALEN_FokI_C+63_Miller2011": {
        "linker_class": "N", "n_residues": 63,
        "reported_efficiency_pct": 85.0, "cut_site_bp_offset": 10,
        "strand": "top", "citation": "Miller et al. 2011, Nat Biotechnol 29:143",
    },
    "TALEN_FokI_NΔ153_C+40": {
        "linker_class": "N", "n_residues": 40,
        "reported_efficiency_pct": 45.0, "cut_site_bp_offset": 8,
        "strand": "top", "citation": "Miller et al. 2011",
    },
    "cTALEN_I-TevI_Beurdeley2013": {
        "linker_class": "N", "n_residues": 59,
        "reported_efficiency_pct": 70.0, "cut_site_bp_offset": 15,
        "strand": "top", "citation": "Beurdeley et al. 2013, Nat Commun 4:1762",
    },
    "scTALEN_FokI_tandem_Sun2014": {
        "linker_class": "F", "n_residues": 95,
        "reported_efficiency_pct": 60.0, "cut_site_bp_offset": 10,
        "strand": "top", "citation": "Sun et al. 2014",
    },
    "TALE_PvuII_Yanik2013": {
        "linker_class": "M", "n_residues": 10,
        "reported_efficiency_pct": 34000.0,  # fold-selectivity, not efficiency %
        "cut_site_bp_offset": 6,
        "strand": "top", "citation": "Yanik et al. 2013",
    },
    "DdCBE_TALE_base_editor_2020": {
        "linker_class": "F", "n_residues": 16,
        "reported_efficiency_pct": 50.0, "cut_site_bp_offset": 4,
        "strand": "top", "citation": "Mok et al. 2020, Nature 583:631",
    },
    "TALE_KRAB_repressor": {
        "linker_class": "F", "n_residues": 10,
        "reported_efficiency_pct": 85.0, "cut_site_bp_offset": 0,
        "strand": "top", "citation": "Garg et al. 2012",
    },
    "TALE_VP64_activator": {
        "linker_class": "F", "n_residues": 8,
        "reported_efficiency_pct": 75.0, "cut_site_bp_offset": 0,
        "strand": "top", "citation": "Zhang et al. 2011",
    },
}


@dataclass
class LinkerDesign:
    """A single linker design recommendation with supporting metrics."""

    name: str
    linker_class: LinkerClass
    n_residues: int
    sequence_motif: str
    p_reach_primary_pct: float        # probability reaching primary GENESIS target
    nearest_A: float                   # minimum approach distance to target
    entropic_cost_kcal_mol: float      # T·ΔS for pre-organization (estimated)
    orientational_concentration: float # κ of the vMF distribution
    notes: str = ""

    @property
    def recommended_sequence(self) -> str:
        """Generate a concrete linker sequence for this design."""
        motifs = {
            LinkerClass.F: "GGS",
            LinkerClass.H: "EAAAK",
            LinkerClass.M: None,
            LinkerClass.N: "GSSG",
            LinkerClass.P: "PG",
        }
        if self.linker_class == LinkerClass.M:
            n_flex = max(2, self.n_residues // 4)
            n_helix = self.n_residues - 2 * n_flex
            flex_part = ("GGS" * 10)[:n_flex]
            helix_part = ("EAAAK" * 10)[:n_helix]
            return flex_part + helix_part + flex_part
        base = motifs.get(self.linker_class, "GGS")
        return (base * 20)[:self.n_residues]


def recommend_linkers(
    tale_structure,
    target: Target,
    reachability_maps: dict,
    linker_library: dict,
    scissile_table: dict,
    top_k: int = 3,
    min_p_reach: float = 5.0,
) -> list[LinkerDesign]:
    """
    Recommend the top-k linker designs for a given TALE structure and target.

    Parameters
    ----------
    tale_structure : TALEStructure
    target : Target — the desired scissile phosphate
    reachability_maps : dict keyed by (class_str, n_residues) → ReachabilityMap
    linker_library : dict keyed by (LinkerClass, n_residues) → LinkerEnsemble
    scissile_table : dict from build_scissile_phosphate_table()
    top_k : int — number of recommendations
    min_p_reach : float — minimum acceptable P(reach) in percent

    Returns
    -------
    list of LinkerDesign, length top_k
    """
    target_key = (target.strand, target.bp_offset)
    if target_key not in scissile_table:
        raise ValueError(f"Target {target_key} not in scissile table. "
                         f"Available: {list(scissile_table.keys())}")

    target_coords = scissile_table[target_key]["coords"]

    candidates = []
    for (cls_str, n), rm in reachability_maps.items():
        p_reach = rm.probability_within(target_coords, target.approach_tolerance_angstrom) * 100
        nearest = rm.nearest_sample_distance(target_coords)

        # Look up orientational concentration from library
        cls = LinkerClass(cls_str)
        ens = linker_library.get((cls, n))
        kappa = ens.orientational_concentration if ens else 0.0

        # Entropic cost estimate: T·ΔS ~ kT * ln(Ω_restricted / Ω_free)
        # Approximate: 0.3 kcal/mol per degree of freedom frozen
        n_frozen_dof = max(0, (kappa - 0.5) * 2)
        entropic_cost = n_frozen_dof * 0.30

        if p_reach >= min_p_reach or nearest < 15.0:
            candidates.append({
                "cls": cls, "n": n, "p_reach": p_reach,
                "nearest": nearest, "kappa": kappa, "entropic": entropic_cost,
            })

    if not candidates:
        return []

    df = pd.DataFrame(candidates)
    # Score: maximize P(reach), then maximize kappa, then minimize n
    df["score"] = (
        df["p_reach"] * 1.0
        + df["kappa"] * 2.0
        - df["n"] * 0.1
    )
    df = df.sort_values("score", ascending=False).reset_index(drop=True)

    # Ensure diversity: pick best F, H, M as the three designs
    designs = []
    used_classes = set()
    preferred_order = [LinkerClass.H, LinkerClass.F, LinkerClass.M,
                       LinkerClass.N, LinkerClass.P]

    for cls_pref in preferred_order:
        sub = df[df["cls"] == cls_pref]
        if sub.empty:
            continue
        best = sub.iloc[0]
        used_classes.add(cls_pref)
        designs.append(LinkerDesign(
            name=f"Design {'RGB'[len(designs)]} ({cls_pref.value}-class, {best['n']} res)",
            linker_class=cls_pref,
            n_residues=int(best["n"]),
            sequence_motif={
                LinkerClass.F: "(GGS)n",
                LinkerClass.H: "(EAAAK)n",
                LinkerClass.M: "(GGS)-(EAAAK)-(GGS)",
                LinkerClass.N: "TALEN C+63 extension",
                LinkerClass.P: "(PG)n",
            }[cls_pref],
            p_reach_primary_pct=float(best["p_reach"]),
            nearest_A=float(best["nearest"]),
            entropic_cost_kcal_mol=float(best["entropic"]),
            orientational_concentration=float(best["kappa"]),
        ))
        if len(designs) == top_k:
            break

    return designs


def genesis_linker_specification(designs: list[LinkerDesign], target: Target) -> dict:
    """
    Format the GENESIS linker specification document (Step 7d of Execution Plan).
    Returns a dict suitable for JSON serialisation.
    """
    return {
        "genesis_version": "v3.2.1",
        "paper": "C5",
        "target_scissile_phosphate": {
            "strand": target.strand,
            "bp_offset": target.bp_offset,
            "approach_tolerance_A": target.approach_tolerance_angstrom,
        },
        "candidate_linkers": [
            {
                "rank": i + 1,
                "name": d.name,
                "linker_class": d.linker_class.value,
                "n_residues": d.n_residues,
                "sequence_motif": d.sequence_motif,
                "recommended_sequence": d.recommended_sequence,
                "p_reach_primary_pct": d.p_reach_primary_pct,
                "nearest_distance_A": d.nearest_A,
                "entropic_cost_kcal_mol": d.entropic_cost_kcal_mol,
                "orientational_concentration_kappa": d.orientational_concentration,
                "notes": d.notes,
            }
            for i, d in enumerate(designs)
        ],
        "frozen_for_C1_C2": {
            "recommended_linker_rank": 1,
            "description": designs[0].name if designs else "None",
            "sequence": designs[0].recommended_sequence if designs else "",
            "constraints": {
                "max_length_residues": designs[0].n_residues if designs else None,
                "class": designs[0].linker_class.value if designs else None,
                "p_reach_pct_at_bp4_top": designs[0].p_reach_primary_pct if designs else 0,
            }
        }
    }
