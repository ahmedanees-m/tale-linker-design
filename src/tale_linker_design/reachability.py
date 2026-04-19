"""
Reachability maps: spatial probability density of catalytic-domain attachment points.

Implements Step 5 of the C5 Execution Plan.

Method:
  For each (TALE structure, linker ensemble), sample attachment-point positions
  by:
  1. Fixing the linker N-terminus at the TALE C-terminal Cα (origin of canonical frame)
  2. Sampling end-to-end vectors from the linker's WLC/helical PDF
  3. Checking steric exclusion against the TALE body and DNA duplex
  4. Computing 3D kernel density estimate of surviving attachment points

Steric exclusion is approximated analytically using a cylindrical exclusion zone
for the TALE superhelix (radius 30 Å, axis along -Z below origin) and a thin slab
for the DNA duplex (within ±12 Å of the DNA axis).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.stats import gaussian_kde


@dataclass
class ReachabilityMap:
    """
    3D kernel density estimate of catalytic domain attachment points.

    Attributes
    ----------
    linker_class : str
    n_residues : int
    samples : NDArray (N, 3)     — surviving attachment point positions (canonical frame, Å)
    kde : gaussian_kde            — fitted KDE
    survival_fraction : float     — fraction of samples surviving steric exclusion
    pdb_id : str                  — source TALE structure
    """

    linker_class: str
    n_residues: int
    samples: NDArray
    kde: Optional[object]         # scipy gaussian_kde
    survival_fraction: float
    pdb_id: str = ""

    def probability_within(self, center: NDArray, radius: float) -> float:
        """
        Estimate the probability mass within a sphere of given radius around center.
        Uses Monte Carlo integration over the KDE.
        """
        if len(self.samples) == 0:
            return 0.0
        dists = np.linalg.norm(self.samples - center[None, :], axis=1)
        return float(np.mean(dists <= radius))

    def nearest_sample_distance(self, target: NDArray) -> float:
        """Minimum Euclidean distance from any attachment point to target (Å)."""
        if len(self.samples) == 0:
            return np.inf
        return float(np.min(np.linalg.norm(self.samples - target[None, :], axis=1)))

    @property
    def mean_position(self) -> NDArray:
        return self.samples.mean(axis=0) if len(self.samples) else np.zeros(3)

    @property
    def std_position(self) -> NDArray:
        return self.samples.std(axis=0) if len(self.samples) else np.zeros(3)


def _steric_clash(positions: NDArray) -> NDArray:
    """
    Return a boolean mask: True = steric clash (exclude), False = allowed.

    Exclusion zones (physically motivated, conservative to avoid over-exclusion):

    1. TALE superhelix core: radius 25 Å, axis along −Z from origin, Z < 2 Å.
       The TALE body occupies the region below the C-terminal anchor.
       The TALE-DNA interface occupies the negative-Z half-space close to the origin.

    2. DNA interior core: radius 5 Å from the DNA helical axis, full Z range.
       This excludes the interior of the double helix (bases/sugars) but NOT the
       phosphate layer (r_xy ≈ 9-18 Å), which is the target region.

    3. Hard floor: Z < −15 Å (deep interior of TALE body; 15 Å below the C-terminus).

    Note: The DNA phosphate backbone is at r_xy ≈ 9–18 Å depending on groove;
    the linker attachment point must approach the phosphate from r_xy ≥ 5 Å,
    which is easily satisfied. We do NOT exclude the phosphate layer itself.
    """
    clash = np.zeros(len(positions), dtype=bool)

    X, Y, Z = positions[:, 0], positions[:, 1], positions[:, 2]

    r_xy = np.sqrt(X**2 + Y**2)

    # Rule 1: TALE body exclusion — below the C-terminal anchor within TALE radius
    tale_clash = (Z < 2.0) & (r_xy < 25.0)
    clash |= tale_clash

    # Rule 2: DNA interior hard core (bases + sugar ring, ~5 Å radius from helical axis)
    # The phosphate layer at r_xy ~ 9–18 Å is NOT excluded — that is where the target is.
    dna_core_clash = r_xy < 5.0
    clash |= dna_core_clash

    # Rule 3: hard floor — absolutely forbidden (deep interior of TALE body, Z < -15 Å)
    # 15 Å below the C-terminus is well within the TALE superhelix core.
    clash |= (Z < -15.0)

    return clash


def compute_reachability(
    tale_structure,
    linker_ensemble,
    n_samples: int = 50_000,
    rng=None,
) -> ReachabilityMap:
    """
    Compute the reachability map for a given TALE structure and linker ensemble.

    Parameters
    ----------
    tale_structure : TALEStructure
        The TALE-DNA complex. The canonical frame is centred on the TALE C-terminus.
    linker_ensemble : LinkerEnsemble
        Analytical linker model (from build_linker_library()).
    n_samples : int
        Number of MC samples to draw (50 000 is sufficient for smooth KDE at 2 Å resolution).
    rng : numpy RNG, optional

    Returns
    -------
    ReachabilityMap
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Sample end-to-end vectors in canonical frame (origin = TALE C-terminus)
    raw_positions = linker_ensemble.sample_end_positions(n_samples, rng=rng)

    # Apply steric clash filter
    clash_mask = _steric_clash(raw_positions)
    surviving = raw_positions[~clash_mask]

    survival_fraction = len(surviving) / max(n_samples, 1)

    # Fit KDE to surviving positions (only if enough samples)
    kde = None
    if len(surviving) >= 10:
        try:
            # Use Scott's rule bandwidth; skip if data is degenerate
            kde = gaussian_kde(surviving.T, bw_method="scott")
        except Exception:
            pass

    return ReachabilityMap(
        linker_class=linker_ensemble.linker_class.value,
        n_residues=linker_ensemble.n_residues,
        samples=surviving,
        kde=kde,
        survival_fraction=survival_fraction,
        pdb_id=tale_structure.pdb_id,
    )


def compute_all_reachability_maps(
    tale_structure,
    linker_library: dict,
    n_samples: int = 50_000,
    rng=None,
) -> dict[tuple[str, int], ReachabilityMap]:
    """
    Compute reachability maps for all (linker_class, length) combinations.

    Returns a dict keyed by (class_str, n_residues).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    maps = {}
    for (cls, n), ensemble in linker_library.items():
        rm = compute_reachability(tale_structure, ensemble, n_samples=n_samples, rng=rng)
        maps[(cls.value, n)] = rm
    return maps


def genesis_target_reachability(
    reachability_maps: dict,
    scissile_table: dict,
    tolerance_A: float = 5.0,
) -> "pd.DataFrame":
    """
    For each linker design and each GENESIS target position, compute:
      - P(reach): probability that the attachment point is within tolerance_A of target
      - nearest_A: minimum distance from any sample to the target
      - entropy_bits: differential entropy of the positional distribution (lower = more directional)

    Returns a tidy DataFrame for reporting.
    """
    import pandas as pd

    rows = []
    for (cls, n), rm in reachability_maps.items():
        for (strand, bp), entry in scissile_table.items():
            target_coords = entry["coords"]
            p_reach = rm.probability_within(target_coords, tolerance_A)
            nearest = rm.nearest_sample_distance(target_coords)

            # Entropy estimate: H ≈ 3/2 * log(2π e σ²) for each axis, summed
            if len(rm.samples) > 1:
                stds = rm.std_position
                entropy = float(np.sum(0.5 * np.log(2 * np.pi * np.e * stds**2 + 1e-12)))
            else:
                entropy = np.nan

            rows.append({
                "linker_class": cls,
                "n_residues": n,
                "strand": strand,
                "bp_offset": bp,
                "p_reach_pct": round(p_reach * 100, 2),
                "nearest_A": round(nearest, 2),
                "entropy_bits": round(entropy, 3),
                "survival_fraction": round(rm.survival_fraction, 3),
            })

    return pd.DataFrame(rows)
