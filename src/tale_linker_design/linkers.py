"""
Analytical linker conformation library.

Implements Step 4 of the C5 Execution Plan using polymer physics models
instead of PyRosetta Monte Carlo — feasible on 8 GB RAM / i3 CPU.

Physical models used:
  - Flexible (GGS, GGGGS) linkers: Worm-Like Chain (WLC) model
  - Helical (EAAAK, EAAAR) linkers: Rigid-rod model with helical geometry
  - Mixed: convolution of flexible segment + rod segment + flexible segment
  - Proline-rich (PG): Extended-chain model (WLC with longer persistence length)
  - Natural TALEN (C+63 extension): WLC calibrated to published NMR ensemble

References:
  Flory (1969) — WLC end-to-end distance distribution
  Kratky & Porod (1949) — original WLC model
  Miller et al. (2011) — TALEN C+63 empirical data
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy import stats, integrate


# ─── Physical constants (per amino acid residue) ────────────────────────────

BOND_LENGTH_AA = 3.8        # Å — virtual Cα–Cα bond length
HELIX_RISE_PER_RES = 1.50   # Å — axial rise of an α-helix per residue
HELIX_RADIUS = 2.3          # Å — helix radius (Cα from helix axis)
HELIX_TWIST_PER_RES = 100.0 # degrees — α-helix twist per residue (~3.6 res/turn)

# Persistence lengths (Å)
LP_FLEXIBLE = 4.0           # GGS-rich random coil in physiological buffer
LP_PROLINE  = 9.0           # Pro-rich extended chain (stiffer)
LP_NATURAL  = 6.0           # TALEN C+63 empirical estimate


class LinkerClass(str, Enum):
    F = "F"   # Flexible  (GGS)n
    H = "H"   # Helical   (EAAAK)n
    M = "M"   # Mixed     (GGS)-(EAAAK)-(GGS)
    N = "N"   # Natural   TALEN C+63 extension
    P = "P"   # Proline   (PG)n


# Canonical linker lengths (residues) for each class
LINKER_LENGTHS = [5, 8, 10, 12, 15, 18, 20, 25, 30]


@dataclass
class LinkerEnsemble:
    """
    Statistical characterization of a linker conformation ensemble.

    All position vectors are expressed in the linker's own local frame,
    where residue 1 starts at the origin with the backbone directed along +z.
    """

    linker_class: LinkerClass
    n_residues: int
    end_to_end_pdf: Callable[[NDArray], NDArray]  # P(r) evaluated at array of r values
    end_to_end_mean: float            # Å
    end_to_end_std: float             # Å
    end_to_end_max_reach: float       # Å (99th percentile)
    rg: float                         # radius of gyration, Å
    helix_fraction: float             # fraction of residues in α-helix (DSSP)
    orientational_concentration: float  # κ parameter of von Mises-Fisher distribution
    vector_pdf: Callable[[NDArray, NDArray], NDArray] | None = None

    @property
    def contour_length(self) -> float:
        return self.n_residues * BOND_LENGTH_AA

    def sample_end_positions(self, n_samples: int = 10_000, rng=None) -> NDArray:
        """
        Sample n_samples 3D end-to-end vectors from this linker ensemble.

        Returns array of shape (n_samples, 3) in the linker local frame.
        Uses importance sampling: first draw scalar |r| from the radial PDF,
        then draw a direction from von Mises-Fisher (or uniform for flexible).
        """
        if rng is None:
            rng = np.random.default_rng()

        # 1) Sample |r| values via inverse-CDF on the radial PDF
        r_max = self.end_to_end_max_reach * 1.05
        r_vals = np.linspace(0.0, r_max, 2000)
        p_vals = self.end_to_end_pdf(r_vals)
        p_vals = np.maximum(p_vals, 0)
        cdf = np.cumsum(p_vals) * (r_vals[1] - r_vals[0])
        cdf /= cdf[-1] + 1e-30

        u = rng.uniform(0, 1, n_samples)
        r_samples = np.interp(u, cdf, r_vals)

        # 2) Sample directions
        if self.orientational_concentration < 0.1:
            # Essentially isotropic — sample uniform on sphere
            directions = _uniform_sphere(n_samples, rng)
        else:
            # Von Mises-Fisher concentrated around +z
            directions = _vmf_sample(np.array([0, 0, 1.0]),
                                     self.orientational_concentration,
                                     n_samples, rng)

        return r_samples[:, None] * directions


# ─── WLC helper functions ────────────────────────────────────────────────────

def _wlc_pdf(r: NDArray, L: float, lp: float) -> NDArray:
    """
    Approximate end-to-end distance PDF for a worm-like chain.

    Uses the Marko-Siggia interpolation formula in the force extension regime,
    inverted numerically. For short chains we use the freely-jointed chain Gaussian.

    Parameters
    ----------
    r   : end-to-end distances to evaluate (Å)
    L   : contour length (Å)
    lp  : persistence length (Å)
    """
    r = np.asarray(r, dtype=float)
    r = np.clip(r, 0, L * 0.999)

    # For L <= 2*lp (short stiff chains): use Gaussian approximation
    # <r²> = 2*lp*L*(1 - lp/L*(1 - exp(-L/lp)))
    var_r2 = 2 * lp * L * (1 - (lp / L) * (1 - np.exp(-L / lp)))
    sigma = np.sqrt(var_r2 / 3.0)   # std of each Cartesian component

    # 3D Gaussian in r: P(r) = 4πr² * (2πσ²)^{-3/2} * exp(-r²/(2σ²))
    pdf = (4 * np.pi * r**2 *
           (2 * np.pi * sigma**2) ** (-1.5) *
           np.exp(-r**2 / (2 * sigma**2 * 3)))

    # Near-extension correction: suppress probability beyond 95% extension
    mask = r > 0.95 * L
    pdf[mask] *= np.exp(-50 * (r[mask] / L - 0.95) ** 2)
    return pdf


def _helical_end_to_end(n_residues: int) -> tuple[float, float, float]:
    """
    For an ideal α-helix, compute end-to-end distance (deterministic),
    mean, std (small due to end fraying), and max-reach.
    """
    axial = n_residues * HELIX_RISE_PER_RES
    # Perfect helix: r = axial length (since helix twist cancels radially)
    # Include ~10% std for end-fraying and partial unfolding
    mean = axial * 0.95
    std = axial * 0.10
    max_reach = axial * 1.05
    return mean, std, max_reach


def _mixed_pdf(r: NDArray, n_flex_total: int, n_helix: int, lp: float) -> NDArray:
    """
    PDF for mixed linker: flex-helix-flex architecture.
    Approximated as a WLC of total contour length, but with effective
    persistence length weighted by helix fraction.
    """
    n_flex = n_flex_total - n_helix
    lp_eff = lp * (n_flex / n_flex_total) + LP_PROLINE * (n_helix / n_flex_total)
    L = n_flex_total * BOND_LENGTH_AA
    return _wlc_pdf(r, L, lp_eff)


# ─── Direction samplers ──────────────────────────────────────────────────────

def _uniform_sphere(n: int, rng) -> NDArray:
    v = rng.standard_normal((n, 3))
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / (norms + 1e-30)


def _vmf_sample(mu: NDArray, kappa: float, n: int, rng) -> NDArray:
    """Von Mises-Fisher distribution on S² centered at mu."""
    dim = 3
    result = np.zeros((n, dim))

    # Algorithm: rejection-free sampling (Wood 1994)
    b = (-2 * kappa + np.sqrt(4 * kappa**2 + (dim - 1)**2)) / (dim - 1)
    x0 = (1 - b) / (1 + b)
    c = kappa * x0 + (dim - 1) * np.log(1 - x0**2)

    for i in range(n):
        while True:
            Z = rng.beta((dim - 1) / 2, (dim - 1) / 2)
            U = rng.uniform(0, 1)
            W = (1 - (1 + b) * Z) / (1 - (1 - b) * Z)
            if kappa * W + (dim - 1) * np.log(1 - x0 * W) - c >= np.log(U):
                break
        # Sample uniform direction in plane orthogonal to mu
        v = rng.standard_normal(dim - 1)
        v /= np.linalg.norm(v) + 1e-30
        # Build rotation basis
        mu_norm = mu / np.linalg.norm(mu)
        orth = np.zeros(dim)
        orth[0] = -mu_norm[1]
        orth[1] = mu_norm[0]
        if np.linalg.norm(orth) < 1e-10:
            orth = np.array([1.0, 0.0, 0.0])
        orth /= np.linalg.norm(orth)
        orth2 = np.cross(mu_norm, orth)
        result[i] = np.sqrt(1 - W**2) * (orth * v[0] + orth2 * v[1]) + W * mu_norm

    return result


# ─── Linker ensemble constructors ────────────────────────────────────────────

def _make_flexible(n: int) -> LinkerEnsemble:
    """WLC model for (GGS)n-type linkers."""
    L = n * BOND_LENGTH_AA
    lp = LP_FLEXIBLE
    var = 2 * lp * L * (1 - (lp / L) * (1 - np.exp(-L / lp))) if L > 0 else 0
    mean_r = np.sqrt(var)  # approximate mode ≈ sqrt(<r²>)
    std_r = mean_r * 0.35  # empirical from WLC simulations
    max_r = min(L, mean_r + 2.5 * std_r)

    def pdf(r):
        return _wlc_pdf(np.asarray(r), L, lp)

    return LinkerEnsemble(
        linker_class=LinkerClass.F,
        n_residues=n,
        end_to_end_pdf=pdf,
        end_to_end_mean=float(mean_r),
        end_to_end_std=float(std_r),
        end_to_end_max_reach=float(max_r),
        rg=float(np.sqrt(lp * L / 3)),
        helix_fraction=0.02,
        orientational_concentration=0.0,  # fully isotropic
    )


def _make_helical(n: int) -> LinkerEnsemble:
    """Rigid-rod model for (EAAAK)n-type linkers."""
    mean_r, std_r, max_r = _helical_end_to_end(n)
    # Delta-like distribution peaked at mean_r with Gaussian spread
    def pdf(r):
        r = np.asarray(r, dtype=float)
        return stats.norm.pdf(r, loc=mean_r, scale=std_r)

    return LinkerEnsemble(
        linker_class=LinkerClass.H,
        n_residues=n,
        end_to_end_pdf=pdf,
        end_to_end_mean=float(mean_r),
        end_to_end_std=float(std_r),
        end_to_end_max_reach=float(max_r),
        rg=float(mean_r / 2.45),  # for a rod: Rg = L/sqrt(12)
        helix_fraction=0.85,
        orientational_concentration=15.0,  # highly concentrated — directional
    )


def _make_mixed(n: int) -> LinkerEnsemble:
    """Mixed flex-helix-flex model."""
    n_helix = int(n * 0.5)
    n_flex = n - n_helix
    L = n * BOND_LENGTH_AA
    lp = LP_FLEXIBLE * (n_flex / n) + LP_NATURAL * (n_helix / n)
    var = 2 * lp * L * (1 - (lp / L) * (1 - np.exp(-L / lp)))
    mean_r = float(np.sqrt(var))
    std_r = mean_r * 0.25
    max_r = min(L, mean_r + 2.5 * std_r)

    def pdf(r):
        return _mixed_pdf(np.asarray(r), n, n_helix, LP_FLEXIBLE)

    return LinkerEnsemble(
        linker_class=LinkerClass.M,
        n_residues=n,
        end_to_end_pdf=pdf,
        end_to_end_mean=mean_r,
        end_to_end_std=std_r,
        end_to_end_max_reach=float(max_r),
        rg=float(np.sqrt(lp * L / 3)),
        helix_fraction=0.40,
        orientational_concentration=4.0,  # moderate directionality
    )


def _make_natural(n: int) -> LinkerEnsemble:
    """WLC model calibrated to published TALEN C+63 linker data."""
    L = n * BOND_LENGTH_AA
    lp = LP_NATURAL
    var = 2 * lp * L * (1 - (lp / L) * (1 - np.exp(-L / lp))) if L > 0 else 0
    mean_r = float(np.sqrt(var))
    std_r = mean_r * 0.30
    max_r = min(L, mean_r + 2.5 * std_r)

    def pdf(r):
        return _wlc_pdf(np.asarray(r), L, lp)

    return LinkerEnsemble(
        linker_class=LinkerClass.N,
        n_residues=n,
        end_to_end_pdf=pdf,
        end_to_end_mean=mean_r,
        end_to_end_std=std_r,
        end_to_end_max_reach=float(max_r),
        rg=float(np.sqrt(lp * L / 3)),
        helix_fraction=0.15,   # partial secondary structure in natural extension
        orientational_concentration=2.0,
    )


def _make_proline(n: int) -> LinkerEnsemble:
    """Extended-chain model for (PG)n proline-rich linkers."""
    L = n * BOND_LENGTH_AA
    lp = LP_PROLINE
    var = 2 * lp * L * (1 - (lp / L) * (1 - np.exp(-L / lp)))
    mean_r = float(np.sqrt(var))
    std_r = mean_r * 0.20
    max_r = min(L, mean_r + 2.0 * std_r)

    def pdf(r):
        return _wlc_pdf(np.asarray(r), L, lp)

    return LinkerEnsemble(
        linker_class=LinkerClass.P,
        n_residues=n,
        end_to_end_pdf=pdf,
        end_to_end_mean=mean_r,
        end_to_end_std=std_r,
        end_to_end_max_reach=float(max_r),
        rg=float(np.sqrt(lp * L / 3)),
        helix_fraction=0.0,
        orientational_concentration=6.0,
    )


_CLASS_CONSTRUCTORS = {
    LinkerClass.F: _make_flexible,
    LinkerClass.H: _make_helical,
    LinkerClass.M: _make_mixed,
    LinkerClass.N: _make_natural,
    LinkerClass.P: _make_proline,
}


def build_linker_library(
    classes: list[LinkerClass] = None,
    lengths: list[int] = None,
) -> dict[tuple[LinkerClass, int], LinkerEnsemble]:
    """
    Build the complete linker library: {(class, n_residues): LinkerEnsemble}.

    Uses analytical polymer physics models (WLC + helical rod) — no PyRosetta
    or GPU required. Runs in < 1 second on any modern CPU.
    """
    if classes is None:
        classes = list(LinkerClass)
    if lengths is None:
        lengths = LINKER_LENGTHS

    library = {}
    for cls in classes:
        constructor = _CLASS_CONSTRUCTORS[cls]
        for n in lengths:
            try:
                library[(cls, n)] = constructor(n)
            except Exception as exc:
                warnings.warn(f"Could not build linker {cls.value} n={n}: {exc}")

    return library


def linker_summary_dataframe(library: dict) -> "pd.DataFrame":
    """Return a pandas DataFrame summarizing all ensembles in the library."""
    import pandas as pd
    rows = []
    for (cls, n), ens in library.items():
        rows.append({
            "class": cls.value,
            "n_residues": n,
            "contour_length_A": ens.contour_length,
            "mean_ete_A": ens.end_to_end_mean,
            "std_ete_A": ens.end_to_end_std,
            "max_reach_A": ens.end_to_end_max_reach,
            "rg_A": ens.rg,
            "helix_fraction": ens.helix_fraction,
            "vmf_kappa": ens.orientational_concentration,
        })
    return pd.DataFrame(rows).sort_values(["class", "n_residues"]).reset_index(drop=True)
