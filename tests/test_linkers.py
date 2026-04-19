"""Unit tests for the linker library (analytical WLC and helical models)."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.linkers import (
    build_linker_library, LinkerClass, BOND_LENGTH_AA, linker_summary_dataframe
)


@pytest.fixture(scope="module")
def library():
    return build_linker_library()


def test_library_size(library):
    """Library must have 5 classes × 9 lengths = 45 ensembles."""
    assert len(library) == 5 * 9


def test_flexible_mean_consistent_with_wlc(library):
    """Flexible linker end-to-end mean should scale roughly as sqrt(N)."""
    ens_10 = library[(LinkerClass.F, 10)]
    ens_20 = library[(LinkerClass.F, 20)]
    ratio = ens_20.end_to_end_mean / ens_10.end_to_end_mean
    # For WLC: <r²> ∝ N for long chains → ratio ≈ sqrt(2) ≈ 1.41
    assert 1.1 < ratio < 2.0, f"Unexpected ratio {ratio}"


def test_helical_mean_scales_linearly(library):
    """Helical linker end-to-end mean should scale linearly with N (rigid rod)."""
    ens_10 = library[(LinkerClass.H, 10)]
    ens_20 = library[(LinkerClass.H, 20)]
    ratio = ens_20.end_to_end_mean / ens_10.end_to_end_mean
    assert 1.7 < ratio < 2.3, f"Expected ~2.0, got {ratio}"


def test_max_reach_never_exceeds_contour(library):
    """Max reach must be ≤ contour length for all ensembles."""
    for (cls, n), ens in library.items():
        L = n * BOND_LENGTH_AA
        assert ens.end_to_end_max_reach <= L * 1.01, (
            f"{cls.value} n={n}: max_reach {ens.end_to_end_max_reach:.1f} > L={L:.1f}"
        )


def test_helical_helix_fraction(library):
    """Helical linkers should have helix_fraction > 0.6."""
    for n in [10, 15, 20]:
        ens = library[(LinkerClass.H, n)]
        assert ens.helix_fraction > 0.6, f"n={n}: helix_fraction={ens.helix_fraction}"


def test_flexible_helix_fraction(library):
    """Flexible linkers should have low helix_fraction (< 0.1)."""
    for n in [10, 15, 20]:
        ens = library[(LinkerClass.F, n)]
        assert ens.helix_fraction < 0.1


def test_pdf_non_negative(library):
    """All PDFs must be non-negative."""
    r = np.linspace(0, 100, 200)
    for (cls, n), ens in library.items():
        vals = ens.end_to_end_pdf(r)
        assert np.all(vals >= -1e-10), f"{cls.value} n={n}: negative PDF values"


def test_sampling_shape(library):
    """sample_end_positions should return (n_samples, 3) array."""
    ens = library[(LinkerClass.F, 15)]
    samples = ens.sample_end_positions(500)
    assert samples.shape == (500, 3)


def test_sampling_within_max_reach(library):
    """Sampled distances should not exceed contour length."""
    ens = library[(LinkerClass.F, 15)]
    samples = ens.sample_end_positions(2000)
    dists = np.linalg.norm(samples, axis=1)
    L = ens.contour_length
    assert dists.max() <= L * 1.05, f"Sample distance {dists.max():.1f} exceeds L={L:.1f}"


def test_summary_dataframe(library):
    """Summary DataFrame should have correct shape."""
    df = linker_summary_dataframe(library)
    assert len(df) == len(library)
    assert "mean_ete_A" in df.columns
