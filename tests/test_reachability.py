"""Unit tests for reachability map computation."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import load_reference
from tale_linker_design.frames import ReferenceFrame, build_scissile_phosphate_table
from tale_linker_design.linkers import build_linker_library, LinkerClass
from tale_linker_design.reachability import compute_reachability, compute_all_reachability_maps


@pytest.fixture(scope="module")
def tale():
    return load_reference("3V6T")


@pytest.fixture(scope="module")
def library():
    return build_linker_library()


@pytest.fixture(scope="module")
def rm_flex_20(tale, library):
    ens = library[(LinkerClass.F, 20)]
    return compute_reachability(tale, ens, n_samples=5_000)


def test_samples_non_empty(rm_flex_20):
    assert len(rm_flex_20.samples) > 0


def test_survival_fraction_positive(rm_flex_20):
    assert 0 < rm_flex_20.survival_fraction <= 1.0


def test_samples_shape(rm_flex_20):
    assert rm_flex_20.samples.ndim == 2
    assert rm_flex_20.samples.shape[1] == 3


def test_no_samples_inside_tale_body(rm_flex_20):
    """Steric filter should have removed all points deep inside the TALE (Z < -15 Å)."""
    hard_floor_violations = rm_flex_20.samples[:, 2] < -15.0
    assert not np.any(hard_floor_violations)


def test_probability_within_reasonable(rm_flex_20):
    """P(within 30 Å of origin) should be meaningful (>10%) for a 20-res flexible linker."""
    p = rm_flex_20.probability_within(np.zeros(3), radius=30.0)
    assert p > 0.05, f"Unexpectedly low reach probability: {p:.3f}"


def test_nearest_distance_non_negative(rm_flex_20):
    target = np.array([5.2, 7.1, 13.6])  # bp+4 top strand coords
    dist = rm_flex_20.nearest_sample_distance(target)
    assert dist >= 0


def test_full_library_maps(tale, library):
    """Compute maps for a subset; verify dict keys are correct."""
    subset = {k: v for k, v in library.items() if k[1] == 10}
    maps = compute_all_reachability_maps(tale, subset, n_samples=1_000)
    for cls in LinkerClass:
        assert (cls.value, 10) in maps
