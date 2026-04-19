"""Unit tests for reference frame computation."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import load_reference
from tale_linker_design.frames import (
    ReferenceFrame, build_scissile_phosphate_table, genesis_target_positions
)


@pytest.fixture(scope="module")
def tale_3v6t():
    return load_reference("3V6T")


@pytest.fixture(scope="module")
def frame(tale_3v6t):
    return ReferenceFrame.from_tale_structure(tale_3v6t)


@pytest.fixture(scope="module")
def scissile_table(tale_3v6t, frame):
    return build_scissile_phosphate_table(tale_3v6t, frame)


def test_origin_at_zero(frame):
    """In the canonical frame, origin should be at (0,0,0)."""
    origin_canonical = frame.to_canonical(frame.origin)
    assert np.allclose(origin_canonical, 0.0, atol=1e-10)


def test_rotation_matrix_orthogonal(frame):
    """Rotation matrix R must be orthogonal (R @ R.T ≈ I)."""
    I = frame.R @ frame.R.T
    assert np.allclose(I, np.eye(3), atol=1e-12)


def test_roundtrip_transform(frame):
    """to_pdb(to_canonical(x)) == x for arbitrary points."""
    test_pts = np.random.default_rng(0).standard_normal((20, 3)) * 30
    for pt in test_pts:
        rt = frame.to_pdb(frame.to_canonical(pt))
        assert np.allclose(rt, pt, atol=1e-10)


def test_scissile_table_has_expected_keys(scissile_table):
    """Both strands, 0–10 bp offsets should be in the table."""
    for strand in ("top", "bottom"):
        for bp in range(11):
            assert (strand, bp) in scissile_table, f"Missing ({strand}, {bp})"


def test_top_strand_phosphates_increasing_z(scissile_table):
    """Top-strand phosphates should have increasing Z as bp increases."""
    z_vals = [scissile_table[("top", bp)]["coords"][2] for bp in range(6)]
    diffs = np.diff(z_vals)
    assert np.all(diffs > 0), f"Z not monotonically increasing: {z_vals}"


def test_primary_genesis_target_reachable_distance(scissile_table):
    """bp+4 top strand should be 13–20 Å from origin (canonical frame)."""
    dist = scissile_table[("top", 4)]["distance_from_origin"]
    assert 10 < dist < 30, f"Unexpected distance {dist:.1f} Å"


def test_genesis_target_positions():
    targets = genesis_target_positions()
    assert len(targets) == 6
    strands = {t.strand for t in targets}
    assert strands == {"top", "bottom"}
