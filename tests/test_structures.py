"""Unit tests for structure loading, chain extraction, and annotation."""

import json
import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from tale_linker_design.structures import (
    load_reference, PRIORITY_STRUCTURES, TALEStructure, _LITERATURE_PHOSPHATES
)


@pytest.fixture(scope="module")
def tale_3v6t():
    return load_reference("3V6T")


@pytest.fixture(scope="module")
def tale_3ugm():
    return load_reference("3UGM")


# ── TALEStructure object ────────────────────────────────────────────────────

def test_load_returns_tale_structure(tale_3v6t):
    assert isinstance(tale_3v6t, TALEStructure)


def test_pdb_id_preserved(tale_3v6t):
    assert tale_3v6t.pdb_id == "3V6T"


def test_n_repeats_correct_3v6t(tale_3v6t):
    assert tale_3v6t.n_repeats == pytest.approx(11.5, abs=0.1)


def test_resolution_correct(tale_3v6t):
    assert tale_3v6t.resolution == pytest.approx(1.85, abs=0.01)


def test_c_terminus_coords_shape(tale_3v6t):
    assert tale_3v6t.c_terminus_coords.shape == (3,)


def test_dna_chain_ids_present(tale_3v6t):
    assert len(tale_3v6t.dna_chain_ids) >= 1


# ── Phosphate coordinate extraction ────────────────────────────────────────

def test_phosphate_coords_present(tale_3v6t):
    assert len(tale_3v6t.dna_phosphate_coords) > 0


def test_phosphate_at_returns_array(tale_3v6t):
    p = tale_3v6t.phosphate_at("top", 0)
    assert p is not None
    assert p.shape == (3,)


def test_phosphate_at_returns_none_for_missing(tale_3v6t):
    result = tale_3v6t.phosphate_at("top", 999)
    assert result is None


def test_both_strands_have_phosphates(tale_3v6t):
    top_keys = [(s, bp) for (s, bp) in tale_3v6t.dna_phosphate_coords if s == "top"]
    bot_keys = [(s, bp) for (s, bp) in tale_3v6t.dna_phosphate_coords if s == "bottom"]
    assert len(top_keys) >= 5
    assert len(bot_keys) >= 5


# ── Multi-structure loading ─────────────────────────────────────────────────

def test_all_priority_structures_load():
    for pdb_id in PRIORITY_STRUCTURES:
        t = load_reference(pdb_id)
        assert isinstance(t, TALEStructure), f"Failed to load {pdb_id}"
        assert t.pdb_id == pdb_id


def test_3ugm_has_more_repeats(tale_3ugm, tale_3v6t):
    assert tale_3ugm.n_repeats > tale_3v6t.n_repeats


# ── Annotation JSON serialisation ──────────────────────────────────────────

def test_annotation_json_serialisable(tale_3v6t):
    ann = tale_3v6t.to_annotation_json()
    json_str = json.dumps(ann)
    assert isinstance(json_str, str)
    assert "3V6T" in json_str


def test_annotation_json_keys(tale_3v6t):
    ann = tale_3v6t.to_annotation_json()
    for key in ("pdb_id", "tale_chain", "dna_strands", "n_repeats", "resolution"):
        assert key in ann, f"Missing key: {key}"


# ── Literature phosphate reference data ────────────────────────────────────

def test_literature_phosphates_coverage():
    ref = _LITERATURE_PHOSPHATES["3V6T"]
    for strand in ("top", "bottom"):
        for bp in range(11):
            assert (strand, bp) in ref, f"Missing ({strand}, {bp})"


def test_literature_phosphate_shapes():
    ref = _LITERATURE_PHOSPHATES["3V6T"]
    for key, coords in ref.items():
        assert coords.shape == (3,), f"Wrong shape for {key}: {coords.shape}"


def test_top_strand_z_increases_with_bp():
    ref = _LITERATURE_PHOSPHATES["3V6T"]
    z_vals = [ref[("top", bp)][2] for bp in range(8)]
    assert all(z_vals[i] < z_vals[i+1] for i in range(len(z_vals)-1)), \
        f"Z not monotonically increasing: {z_vals}"
