import numpy as np
import pytest

from tale_linker_design.linkers import build_linker_library, LinkerClass
from tale_linker_design.reachability import compute_all_reachability_maps

# Mock tale structure just to have a target vector
class MockTALE:
    pdb_id = "MOCK"
    c_terminus_coords = np.array([0.0, 0.0, 0.0])
    dna_phosphate_coords = {("top", 4): np.array([10.0, 0.0, 10.0])}

    def get_steric_centers(self):
        return np.array([[5.0, 5.0, 5.0]])

def test_reachability_determinism():
    """Verify that multiple runs with the same seed produce bitwise-identical results."""
    tale = MockTALE()
    library = build_linker_library(classes=[LinkerClass.H], lengths=[10, 15])
    
    rng1 = np.random.default_rng(42)
    maps1 = compute_all_reachability_maps(tale, library, n_samples=5000, rng=rng1)
    
    rng2 = np.random.default_rng(42)
    maps2 = compute_all_reachability_maps(tale, library, n_samples=5000, rng=rng2)

    for key in maps1:
        np.testing.assert_array_equal(maps1[key].samples, maps2[key].samples)
