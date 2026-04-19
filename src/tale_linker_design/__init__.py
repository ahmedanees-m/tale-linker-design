"""
tale_linker_design — Geometric framework for TALE-fusion catalytic domain linker design.

Public API:
    load_reference(pdb_id)          → TALEStructure
    Target(strand, bp_offset)       → Target
    recommend_linkers(tale, target) → list[LinkerDesign]
    plot_reachability(tale, ...)    → matplotlib Figure
"""

from .structures import TALEStructure, load_reference
from .frames import Target, ReferenceFrame
from .linkers import LinkerClass, build_linker_library, LinkerEnsemble
from .reachability import ReachabilityMap, compute_reachability
from .design import recommend_linkers, LinkerDesign
from .visualize import plot_reachability, plot_linker_distributions

__version__ = "0.1.0"
__all__ = [
    "TALEStructure",
    "load_reference",
    "Target",
    "ReferenceFrame",
    "LinkerClass",
    "build_linker_library",
    "LinkerEnsemble",
    "ReachabilityMap",
    "compute_reachability",
    "recommend_linkers",
    "LinkerDesign",
    "plot_reachability",
    "plot_linker_distributions",
]
