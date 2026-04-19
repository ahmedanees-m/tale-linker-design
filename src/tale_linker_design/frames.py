"""
Reference frame definition and scissile-phosphate coordinate system.

Implements Step 3 of the C5 Execution Plan.

The canonical frame is defined as:
  - Origin:  Cα of the last residue of the last full TALE repeat
  - Z-axis:  DNA helical axis pointing 3' (toward fusion domain side)
  - X-axis:  Perpendicular to Z, pointing into the major groove
  - Y-axis:  Z × X (right-handed)

All distances are in Ångströms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Target:
    """
    A specific scissile phosphate position to be reached by the catalytic domain.

    Parameters
    ----------
    strand : "top" or "bottom"
        "top"  = sense strand (same as TALE footprint strand, runs 5'→3' in +Z direction)
        "bottom" = antisense strand
    bp_offset : int
        Number of base pairs 3' of the TALE footprint end (0 = immediately adjacent).
    approach_tolerance_angstrom : float
        Acceptable distance from target for a linker attachment point to be counted
        as "reaching" the target. Default 5.0 Å (catalytic domain will itself extend
        an additional ~3 Å to place the active-site residue at the phosphate).
    """

    strand: str
    bp_offset: int
    approach_tolerance_angstrom: float = 5.0

    def __post_init__(self):
        if self.strand not in ("top", "bottom"):
            raise ValueError(f"strand must be 'top' or 'bottom', got {self.strand!r}")
        if self.bp_offset < 0:
            raise ValueError("bp_offset must be >= 0")


@dataclass
class ReferenceFrame:
    """
    Canonical coordinate frame for a TALE-DNA complex.

    The frame is expressed as a rotation matrix R (3×3) and translation t (3,)
    that map PDB coordinates → canonical frame coordinates:

        x_canonical = R @ (x_pdb - t)
    """

    origin: NDArray          # (3,) Cα of last TALE repeat in PDB frame
    z_axis: NDArray          # (3,) unit vector, DNA helical axis pointing 3'
    x_axis: NDArray          # (3,) unit vector, toward major groove
    y_axis: NDArray          # (3,) = z × x
    R: NDArray               # (3,3) rotation matrix [x_axis, y_axis, z_axis] rows
    pdb_id: str = ""

    @classmethod
    def from_tale_structure(cls, tale_structure) -> "ReferenceFrame":
        """
        Compute the canonical frame from a TALEStructure.

        The Z-axis is estimated from the average displacement of successive
        phosphate atoms along the top strand; X-axis points from the DNA
        helical axis toward the major groove.
        """
        phosphates = tale_structure.dna_phosphate_coords
        origin = tale_structure.c_terminus_coords.copy()

        # Collect top-strand phosphates in order to estimate helical axis
        top_pts = []
        for bp in range(8):
            key = ("top", bp)
            if key in phosphates:
                top_pts.append(phosphates[key])

        if len(top_pts) >= 2:
            # Helical axis ≈ direction of linear trend of successive phosphates
            pts = np.array(top_pts)
            # Simple PCA: first principal component of consecutive differences
            diffs = np.diff(pts, axis=0)
            z_raw = diffs.mean(axis=0)
            z_axis = z_raw / (np.linalg.norm(z_raw) + 1e-12)
        else:
            # Fallback: canonical B-DNA helical axis is along z in original frame
            z_axis = np.array([0.0, 0.0, 1.0])

        # X-axis: from origin toward the closest top-strand phosphate,
        # then Gram-Schmidt orthogonalised against z_axis
        if ("top", 0) in phosphates:
            x_raw = phosphates[("top", 0)] - origin
        else:
            x_raw = np.array([1.0, 0.0, 0.0])
        x_raw -= np.dot(x_raw, z_axis) * z_axis
        x_axis = x_raw / (np.linalg.norm(x_raw) + 1e-12)

        y_axis = np.cross(z_axis, x_axis)
        y_axis /= np.linalg.norm(y_axis) + 1e-12

        R = np.vstack([x_axis, y_axis, z_axis])   # rows are basis vectors
        return cls(
            origin=origin,
            z_axis=z_axis,
            x_axis=x_axis,
            y_axis=y_axis,
            R=R,
            pdb_id=tale_structure.pdb_id,
        )

    def to_canonical(self, coords_pdb: NDArray) -> NDArray:
        """Transform PDB-frame coordinates to canonical frame."""
        return self.R @ (coords_pdb - self.origin)

    def to_pdb(self, coords_canonical: NDArray) -> NDArray:
        """Transform canonical-frame coordinates back to PDB frame."""
        return self.R.T @ coords_canonical + self.origin


def build_scissile_phosphate_table(
    tale_structure,
    frame: Optional[ReferenceFrame] = None,
    bp_range: int = 11,
) -> dict:
    """
    Build a lookup table of scissile-phosphate positions in the canonical frame.

    Parameters
    ----------
    tale_structure : TALEStructure
    frame : ReferenceFrame (computed if not provided)
    bp_range : int, number of bp positions (0..bp_range-1) on each strand

    Returns
    -------
    dict with keys (strand, bp_offset) mapping to:
        {
          "coords": np.ndarray(3,),   # canonical-frame position
          "attack_vector": np.ndarray(3,),  # in-line nucleophilic attack direction
          "distance_from_origin": float,   # Å
        }
    """
    if frame is None:
        frame = ReferenceFrame.from_tale_structure(tale_structure)

    table = {}
    phosphates = tale_structure.dna_phosphate_coords

    for strand in ("top", "bottom"):
        for bp in range(bp_range):
            key = (strand, bp)
            if key not in phosphates:
                continue

            coords_pdb = phosphates[key]
            coords_can = frame.to_canonical(coords_pdb)
            dist = float(np.linalg.norm(coords_can))

            # In-line attack vector: anti-parallel to P–O3' bond direction
            # Approximated from B-DNA geometry (helical twist 36°/bp, rise 3.4 Å)
            sign = +1 if strand == "top" else -1
            theta = np.radians(36 * bp)
            attack = np.array([
                sign * np.sin(theta) * 0.57,
                sign * np.cos(theta) * 0.57,
                -sign * 0.82
            ])
            attack /= np.linalg.norm(attack)

            table[key] = {
                "coords": coords_can,
                "attack_vector": attack,
                "distance_from_origin": dist,
            }

    return table


def genesis_target_positions() -> list[Target]:
    """
    Return the six GENESIS-specific target scissile phosphate positions.

    Primary: top-strand bp +4 (optimal active-site geometry for GENESIS subunit T)
    Secondary: top-strand bp +3, +5; bottom-strand bp +3, +4, +5
    """
    return [
        Target(strand="top",    bp_offset=4, approach_tolerance_angstrom=5.0),  # primary
        Target(strand="top",    bp_offset=3, approach_tolerance_angstrom=5.0),
        Target(strand="top",    bp_offset=5, approach_tolerance_angstrom=5.0),
        Target(strand="bottom", bp_offset=3, approach_tolerance_angstrom=5.0),
        Target(strand="bottom", bp_offset=4, approach_tolerance_angstrom=5.0),
        Target(strand="bottom", bp_offset=5, approach_tolerance_angstrom=5.0),
    ]
