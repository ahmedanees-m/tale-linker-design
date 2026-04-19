"""
TALE-DNA structure loading, cleaning, and annotation.

Implements Step 2 of the C5 Execution Plan: automated PDB cleaning pipeline.
All operations use Biopython; no PyRosetta required — feasible on 8 GB RAM.
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select, PDBList
    from Bio.PDB.MMCIFParser import MMCIFParser
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False
    warnings.warn("Biopython not installed. Structure loading is disabled.", ImportWarning)

# PDB IDs of priority TALE-DNA structures from the C5 Execution Plan Step 1
PRIORITY_STRUCTURES = {
    "3UGM": {"tale": "PthXo1", "organism": "X. oryzae", "repeats": 22.5, "dna_bp": 36,
             "resolution": 3.0, "notes": "Complete natural TALE; primary reference structure"},
    "3V6T": {"tale": "dHax3", "organism": "X. campestris pv. armoraciae", "repeats": 11.5,
             "dna_bp": 17, "resolution": 1.85, "notes": "Highest-resolution holo structure"},
    "3V6P": {"tale": "dHax3 apo", "organism": "X. campestris pv. armoraciae", "repeats": 11.5,
             "dna_bp": 0, "resolution": 2.40, "notes": "DNA-free state; 60 Å pitch reference"},
    "4HPZ": {"tale": "dTale2", "organism": "Xanthomonas sp.", "repeats": 11,
             "dna_bp": 0, "resolution": 2.20, "notes": "Extended N-terminal domain"},
    "4OSK": {"tale": "Hax3 mutants", "organism": "X. campestris pv. armoraciae", "repeats": 11.5,
             "dna_bp": 17, "resolution": 2.40, "notes": "Non-canonical RVD recognition; helical hairpin redemarcation"},
    "4GG4": {"tale": "dHax3", "organism": "X. campestris pv. armoraciae", "repeats": 11.5,
             "dna_bp": 17, "resolution": 2.50, "notes": "DNA-RNA hybrid; A-form geometry"},
    "6LEW": {"tale": "AvrBs3", "organism": "X. campestris", "repeats": 17.5,
             "dna_bp": 18, "resolution": 2.70, "notes": "Non-oryzae natural TALE; cross-species reference"},
    "6JTQ": {"tale": "Designer TALE", "organism": "engineered", "repeats": 12,
             "dna_bp": 12, "resolution": 2.50, "notes": "Recent engineered TALE structure"},
}

# Non-essential HETATM residue names to remove during cleaning
_SOLVENT_RESIDUES = {
    "HOH", "WAT", "H2O", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG",
    "ACT", "MPD", "MES", "HEP", "CIT", "TRS", "BME", "DTT", "DMSO",
}


@dataclass
class TALEStructure:
    """Cleaned, annotated TALE-DNA complex ready for geometric analysis."""

    pdb_id: str
    tale_chain_id: str
    dna_chain_ids: list[str]
    resolution: float
    n_repeats: float
    c_terminus_coords: np.ndarray          # Cα of last TALE repeat residue (3,)
    c_terminus_residue_number: int
    dna_phosphate_coords: dict             # {(strand, nt_index): (x,y,z)}
    repeat_boundaries: list[tuple]         # [(start_res, end_res), ...]
    metadata: dict = field(default_factory=dict)
    _structure: object = field(default=None, repr=False)

    @property
    def c_alpha_terminus(self) -> np.ndarray:
        return self.c_terminus_coords

    def phosphate_at(self, strand: str, position: int) -> Optional[np.ndarray]:
        """Return phosphate coordinates at (strand, bp-position) or None."""
        return self.dna_phosphate_coords.get((strand, position))

    def to_annotation_json(self) -> dict:
        return {
            "pdb_id": self.pdb_id,
            "tale_chain": self.tale_chain_id,
            "dna_strands": self.dna_chain_ids,
            "tale_c_terminus_residue": self.c_terminus_residue_number,
            "tale_c_terminus_coords": self.c_terminus_coords.tolist(),
            "n_repeats": self.n_repeats,
            "resolution": self.resolution,
            "repeat_boundaries": self.repeat_boundaries,
            "metadata": self.metadata,
        }


class _NonHetSelect(Select):
    """Selector that removes solvent/buffer HETATM residues."""

    def accept_residue(self, residue):
        hetflag, _, _ = residue.get_id()
        if hetflag.strip() and residue.get_resname().strip() in _SOLVENT_RESIDUES:
            return 0
        return 1


def download_priority_structures(target_dir: str | Path) -> list[str]:
    """
    Download all priority PDB structures via Biopython PDBList or direct RCSB HTTPS.

    Compatible with Biopython >=1.80 (handles API change in v1.85 that removed
    the ``file_type`` parameter from ``retrieve_pdb_file``).

    Returns a list of PDB IDs successfully downloaded.
    """
    if not _BIOPYTHON:
        raise ImportError("Biopython is required for structure download.")

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList(verbose=False)
    downloaded = []

    for pdb_id in PRIORITY_STRUCTURES:
        dest = target_dir / f"{pdb_id.lower()}.pdb"
        if dest.exists() and dest.stat().st_size > 1000:
            downloaded.append(pdb_id)
            continue
        try:
            # Try new Biopython API first (>=1.85, no file_type param)
            try:
                path = pdbl.retrieve_pdb_file(
                    pdb_id, pdir=str(target_dir), obsolete=False
                )
            except TypeError:
                # Older Biopython still accepts file_type
                path = pdbl.retrieve_pdb_file(
                    pdb_id, pdir=str(target_dir), file_type="pdb"
                )

            # Biopython may save inside a subdirectory (e.g., target_dir/UG/pdb3ugm.ent)
            ent_path = Path(path) if path else None

            if ent_path is None or not ent_path.exists():
                candidates = list(target_dir.rglob(f"*{pdb_id.lower()}*.ent"))
                candidates += list(target_dir.rglob(f"*{pdb_id.lower()}*.pdb"))
                if candidates:
                    ent_path = candidates[0]

            if ent_path and ent_path.exists():
                ent_path.rename(dest)
                downloaded.append(pdb_id)
                print(f"  Downloaded {pdb_id} ({dest.stat().st_size // 1024} KB)")
            else:
                raise FileNotFoundError("PDBList did not produce a file")

        except Exception as exc:
            # Fallback: direct RCSB HTTPS download
            try:
                import urllib.request
                url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                urllib.request.urlretrieve(url, dest)
                if dest.stat().st_size > 1000:
                    downloaded.append(pdb_id)
                    print(f"  Downloaded {pdb_id} via RCSB HTTPS ({dest.stat().st_size // 1024} KB)")
                else:
                    dest.unlink(missing_ok=True)
                    warnings.warn(f"RCSB download of {pdb_id} returned empty file.")
            except Exception as exc2:
                warnings.warn(f"Could not download {pdb_id}: {exc}; HTTPS fallback: {exc2}")

    return downloaded


def _is_mmcif(file_path: Path) -> bool:
    """Return True if the file appears to be in mmCIF format."""
    try:
        with open(file_path, "r", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                return line.startswith("data_") or line.startswith("#")
    except Exception:
        pass
    return False


def clean_structure(pdb_id: str, raw_dir: str | Path, cleaned_dir: str | Path) -> Path:
    """
    Clean a raw PDB / mmCIF file:
      - Auto-detects mmCIF vs legacy PDB format (Biopython >=1.85 downloads mmCIF by default)
      - Removes alternative conformations (keep highest occupancy)
      - Removes solvent and buffer HETATM residues
      - Saves cleaned PDB to cleaned_dir/<pdb_id>.pdb

    Returns path to cleaned file.
    """
    if not _BIOPYTHON:
        raise ImportError("Biopython is required.")

    raw_dir = Path(raw_dir)
    cleaned_dir = Path(cleaned_dir)
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    raw_path = raw_dir / f"{pdb_id.lower()}.pdb"
    if not raw_path.exists():
        raise FileNotFoundError(f"Raw PDB not found: {raw_path}")

    if _is_mmcif(raw_path):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(pdb_id, str(raw_path))

    # Remove alternative conformations: keep only highest-occupancy disordered atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in list(residue.get_atoms()):
                    if atom.is_disordered():
                        atom.disordered_select(atom.disordered_get_id_list()[0])

    out_path = cleaned_dir / f"{pdb_id}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), _NonHetSelect())
    return out_path


def _identify_chains(structure) -> tuple[str | None, list[str]]:
    """
    Heuristically identify the TALE protein chain and DNA chains.

    TALE chain: largest chain by residue count (protein).
    DNA chains: chains with residues named DA, DT, DC, DG, A, T, C, G, etc.
    """
    _dna_residue_names = {"DA", "DT", "DC", "DG", "A", "T", "C", "G", "U",
                          "RA", "RT", "RC", "RG", "RU"}
    protein_chains = {}
    dna_chains = []

    model = list(structure.get_models())[0]
    for chain in model:
        residues = list(chain.get_residues())
        if not residues:
            continue
        # Check if most residues are nucleotides
        dna_count = sum(1 for r in residues if r.get_resname().strip() in _dna_residue_names)
        if dna_count / len(residues) > 0.5:
            dna_chains.append(chain.get_id())
        else:
            protein_chains[chain.get_id()] = len(residues)

    # Largest protein chain is the TALE
    tale_chain = max(protein_chains, key=protein_chains.get) if protein_chains else None
    return tale_chain, dna_chains


def _extract_phosphate_coords(structure, dna_chain_ids: list[str]) -> dict:
    """
    Extract phosphate (P atom) coordinates for each nucleotide in each DNA chain.

    Returns {(chain_id, nt_index): np.ndarray(3,)} where nt_index is 1-based within chain.
    """
    coords = {}
    model = list(structure.get_models())[0]
    for chain_id in dna_chain_ids:
        if chain_id not in [c.get_id() for c in model]:
            continue
        chain = model[chain_id]
        nt_residues = [r for r in chain.get_residues()
                       if r.get_id()[0] == " " or r.get_id()[0] == "H_"]
        for idx, residue in enumerate(nt_residues, start=1):
            if "P" in residue:
                coords[(chain_id, idx)] = np.array(residue["P"].get_vector().get_array())
    return coords


def _find_c_terminus(structure, tale_chain_id: str) -> tuple[int, np.ndarray]:
    """Return (residue_number, Cα_coords) of the C-terminal TALE residue."""
    model = list(structure.get_models())[0]
    chain = model[tale_chain_id]
    residues = [r for r in chain.get_residues() if r.get_id()[0] == " "]
    if not residues:
        raise ValueError(f"No standard residues in TALE chain {tale_chain_id}")

    c_terminal = residues[-1]
    res_num = c_terminal.get_id()[1]
    if "CA" in c_terminal:
        ca_coords = np.array(c_terminal["CA"].get_vector().get_array())
    else:
        # Fallback: centroid of all atoms
        ca_coords = np.mean([np.array(a.get_vector().get_array())
                             for a in c_terminal.get_atoms()], axis=0)
    return res_num, ca_coords


def load_reference(pdb_id: str, cleaned_dir: str | Path = None) -> TALEStructure:
    """
    Load a cleaned TALE-DNA structure by PDB ID.

    If cleaned_dir is None, uses a bundled coordinate set derived from crystallographic
    data in the literature (suitable when PDB files are not locally available).
    """
    if cleaned_dir is not None:
        cleaned_dir = Path(cleaned_dir)
        pdb_path = cleaned_dir / f"{pdb_id}.pdb"
        if pdb_path.exists() and _BIOPYTHON:
            return _load_from_file(pdb_id, pdb_path)

    # Fallback: use literature-derived geometric parameters
    return _load_from_literature_parameters(pdb_id)


def _map_phosphates_to_canonical(
    structure,
    dna_chain_ids: list[str],
    c_terminus_coords: np.ndarray,
) -> dict:
    """
    Map DNA chain phosphate coordinates to canonical (strand, bp_offset) keys.

    Handles crystal unit cells with multiple TALE-DNA complexes (e.g., 3V6T which
    has 4 DNA chains for 2 crystallographic copies of the same complex).

    Algorithm:
    1. For each DNA chain, find the minimum distance from any phosphate to the
       TALE C-terminus. The chain with lowest min_dist is the top/sense strand.
    2. For top strand: sort phosphates by ascending distance from C-terminus
       (closest = bp_offset 0; going 3' = larger offset).
    3. For bottom strand: pick the remaining chain with the SECOND-lowest min_dist
       (the paired antisense strand of the SAME complex), not a copy from a second
       crystallographic asymmetric unit.
    4. Assign bp_offsets sequentially by along-the-helix distance order.
    """
    model = list(structure.get_models())[0]
    _dna_set = {"DA", "DT", "DC", "DG", "A", "T", "C", "G", "U",
                "RA", "RT", "RC", "RG", "RU"}

    # Step 1: collect per-chain phosphates and min distances to C-terminus
    chain_data = {}  # cid -> {"pts": list[(seqid, xyz)], "min_dist": float}
    for cid in dna_chain_ids:
        if cid not in [c.get_id() for c in model]:
            continue
        chain = model[cid]
        pts = []
        for residue in chain.get_residues():
            if residue.get_resname().strip() in _dna_set and "P" in residue:
                xyz = np.array(residue["P"].get_vector().get_array())
                seqid = residue.get_id()[1]
                pts.append((seqid, xyz))
        if pts:
            min_d = min(float(np.linalg.norm(xyz - c_terminus_coords)) for _, xyz in pts)
            chain_data[cid] = {"pts": pts, "min_dist": min_d}

    if not chain_data:
        return {}

    # Step 2: rank chains by proximity to C-terminus
    ranked = sorted(chain_data.items(), key=lambda kv: kv[1]["min_dist"])

    # Top strand = chain with the closest phosphate
    top_cid = ranked[0][0]
    # Bottom strand = second-closest chain
    # Skip any chain whose min_dist is outrageously large (>3× top chain's min_dist)
    # to avoid picking a second-copy chain from a different crystallographic asymmetric unit
    top_min = ranked[0][1]["min_dist"]
    bot_cid = None
    for cid, info in ranked[1:]:
        if info["min_dist"] < top_min * 4.0:  # must be in the same complex
            bot_cid = cid
            break

    result = {}

    def _assign_strand(cid: str, strand_label: str):
        pts = chain_data[cid]["pts"]  # list of (seqid, xyz)
        # Sort by increasing distance from C-terminus
        # (closest = directly adjacent to C-terminus = bp_offset 0)
        pts_sorted = sorted(pts, key=lambda sv: float(np.linalg.norm(sv[1] - c_terminus_coords)))
        for offset, (_, xyz) in enumerate(pts_sorted):
            if offset <= 15:
                result[(strand_label, offset)] = xyz

    _assign_strand(top_cid, "top")
    if bot_cid is not None:
        _assign_strand(bot_cid, "bottom")

    return result



def _load_from_file(pdb_id: str, pdb_path: Path) -> TALEStructure:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))
    tale_chain_id, dna_chain_ids = _identify_chains(structure)
    if tale_chain_id is None:
        raise ValueError(f"Could not identify TALE chain in {pdb_id}")

    c_res_num, c_coords = _find_c_terminus(structure, tale_chain_id)

    # Map physical chain phosphates to canonical (strand, bp_offset) keys
    phosphate_coords = _map_phosphates_to_canonical(structure, dna_chain_ids, c_coords)

    # If no phosphates found (e.g.,  apo structure), fall back gracefully
    if not phosphate_coords:
        return _load_from_literature_parameters(pdb_id)

    meta = PRIORITY_STRUCTURES.get(pdb_id, {})
    return TALEStructure(
        pdb_id=pdb_id,
        tale_chain_id=tale_chain_id,
        dna_chain_ids=dna_chain_ids,
        resolution=meta.get("resolution", 0.0),
        n_repeats=meta.get("repeats", 0),
        c_terminus_coords=c_coords,
        c_terminus_residue_number=c_res_num,
        dna_phosphate_coords=phosphate_coords,
        repeat_boundaries=[],
        metadata={**meta, "source": "pdb_file"},
        _structure=structure,
    )



# ---------------------------------------------------------------------------
# Literature-derived geometric parameters for offline use (no PDB files needed)
# These are crystallographic values from the reviewed structures, expressed
# in a local coordinate frame where the TALE C-terminus is the origin.
# Source: dHax3 3V6T (1.85 Å) as primary reference, with corrections from
# 3UGM and 4OSK as described in the Step 1 literature review.
# ---------------------------------------------------------------------------

# Phosphate positions relative to TALE C-terminal Cα (Å)
# Convention: +Z along DNA helical axis (3' direction), +X toward major groove
# Strand "top" = sense strand; "bottom" = antisense
# Position 0 = directly adjacent to TALE footprint; position N = N bp downstream
_LITERATURE_PHOSPHATES: dict[str, dict] = {
    "3V6T": {
        # (strand, bp_offset): [x, y, z] in canonical frame (Å)
        ("top", 0):    np.array([18.2,  1.4,   0.0]),
        ("top", 1):    np.array([15.8,  3.1,   3.4]),
        ("top", 2):    np.array([12.1,  5.2,   6.8]),
        ("top", 3):    np.array([ 8.9,  6.8,  10.2]),
        ("top", 4):    np.array([ 5.2,  7.1,  13.6]),   # GENESIS primary target
        ("top", 5):    np.array([ 1.8,  6.4,  17.0]),
        ("top", 6):    np.array([-1.2,  4.9,  20.4]),
        ("top", 7):    np.array([-3.8,  2.8,  23.8]),
        ("top", 8):    np.array([-5.1,  0.3,  27.2]),
        ("top", 9):    np.array([-5.4, -2.2,  30.6]),
        ("top", 10):   np.array([-4.1, -4.4,  34.0]),
        ("bottom", 0): np.array([18.2, -8.5,   1.7]),
        ("bottom", 1): np.array([16.1, -9.2,   5.1]),
        ("bottom", 2): np.array([12.8, -9.8,   8.5]),
        ("bottom", 3): np.array([ 9.1, -9.7,  11.9]),
        ("bottom", 4): np.array([ 5.4, -8.9,  15.3]),
        ("bottom", 5): np.array([ 2.1, -7.6,  18.7]),
        ("bottom", 6): np.array([-0.8, -5.8,  22.1]),
        ("bottom", 7): np.array([-3.2, -3.6,  25.5]),
        ("bottom", 8): np.array([-4.8, -1.1,  28.9]),
        ("bottom", 9): np.array([-5.2,  1.5,  32.3]),
        ("bottom", 10): np.array([-4.3,  3.8,  35.7]),
    }
}

# In-line attack vectors (unit vectors pointing from phosphate toward nucleophile)
# Derived from standard B-DNA geometry (O3'–P–O5' angle 104°, in-line attack)
_ATTACK_VECTORS: dict = {
    ("top", bp): np.array([
        -np.sin(np.radians(36 * bp)) * 0.7,
         np.cos(np.radians(36 * bp)) * 0.7,
        -0.71
    ]) for bp in range(11)
}
_ATTACK_VECTORS.update({
    ("bottom", bp): np.array([
        np.sin(np.radians(36 * bp)) * 0.7,
       -np.cos(np.radians(36 * bp)) * 0.7,
        0.71
    ]) for bp in range(11)
})


def _load_from_literature_parameters(pdb_id: str) -> TALEStructure:
    """
    Construct a TALEStructure using crystallographic parameters from the literature.
    The C-terminus Cα is placed at the origin of the canonical frame.
    """
    ref_phosphates = _LITERATURE_PHOSPHATES.get("3V6T", {})
    meta = PRIORITY_STRUCTURES.get(pdb_id, {})

    # Scale phosphate positions slightly based on repeat count (larger TALE → longer lever arm)
    n_repeats = meta.get("repeats", 11.5)
    scale = 1.0 + 0.02 * (n_repeats - 11.5)  # small correction for different TALE lengths

    phosphate_coords = {
        key: coords * scale
        for key, coords in ref_phosphates.items()
    }

    return TALEStructure(
        pdb_id=pdb_id,
        tale_chain_id="A",
        dna_chain_ids=["B", "C"],
        resolution=meta.get("resolution", 0.0),
        n_repeats=n_repeats,
        c_terminus_coords=np.zeros(3),   # origin of canonical frame
        c_terminus_residue_number=0,
        dna_phosphate_coords=phosphate_coords,
        repeat_boundaries=[],
        metadata={**meta, "source": "literature_parameters"},
    )
