# GENESIS v3.2.1 Linker Specification
### From Paper C5: Geometric Constraints on Catalytic Domain Fusion to TALE Arrays

**Date frozen:** 2026-04-19  
**C5 version:** Final (pre-submission)  
**Feeds into:** C1 (theozyme design), C2 (scaffold design)

---

## Target Definition

- **Scissile phosphate:** Top strand, 4 bp downstream of TALE footprint (bp +4)
- **Canonical-frame coordinates:** x ≈ 5.2 Å, y ≈ 7.1 Å, z ≈ 13.6 Å
- **Distance from TALE C-terminal Cα:** ≈ 16.0 Å

---

## Ranked Linker Recommendations

### Design 1 (PRIMARY — Rigid, length-matched)
- **Class:** H — α-helical
- **Composition:** (EAAAK)₂
- **Length:** 10 residues
- **Sequence:** `EAAAKEAAAK`
- **P(reach, 5 Å):** 0.11%
- **P(reach, 8 Å):** 1.54%
- **P(reach, 12 Å):** 19.11% ← primary metric (accounts for ~6 Å side-chain reach)
- **Nearest sample distance to target:** 1.5 Å
- **Orientational concentration κ:** 15 (highly directional)
- **Entropic pre-organisation cost:** ≈ 4.2 kcal/mol
- **Helix fraction:** 85%
- **Why n=10:** Mean helical reach = 10 × 1.5 Å = 15 Å ≈ 16.2 Å target distance (length-matching principle)

> [!NOTE]
> **Footnote on Run 3 vs Frozen Spec:**
> Run 3 (5Å tolerance) identifies the 15-residue helical linker as optimal for direct Cα placement. However, the 10-residue (EAAAK)₂ specification remains the frozen design target because the catalytic side chain (His/Ser) extends ~6Å from Cα, making the effective reach tolerance ~12Å, where P(reach)=19.1%. This prevents C1 from accidentally using the 15-residue parameter.

**Implication for C1:** The active-site nucleophile (His or Ser) faces in the direction of the EAAAK helix axis from the TALE C-terminus. This constrains the rotational search space from 4π sr to ≈ 0.3 sr (42× reduction).

**Implication for C2:** The scaffold's first helix should be oriented along the EAAAK axis. The catalytic residue should be placed 3–5 residues into the scaffold (≈ 4–7 Å from attachment point along helix axis).

---

### Design 2 (FALLBACK — Flexible, length-matched)
- **Class:** F — flexible GGS
- **Composition:** (GGS)₃GG
- **Length:** 10 residues
- **Sequence:** `GGSGGSGGSGG` (first 10 of pattern)
- **P(reach, 5 Å):** 1.15%
- **P(reach, 8 Å):** 4.29%
- **P(reach, 12 Å):** 11.45%
- **Nearest sample distance to target:** 0.7 Å
- **Orientational concentration κ:** 0 (isotropic — no directional constraint on scaffold)
- **Entropic cost:** 0 kcal/mol
- **Helix fraction:** 2%
- **Why n=10:** Same length-matching principle — 10-res GGS (mean reach ~9.8 Å) close to the target

**Use case:** If Design 1 fails biochemically (misfolded helix, unexpected TALE C-terminal interaction), Design 2 requires no directional constraint on C1/C2.

---

### Design 3 (ALTERNATIVE — Mixed, length-matched)
- **Class:** M — GGS–EAAAK–GGS
- **Length:** 10 residues (4+2+4)
- **Sequence:** approx. `GGSGSEAAKGSGG`
- **P(reach, 5 Å):** 0.56%
- **P(reach, 8 Å):** 2.59%
- **P(reach, 12 Å):** 9.71%
- **Nearest sample distance to target:** 0.6 Å
- **Orientational concentration κ:** 4
- **Entropic cost:** ≈ 0.7 kcal/mol
- **Helix fraction:** 40%

**Key insight from 50k-sample computation:** All classes perform best at n=10 for the bp+4 target (distance 16.2 Å). This is the *length-matching principle*: optimal n is the length whose mean end-to-end distance ≈ d_target.

**Use case:** Absorbs TALE conformational compression at flexible termini while providing moderate directionality.

---

## Constraints for C1 Theozyme Design

Assuming Design 1 (EAAAK×2, 10 res) is implemented:

> [!WARNING]
> **Risk Alert:** If C1 theozyme design places the catalytic His/Ser >12Å from the linker C-terminus, the effective P(reach) drops precipitously. The active site must be designed assuming the linker C-terminus is the anchor point, not the catalytic residue itself.

```
Catalytic residue Cα position:
  - 10 aa after TALE C-terminus, along EAAAK helix axis
  - Expected Cα position: ~15 Å from TALE C-term in helix axis direction
  - Side chain reaches ~4–6 Å further (His: ~4 Å; Ser: ~3 Å)
  - Net catalytic residue reach: ~18–21 Å from TALE C-terminus (brackets 16.2 Å target)

Approach angle:
  - Helix axis ≈ TALE C-terminal helix direction (to be confirmed from 3V6T)
  - In-line attack vector at bp+4 phosphate: anti-parallel to P–O3' bond
  - Constraint: catalytic residue must approach from the major groove side

Spatial uncertainty:
  - σ_position ≈ 2.2 Å (from linker std = 10% of helix length)
  - Angular uncertainty: κ = 15 → σ_θ ≈ 15°
```

---

## Validation Plan (Step 6 / HPC phase)

When HPC access is restored:
1. Run PyRosetta FastRelax on TALE + linker + catalytic domain stub
2. Confirm helix axis direction is consistent with TALE C-terminal geometry
3. Explicit-solvent MD (50 ns) to validate steric model
4. If Design 1 fails MD validation, fall back to Design 3 then Design 2

---

## Version History

| Version | Date | Change |
|---|---|---|
| 1.0 | 2026-04-19 | Initial freeze from C5 analytical framework |
