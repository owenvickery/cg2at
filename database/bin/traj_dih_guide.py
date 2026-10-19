"""
traj_dih_guide.py
=================
Generates trajectory-informed dihedral restraints for CG2AT.

Extracts CG-backbone-derived φ/ψ dihedral distributions from a new
Martini 3 CG trajectory (which has explicit N, CA, C, O backbone beads),
computes per-dihedral Trajectory Dihedral Confidence (TDC) scores, and
produces GROMACS dihedral restraints that improve Ramachandran quality
after CG→AA backmapping — without biasing toward any reference structure.

IMPORTANT: This script requires the NEW Martini 3 four-bead backbone
representation (N, CA, C, O interaction sites per residue).
Old Martini 3 single-BB-bead trajectories are NOT supported.

Note on coordinate precision: Although N, CA, C, O are named as backbone
atoms, they are CG interaction sites — not exact atomistic coordinates.
The derived φ/ψ angles are therefore described as "CG-backbone-derived"
or "four-bead-backbone-derived" restraints, not fully atomistic φ/ψ.

Workflow
--------
    # Step 1 — run CG2AT de novo (no -a flag)
    cg2at.py -c last_frame.gro -fg martini_3-2_nextgen_protein_charmm36 \\
             -ff charmm36-cg2at -w tip3p

    # Step 2 — generate dihedral restraints
    python traj_dih_guide.py --restraints \\
        -tpr  protein_CG.tpr \\
        -xtc  trajectory.xtc \\
        -aa   CG2AT_*/FINAL/final_cg2at_de_novo.pdb \\
        -o    dihedral_restraints.itp \\
        -step 10

    # Step 3 — add to PROTEIN_0.itp (instructions printed automatically)
    # Step 4 — run short NVT MD with restraints (instructions printed)

    # Optional: generate backbone reference PDB + TDC plots
    python traj_dih_guide.py \\
        -tpr  protein_CG.tpr \\
        -xtc  trajectory.xtc \\
        -c    last_frame.gro \\
        -o    traj_backbone_ref.pdb \\
        -step 10 -prefix myprotein

Scientific basis
----------------
TDC (Trajectory Dihedral Confidence) — per dihedral
    Analogous to pLDDT but measuring dihedral consistency.
    Computed independently for φ and ψ:

        TDC_phi[i] = 100 × (1 - σ_phi / σ_max)
        TDC_psi[i] = 100 × (1 - σ_psi / σ_max)

    where σ_max = 82.7° = std of uniform circular distribution.

    This ensures stable φ can be restrained independently of flexible ψ
    (and vice versa), avoiding over-restraining mixed-flexibility residues.

Force constant — per dihedral
    k_phi = k_max × exp(-σ_phi / σ_inflection)
    k_psi = k_max × exp(-σ_psi / σ_inflection)

    Default k_max = 500 kJ/mol/rad² (conservative — avoids over-restraining
    when many φ/ψ restraints are simultaneously active).
    Use --kmax 1000 for stronger restraints if needed.

Multimodality filter
    Restraints are skipped if the dihedral distribution is bimodal.
    Circular mean of a bimodal distribution is meaningless and can
    force the structure into an incorrect conformation.
    A restraint is written only if:
        - N samples >= MIN_SAMPLES_UNIMODAL (default 30)
        - circular std <= MAX_STD_UNIMODAL (default 35°)
        - dominant basin fraction >= MIN_BASIN_FRACTION (default 0.70)

Adaptive tolerance (dphi) — per dihedral
    dphi = clip(σ × DPHI_SCALE, DPHI_MIN, DPHI_MAX)

    Rigid dihedral  (σ=10°) → dphi=5.0°  (tight)
    Medium          (σ=25°) → dphi=7.5°
    Flexible        (σ=50°) → dphi=15.0° (loose)

    Override with --dphi flag for fixed tolerance (e.g. --dphi 10.0)

Residue mapping
    Restraint atom indices are derived from the AA PDB structure.
    Residues are matched by (resname, chain, sequential position)
    for robustness across multi-chain systems and non-standard numbering.

Requirements
    MDAnalysis >= 2.0
    numpy, scipy, matplotlib

Developer : Hafez Razmazma
Contact: hafez.razmazma@warwick.ac.uk
"""

import os
import sys
import warnings
import argparse
import numpy as np
from pathlib import Path

# ── Matplotlib backend — must be set BEFORE importing pyplot ──────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
except Exception:
    pass
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
from scipy.stats import circmean, circstd

# ══════════════════════════════════════════════════════════════════════════════
# Constants — adjust here to tune behaviour
# ══════════════════════════════════════════════════════════════════════════════
MAX_CIRCULAR_STD     = 82.7   # degrees — std of uniform circular distribution
TDC_HIGH             = 70.0   # TDC > 70  → high confidence
TDC_LOW              = 40.0   # TDC < 40  → low confidence (dihedral skipped)
SIG_INFLECTION       = 30.0   # degrees — force constant decay scale
SIG_SCALE            = 10.0   # degrees — sigmoid width for weight
DPHI_MIN             = 5.0    # degrees — minimum adaptive dphi
DPHI_MAX             = 20.0   # degrees — maximum adaptive dphi
DPHI_SCALE           = 0.3    # fraction of std used for dphi
K_MAX                = 500.0  # kJ/mol/rad² — conservative default (was 1000)
                               # Use --kmax 1000 for stronger restraints

# Multimodality filter parameters
# These are applied AFTER TDC filtering — only for high/medium confidence
# residues that might still have bimodal distributions.
#
# IMPORTANT — CG-appropriate thresholds:
#   CG trajectories naturally show broader distributions than AA.
#   MAX_STD_UNIMODAL is intentionally set higher than for AA data.
#   The primary quality gate is TDC_LOW — the multimodality filter
#   is a secondary check for genuine two-basin distributions only.
MIN_SAMPLES_UNIMODAL = 30     # minimum frames to attempt unimodality check
MAX_STD_UNIMODAL     = 60.0   # degrees — max std to consider unimodal
                               # (set high for CG data; AA would use ~35°)
MIN_BASIN_FRACTION   = 0.70   # fraction of samples in dominant basin

# Required backbone atom names for new Martini 3 four-bead backbone
REQUIRED_BACKBONE    = ["N", "CA", "C"]


# ══════════════════════════════════════════════════════════════════════════════
# 0.  BACKBONE VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def validate_backbone(tpr_file, xtc_file):
    """
    Verify the trajectory uses the new Martini 3 four-bead backbone.
    Exits with a clear error if old BB-bead representation is detected.

    This script requires: N, CA, C (and optionally O) per residue.
    Old Martini 3 with a single BB bead per residue is NOT supported.
    """
    u       = mda.Universe(tpr_file, xtc_file)
    protein = u.select_atoms("protein")

    if len(protein) == 0:
        sys.exit("ERROR: No protein atoms found in trajectory.")

    bead_names = set(protein.names)

    # Check for old Martini 3 single-BB backbone
    if 'BB' in bead_names and 'CA' not in bead_names:
        sys.exit(
            "\nERROR: Old Martini 3 single-BB backbone detected.\n"
            "\n"
            "This script requires the NEW Martini 3 four-bead backbone\n"
            "representation with explicit N, CA, C, O interaction sites.\n"
            "\n"
            "Detected bead names: " + str(sorted(bead_names)) + "\n"
            "\n"
            "If you are using old Martini 3 (BB bead), please:\n"
            "  1. Re-simulate with the new Martini 3 force field, OR\n"
            "  2. Use CG2AT de novo mode without TDC restraints.\n"
        )

    # Check that required backbone atoms exist
    for atom_name in REQUIRED_BACKBONE:
        if len(protein.select_atoms(f"name {atom_name}")) == 0:
            sys.exit(
                f"\nERROR: Required backbone site '{atom_name}' not found.\n"
                "\n"
                "traj_dih_guide.py requires the new Martini 3 four-bead\n"
                "backbone representation with N, CA, C, O interaction sites.\n"
                "\n"
                f"Found bead names: {sorted(bead_names)}\n"
                "\n"
                "Old Martini 3 (BB bead only) is not supported.\n"
            )

    print(f"  Backbone validation: four-bead representation confirmed ✅")
    print(f"  Backbone sites found: "
          f"{sorted(b for b in bead_names if b in ['N','CA','C','O'])}")


# ══════════════════════════════════════════════════════════════════════════════
# 1.  PHI / PSI EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def extract_phi_psi(tpr_file, xtc_file, start=0, stop=None, step=1):
    """
    Extract CG-backbone-derived φ/ψ distributions from a new Martini 3
    CG trajectory with explicit N, CA, C, O backbone interaction sites.

    NOTE: φ/ψ angles are computed from CG bead coordinates (N, CA, C, O),
    not from exact atomistic positions. These are described as
    "four-bead-backbone-derived" or "CG-backbone-derived" dihedrals.

    phi = C(i-1) - N(i)  - CA(i) - C(i)
    psi = N(i)   - CA(i) - C(i)  - N(i+1)

    Returns
    -------
    distributions : dict
        {resid: {'resname': str, 'phi': np.ndarray, 'psi': np.ndarray}}
    """
    u = mda.Universe(tpr_file, xtc_file)
    protein  = u.select_atoms("protein")
    residues = protein.residues
    n_res    = len(residues)

    # Check for old BB-bead model
    bead_names = set(protein.names)
    if 'BB' in bead_names and 'CA' not in bead_names:
        sys.exit(
            "\nERROR: Old Martini 3 single-BB backbone detected.\n"
            "This script requires the new Martini 3 four-bead backbone "
            "(N, CA, C, O).\n"
        )

    distributions = {
        res.resid: {'resname': res.resname, 'phi': [], 'psi': []}
        for res in residues
    }

    n_proc = len(u.trajectory[start:stop:step])
    print(f"Found {n_res} residues in trajectory")
    print(f"Trajectory frames: {len(u.trajectory)}")
    print(f"Using every {step} frame(s) — processing {n_proc} frames")
    print(f"Calculating phi/psi angles ...")

    frame_count = 0
    for ts in u.trajectory[start:stop:step]:
        for i, res in enumerate(residues):

            # ── PHI ──────────────────────────────────────────────────────
            if i > 0:
                prev = residues[i - 1]
                try:
                    C_prev  = prev.atoms.select_atoms("name C")[0].position
                    N_curr  = res.atoms.select_atoms("name N")[0].position
                    CA_curr = res.atoms.select_atoms("name CA")[0].position
                    C_curr  = res.atoms.select_atoms("name C")[0].position
                    phi = calc_dihedrals(
                        C_prev[np.newaxis],  N_curr[np.newaxis],
                        CA_curr[np.newaxis], C_curr[np.newaxis]
                    )[0]
                    distributions[res.resid]['phi'].append(np.degrees(phi))
                except (IndexError, Exception):
                    pass

            # ── PSI ──────────────────────────────────────────────────────
            if i < n_res - 1:
                nxt = residues[i + 1]
                try:
                    N_curr  = res.atoms.select_atoms("name N")[0].position
                    CA_curr = res.atoms.select_atoms("name CA")[0].position
                    C_curr  = res.atoms.select_atoms("name C")[0].position
                    N_next  = nxt.atoms.select_atoms("name N")[0].position
                    psi = calc_dihedrals(
                        N_curr[np.newaxis],  CA_curr[np.newaxis],
                        C_curr[np.newaxis],  N_next[np.newaxis]
                    )[0]
                    distributions[res.resid]['psi'].append(np.degrees(psi))
                except (IndexError, Exception):
                    pass

        frame_count += 1
        if frame_count % 100 == 0:
            print(f"  Processed {frame_count} frames ...")

    for resid in distributions:
        distributions[resid]['phi'] = np.array(distributions[resid]['phi'])
        distributions[resid]['psi'] = np.array(distributions[resid]['psi'])

    print(f"Done — {frame_count} frames processed for {n_res} residues.\n")
    return distributions


# ══════════════════════════════════════════════════════════════════════════════
# 2.  MULTIMODALITY FILTER
# ══════════════════════════════════════════════════════════════════════════════

def _is_unimodal(values_deg,
                 min_samples=MIN_SAMPLES_UNIMODAL,
                 max_std=MAX_STD_UNIMODAL,
                 dominant_fraction=MIN_BASIN_FRACTION):
    """
    Check if a circular dihedral distribution is unimodal (single basin).

    Returns True if the distribution is suitable for restraining:
        - enough samples
        - circular std below threshold (not hopelessly broad)
        - dominant basin contains sufficient fraction of samples

    The multimodality check uses a two-step approach:
    1. Quick std check — if std <= max_std, likely unimodal
    2. Basin fraction check — what fraction of samples fall within
       ±60° of the circular mean? If < dominant_fraction → bimodal.

    Parameters
    ----------
    values_deg       : array of dihedral values in degrees
    min_samples      : minimum sample count (default 30)
    max_std          : maximum circular std to attempt basin check (default 60°)
                       NOTE: set higher for CG data than for AA data
    dominant_fraction: minimum fraction in dominant basin (default 0.70)

    Returns
    -------
    bool : True if distribution appears unimodal
    str  : reason string if rejected (or '' if accepted)
    """
    n = len(values_deg)

    if n < min_samples:
        return False, f"insufficient samples ({n} < {min_samples})"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rad  = np.radians(values_deg)
        std  = np.degrees(circstd(rad))
        mean = np.degrees(circmean(rad))

    # If std is above max_std → distribution is very broad.
    # For CG data, broad distributions are common and acceptable;
    # they will receive weak restraints (low k) via TDC scoring.
    # Only skip here if truly unreasonably broad AND bimodal.
    if std > max_std:
        # Additional check: is it genuinely bimodal?
        # Compute fraction within ±60° of mean
        mean_wrapped = ((mean + 180) % 360) - 180
        diffs = np.abs(((values_deg - mean_wrapped + 180) % 360) - 180)
        fraction_in_basin = np.sum(diffs < 60.0) / n

        if fraction_in_basin < dominant_fraction:
            return False, (f"bimodal distribution "
                           f"(std={std:.1f}°>{max_std:.0f}°, "
                           f"basin={fraction_in_basin:.2f}<{dominant_fraction:.2f})")
        # Broad but unimodal — accept (weak restraint via TDC)
        return True, ''

    # std <= max_std: check basin fraction as additional verification
    mean_wrapped = ((mean + 180) % 360) - 180
    diffs = np.abs(((values_deg - mean_wrapped + 180) % 360) - 180)
    fraction_in_basin = np.sum(diffs < 60.0) / n

    if fraction_in_basin < dominant_fraction:
        return False, (f"bimodal distribution "
                       f"(std={std:.1f}°, "
                       f"basin={fraction_in_basin:.2f}<{dominant_fraction:.2f})")

    return True, ''


# ══════════════════════════════════════════════════════════════════════════════
# 3.  TDC SCORING — per dihedral
# ══════════════════════════════════════════════════════════════════════════════

def _circ_stats(values_deg, min_samples=10):
    """Return (mean_deg, std_deg) for a circular distribution, [-180,180]."""
    if len(values_deg) < min_samples:
        return np.nan, np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rad  = np.radians(values_deg)
        mean = np.degrees(circmean(rad))
        std  = np.degrees(circstd(rad))
    mean = ((mean + 180) % 360) - 180
    return mean, std


def compute_TDC(distributions, min_samples=10):
    """
    Compute per-dihedral Trajectory Dihedral Confidence (TDC) scores.

    TDC is computed independently for φ and ψ:

        TDC_phi[i] = 100 × (1 - σ_phi / σ_max)
        TDC_psi[i] = 100 × (1 - σ_psi / σ_max)

    This allows stable φ to be restrained independently of flexible ψ,
    avoiding over-restraining of mixed-flexibility residues.

    The residue-level TDC reported in the table is the mean of available
    per-dihedral TDC values (for display purposes only).

    Returns
    -------
    global_TDC    : float   mean TDC across all residues
    global_weight : float   mean weight (terminal residues excluded)
    per_res_TDC   : dict    {resid: float}  residue-level TDC (display)
    per_res_weight: dict    {resid: float}  residue-level weight (display)
    per_res_stats : dict    {resid: dict}   per-dihedral means and stds
    """
    per_res_TDC    = {}
    per_res_weight = {}
    per_res_stats  = {}

    for resid, data in distributions.items():
        phi_arr = np.array(data['phi'])
        psi_arr = np.array(data['psi'])

        phi_mean, phi_std = _circ_stats(phi_arr, min_samples)
        psi_mean, psi_std = _circ_stats(psi_arr, min_samples)

        is_terminal = np.isnan(phi_std) or np.isnan(psi_std)

        # Per-dihedral TDC and k
        tdc_phi = float(np.clip(
            100.0 * (1.0 - phi_std / MAX_CIRCULAR_STD), 0, 100
        )) if not np.isnan(phi_std) else np.nan

        tdc_psi = float(np.clip(
            100.0 * (1.0 - psi_std / MAX_CIRCULAR_STD), 0, 100
        )) if not np.isnan(psi_std) else np.nan

        # Residue-level TDC — mean of available per-dihedral TDCs
        valid_tdcs = [t for t in [tdc_phi, tdc_psi] if not np.isnan(t)]
        if not valid_tdcs:
            per_res_TDC[resid]    = 0.0
            per_res_weight[resid] = 0.0
            per_res_stats[resid]  = {
                'phi_mean': np.nan, 'phi_std': np.nan,
                'psi_mean': np.nan, 'psi_std': np.nan,
                'tdc_phi' : np.nan, 'tdc_psi': np.nan,
            }
            continue

        tdc_residue = float(np.mean(valid_tdcs))

        # Sigmoid weight (for display — actual restraint uses per-dihedral k)
        valid_stds  = [s for s in [phi_std, psi_std] if not np.isnan(s)]
        uncertainty = float(np.mean(valid_stds))
        weight = 0.0 if is_terminal else \
                 1.0 / (1.0 + np.exp((uncertainty - SIG_INFLECTION) / SIG_SCALE))

        per_res_TDC[resid]    = tdc_residue
        per_res_weight[resid] = float(weight)
        per_res_stats[resid]  = {
            'phi_mean': phi_mean, 'phi_std': phi_std,
            'psi_mean': psi_mean, 'psi_std': psi_std,
            'tdc_phi' : tdc_phi,  'tdc_psi': tdc_psi,
        }

    global_TDC    = float(np.nanmean(list(per_res_TDC.values())))
    global_weight = float(np.nanmean(
        [w for w in per_res_weight.values() if w > 0]
    )) if any(w > 0 for w in per_res_weight.values()) else float('nan')

    return global_TDC, global_weight, per_res_TDC, per_res_weight, per_res_stats


# ══════════════════════════════════════════════════════════════════════════════
# 4.  PRINT TDC TABLE
# ══════════════════════════════════════════════════════════════════════════════

def print_TDC_table(distributions, per_res_TDC, per_res_weight, per_res_stats,
                    global_TDC, global_weight):
    """Print per-residue TDC table to stdout."""
    header = (f"\n{'Res':>6}  {'Name':>5}   {'TDC':>6}  {'Weight':>7}"
              f"   {'φ mean':>8}  {'φ std':>7}  {'ψ mean':>8}  {'ψ std':>7}"
              f"  {'Confidence':>10}")
    sep = '-' * 84
    print(header)
    print(sep)

    for resid in sorted(distributions.keys()):
        st  = per_res_stats.get(resid, {})
        tdc = per_res_TDC.get(resid, 0.0)
        w   = per_res_weight.get(resid, 0.0)

        conf = ('HIGH'   if tdc > TDC_HIGH else
                'MEDIUM' if tdc > TDC_LOW  else
                'LOW')

        def _fmt(v, w=8, fallback='NaN'):
            return f"{v:{w}.1f}" if not np.isnan(v) else f"{fallback:>{w}}"

        w_str = f"{w:>7.2f}" if not np.isnan(w) else f"{'nan':>7}"

        print(f"{resid:>6}  {distributions[resid]['resname']:>5}  "
              f"{tdc:>7.1f}  {w_str}"
              f"  {_fmt(st.get('phi_mean', np.nan))}"
              f"  {_fmt(st.get('phi_std',  np.nan), 7)}"
              f"  {_fmt(st.get('psi_mean', np.nan))}"
              f"  {_fmt(st.get('psi_std',  np.nan), 7)}"
              f"  {conf:>10}")

    print(sep)
    w_str = f"{global_weight:.2f}" if not np.isnan(global_weight) else "nan"
    print(f"Global TDC: {global_TDC:.1f}  Global weight: {w_str}\n")


# ══════════════════════════════════════════════════════════════════════════════
# 5.  RESIDUE MAPPING (robust — by resname + chain + sequential position)
# ══════════════════════════════════════════════════════════════════════════════

def build_residue_map(aa_pdb, distributions):
    """
    Build a robust mapping from AA residue sequence position to CG resid.

    Maps by (resname, sequential_position_within_chain) rather than
    raw resid numbers, which may differ between CG and AA structures
    due to non-standard numbering or multi-chain systems.

    Returns
    -------
    aa_to_cg : list of (aa_residue, cg_resid) pairs, in sequence order
    n_mismatch : int  number of resname mismatches (warning if > 0)
    """
    u        = mda.Universe(aa_pdb)
    protein  = u.select_atoms("protein")
    aa_residues = list(protein.residues)

    cg_resids   = sorted(distributions.keys())
    n_aa        = len(aa_residues)
    n_cg        = len(cg_resids)
    n_map       = min(n_aa, n_cg)

    aa_to_cg   = []
    n_mismatch = 0

    for i in range(n_map):
        aa_res    = aa_residues[i]
        cg_resid  = cg_resids[i]
        cg_resname = distributions[cg_resid]['resname']

        if aa_res.resname != cg_resname:
            n_mismatch += 1

        aa_to_cg.append((aa_res, cg_resid))

    if n_mismatch > 0:
        print(f"  WARNING: {n_mismatch} resname mismatches between AA and CG "
              f"— check residue order!")
    if n_aa != n_cg:
        print(f"  WARNING: AA has {n_aa} residues, CG has {n_cg} — "
              f"using first {n_map}")

    return aa_to_cg, n_mismatch


# ══════════════════════════════════════════════════════════════════════════════
# 6.  DIHEDRAL RESTRAINTS ITP — per-dihedral TDC + multimodality filter
# ══════════════════════════════════════════════════════════════════════════════

def write_dihedral_restraints(distributions, per_res_TDC, per_res_stats=None,
                               aa_pdb=None, output_itp=None,
                               k_max=K_MAX, sigma=SIG_INFLECTION, tdc_low=TDC_LOW,
                               fixed_dphi=None, dphi_fixed=None,
                               min_samples=MIN_SAMPLES_UNIMODAL):
    """
    Write GROMACS [ dihedral_restraints ] ITP using AA atom indices.

    Per-dihedral TDC:
        k_phi and k_psi are computed independently from σ_phi and σ_psi.
        A stable φ can be strongly restrained even if ψ is flexible.

    Multimodality filter:
        Bimodal distributions are skipped — their circular mean is
        meaningless and would force the structure into an incorrect
        conformation.

    Atom indices are read from the de novo CG2AT AA output (aa_pdb).
    Residues are matched by sequential position for robustness.

    Parameters
    ----------
    distributions : dict   phi/psi distributions from extract_phi_psi()
    per_res_TDC   : dict   TDC scores from compute_TDC()
    per_res_stats : dict   per-residue statistics from compute_TDC()
    aa_pdb        : str    path to AA structure (CG2AT de novo output)
    output_itp    : str    output .itp file path
    k_max         : float  maximum force constant (kJ/mol/rad²)
                           Default 500 — conservative for simultaneous restraints
    sigma         : float  decay scale for force constant (degrees)
    tdc_low       : float  minimum residue TDC to attempt restraints
    fixed_dphi    : float  if set, use fixed dphi instead of adaptive
    dphi_fixed    : float  alias for fixed_dphi (backward compatibility)
    min_samples   : int    minimum frames for multimodality check
    """
    # Accept both parameter name variants for backward compatibility
    _fixed_dphi = fixed_dphi if fixed_dphi is not None else dphi_fixed

    # ── Build robust residue map ───────────────────────────────────────────────
    aa_to_cg, n_mismatch = build_residue_map(aa_pdb, distributions)
    n_res = len(aa_to_cg)

    u_aa     = mda.Universe(aa_pdb)
    protein  = u_aa.select_atoms("protein")
    aa_residues_all = list(protein.residues)

    lines = [
        '; CG-backbone-derived φ/ψ dihedral restraints\n',
        '; Generated by traj_dih_guide.py (four-bead-backbone representation)\n',
        ';\n',
        '; NOTE: φ/ψ angles derived from new Martini 3 N/CA/C/O backbone\n',
        ';       interaction sites, not exact atomistic coordinates.\n',
        ';\n',
        '; Per-dihedral TDC and force constants:\n',
        f';   k_phi = k_max × exp(-σ_phi / {sigma:.0f}°)   [per-dihedral]\n',
        f';   k_psi = k_max × exp(-σ_psi / {sigma:.0f}°)   [per-dihedral]\n',
        f';   k_max = {k_max:.0f} kJ/mol/rad²\n',
        ';\n',
        '; Multimodality filter: bimodal distributions are skipped.\n',
        f';   min_samples={min_samples}, '
        f'max_std={MAX_STD_UNIMODAL:.1f}°, '
        f'basin_fraction={MIN_BASIN_FRACTION:.1f}\n',
        ';\n',
        f';   dphi = {"fixed " + str(_fixed_dphi) + "°" if _fixed_dphi else "adaptive clip(σ×" + str(DPHI_SCALE) + ", " + str(DPHI_MIN) + "°, " + str(DPHI_MAX) + "°)"}\n',
        f';   TDC < {tdc_low:.0f} (LOW confidence) → residue skipped\n',
        ';\n',
        '[ dihedral_restraints ]\n',
        '; ai    aj    ak    al   type    phi(deg)  dphi(deg)  kfac(kJ/mol/rad2)\n',
    ]

    n_written  = 0
    n_bimodal  = 0
    n_low_tdc  = 0
    n_terminal = 0

    for i, (aa_res, cg_resid) in enumerate(aa_to_cg):
        resname = aa_res.resname
        tdc     = per_res_TDC.get(cg_resid, 0.0)

        # Skip low-confidence residues (per-residue gate)
        if tdc < tdc_low:
            n_low_tdc += 1
            continue

        phi_arr = np.array(distributions[cg_resid]['phi'])
        psi_arr = np.array(distributions[cg_resid]['psi'])

        st = per_res_stats.get(cg_resid, {}) if per_res_stats else {}

        # ── PHI: C(i-1) - N(i) - CA(i) - C(i) ──────────────────────────
        if i > 0 and len(phi_arr) >= 10:
            prev_aa = aa_residues_all[i - 1]
            try:
                ai = prev_aa.atoms.select_atoms("name C")[0].index  + 1
                aj = aa_res.atoms.select_atoms("name N")[0].index   + 1
                ak = aa_res.atoms.select_atoms("name CA")[0].index  + 1
                al = aa_res.atoms.select_atoms("name C")[0].index   + 1

                phi_mean = st.get('phi_mean') if st else None
                phi_std  = st.get('phi_std')  if st else None

                if phi_mean is None or np.isnan(phi_mean):
                    phi_mean, phi_std = _circ_stats(phi_arr)

                if np.isnan(phi_mean):
                    n_terminal += 1
                else:
                    # Per-dihedral TDC filter
                    tdc_phi = float(np.clip(
                        100.0 * (1.0 - phi_std / MAX_CIRCULAR_STD), 0, 100))

                    if tdc_phi < tdc_low:
                        pass  # skip this dihedral
                    else:
                        # Multimodality filter
                        uni_ok, uni_reason = _is_unimodal(phi_arr, min_samples)
                        if not uni_ok:
                            lines.append(
                                f'; SKIP PHI {resname}{aa_res.resid}'
                                f' — {uni_reason}\n')
                            n_bimodal += 1
                        else:
                            k_phi   = k_max * np.exp(-phi_std / sigma)
                            dphi_v  = (_fixed_dphi if _fixed_dphi is not None
                                       else float(np.clip(
                                           phi_std * DPHI_SCALE *
                                           (1.5 if resname == 'GLY' else 1.0),
                                           DPHI_MIN, DPHI_MAX)))
                            lines.append(
                                f'; PHI {resname}{aa_res.resid}  '
                                f'mean={phi_mean:.1f}°  std={phi_std:.1f}°  '
                                f'dphi={dphi_v:.1f}°  k={k_phi:.1f}  '
                                f'TDC_phi={tdc_phi:.1f}\n')
                            lines.append(
                                f'  {ai:<6} {aj:<6} {ak:<6} {al:<6}  1  '
                                f'{phi_mean:>8.2f}  {dphi_v:>5.2f}  {k_phi:>8.2f}\n')
                            n_written += 1
            except IndexError:
                pass

        # ── PSI: N(i) - CA(i) - C(i) - N(i+1) ──────────────────────────
        if i < n_res - 1 and len(psi_arr) >= 10:
            next_aa = aa_residues_all[i + 1]
            try:
                ai = aa_res.atoms.select_atoms("name N")[0].index   + 1
                aj = aa_res.atoms.select_atoms("name CA")[0].index  + 1
                ak = aa_res.atoms.select_atoms("name C")[0].index   + 1
                al = next_aa.atoms.select_atoms("name N")[0].index  + 1

                psi_mean = st.get('psi_mean') if st else None
                psi_std  = st.get('psi_std')  if st else None

                if psi_mean is None or np.isnan(psi_mean):
                    psi_mean, psi_std = _circ_stats(psi_arr)

                if np.isnan(psi_mean):
                    n_terminal += 1
                else:
                    tdc_psi = float(np.clip(
                        100.0 * (1.0 - psi_std / MAX_CIRCULAR_STD), 0, 100))

                    if tdc_psi < tdc_low:
                        pass
                    else:
                        uni_ok, uni_reason = _is_unimodal(psi_arr, min_samples)
                        if not uni_ok:
                            lines.append(
                                f'; SKIP PSI {resname}{aa_res.resid}'
                                f' — {uni_reason}\n')
                            n_bimodal += 1
                        else:
                            k_psi   = k_max * np.exp(-psi_std / sigma)
                            dphi_v  = (_fixed_dphi if _fixed_dphi is not None
                                       else float(np.clip(
                                           psi_std * DPHI_SCALE *
                                           (1.5 if resname == 'GLY' else 1.0),
                                           DPHI_MIN, DPHI_MAX)))
                            lines.append(
                                f'; PSI {resname}{aa_res.resid}  '
                                f'mean={psi_mean:.1f}°  std={psi_std:.1f}°  '
                                f'dphi={dphi_v:.1f}°  k={k_psi:.1f}  '
                                f'TDC_psi={tdc_psi:.1f}\n')
                            lines.append(
                                f'  {ai:<6} {aj:<6} {ak:<6} {al:<6}  1  '
                                f'{psi_mean:>8.2f}  {dphi_v:>5.2f}  {k_psi:>8.2f}\n')
                            n_written += 1
            except IndexError:
                pass

    with open(output_itp, 'w') as f:
        f.writelines(lines)

    # ── Atom index check ──────────────────────────────────────────────────────
    data_lines = [l for l in lines
                  if l.strip() and l.strip()[0].isdigit()]
    max_idx    = (max(int(l.split()[3]) for l in data_lines)
                  if data_lines else 0)
    n_atoms    = len(u_aa.atoms)

    print(f"Dihedral restraints written → {output_itp}")
    print(f"  Restraints written : {n_written}")
    print(f"  Skipped (bimodal)  : {n_bimodal}  "
          f"(multimodality filter — circular mean would be meaningless)")
    print(f"  Skipped (low TDC)  : {n_low_tdc}  "
          f"(TDC < {tdc_low:.0f})")
    print(f"  Skipped (terminal) : {n_terminal}")
    print(f"  Atom indices from  : {aa_pdb}")
    print(f"  Max atom index     : {max_idx}  (AA has {n_atoms} atoms)")
    if max_idx <= n_atoms:
        print(f"  Atom indices OK ✅\n")
    else:
        print(f"  WARNING: max index {max_idx} > {n_atoms} — check mapping!\n")

    return n_written


# ══════════════════════════════════════════════════════════════════════════════
# 7.  RESTRAINT SATISFACTION REPORT
# ══════════════════════════════════════════════════════════════════════════════

def report_restraint_satisfaction(itp_file, final_pdb, output_txt=None):
    """
    Check how well the final structure satisfies the dihedral restraints.

    For each restrained dihedral:
        - reads target φ/ψ and dphi from the ITP file
        - measures actual φ/ψ in the final PDB
        - reports |deviation| and satisfied/violated status

    Parameters
    ----------
    itp_file   : str  path to dihedral_restraints.itp
    final_pdb  : str  path to final_cg2at_restrained.pdb
    output_txt : str  optional output report file

    Returns
    -------
    dict with satisfaction statistics
    """
    import re

    # ── Parse ITP file ─────────────────────────────────────────────────────────
    restraints = []  # list of (ai, aj, ak, al, target_phi, dphi, kfac, comment)
    last_comment = ''

    with open(itp_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(';'):
                last_comment = line
            elif line.strip() and not line.strip().startswith('['):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        ai, aj, ak, al = int(parts[0]), int(parts[1]), \
                                         int(parts[2]), int(parts[3])
                        target = float(parts[5])
                        dphi   = float(parts[6])
                        kfac   = float(parts[7])
                        restraints.append({
                            'ai': ai, 'aj': aj, 'ak': ak, 'al': al,
                            'target': target, 'dphi': dphi, 'kfac': kfac,
                            'label': last_comment.lstrip('; ').strip(),
                        })
                    except (ValueError, IndexError):
                        pass

    if not restraints:
        print("  No restraints found in ITP file.")
        return {}

    # ── Load final structure ───────────────────────────────────────────────────
    u    = mda.Universe(final_pdb)
    pos  = u.atoms.positions  # shape (n_atoms, 3)

    def get_dihedral(ai, aj, ak, al):
        """Compute dihedral angle from 1-based atom indices."""
        try:
            p = [pos[idx - 1] for idx in [ai, aj, ak, al]]
            angle = calc_dihedrals(
                p[0][np.newaxis], p[1][np.newaxis],
                p[2][np.newaxis], p[3][np.newaxis]
            )[0]
            return np.degrees(angle)
        except (IndexError, Exception):
            return np.nan

    # ── Check satisfaction ─────────────────────────────────────────────────────
    satisfied = 0
    violated  = []
    results   = []

    for r in restraints:
        actual = get_dihedral(r['ai'], r['aj'], r['ak'], r['al'])
        if np.isnan(actual):
            continue

        # Circular deviation
        diff = abs(((actual - r['target'] + 180) % 360) - 180)
        ok   = diff <= r['dphi']

        results.append({
            'label'    : r['label'],
            'target'   : r['target'],
            'actual'   : actual,
            'deviation': diff,
            'dphi'     : r['dphi'],
            'kfac'     : r['kfac'],
            'satisfied': ok,
        })

        if ok:
            satisfied += 1
        else:
            violated.append(results[-1])

    total    = len(results)
    pct_ok   = 100.0 * satisfied / total if total > 0 else 0.0
    n_viol   = len(violated)

    # Sort violations by deviation (largest first)
    violated.sort(key=lambda x: -x['deviation'])

    # ── Print report ───────────────────────────────────────────────────────────
    sep = '─' * 65
    print(f"\n  {sep}")
    print(f"  Restraint satisfaction report")
    print(f"  {sep}")
    print(f"  Total restraints checked : {total}")
    print(f"  Satisfied (|Δ| ≤ dphi)  : {satisfied}  ({pct_ok:.1f}%)")
    print(f"  Violated                 : {n_viol}  ({100-pct_ok:.1f}%)")

    if violated:
        print(f"\n  Largest violations (top {min(10, n_viol)}):")
        print(f"  {'Label':<30} {'Target':>8} {'Actual':>8} "
              f"{'|Δ|':>6} {'dphi':>6}  Status")
        print(f"  {'─'*30} {'─'*8} {'─'*8} {'─'*6} {'─'*6}  {'─'*8}")
        for r in violated[:10]:
            print(f"  {r['label']:<30} {r['target']:>8.1f}° "
                  f"{r['actual']:>8.1f}° {r['deviation']:>6.1f}° "
                  f"{r['dphi']:>6.1f}°  ❌")

    print(f"  {sep}\n")

    stats = {
        'total'    : total,
        'satisfied': satisfied,
        'violated' : n_viol,
        'pct_ok'   : pct_ok,
        'violations': violated,
    }

    # ── Write report file ──────────────────────────────────────────────────────
    if output_txt:
        with open(output_txt, 'w') as f:
            f.write(f"Restraint satisfaction report\n")
            f.write(f"ITP:   {itp_file}\n")
            f.write(f"PDB:   {final_pdb}\n\n")
            f.write(f"Total    : {total}\n")
            f.write(f"Satisfied: {satisfied} ({pct_ok:.1f}%)\n")
            f.write(f"Violated : {n_viol}\n\n")
            if violated:
                f.write(f"{'Label':<30} {'Target':>8} {'Actual':>8} "
                        f"{'|Δ|':>6} {'dphi':>6}\n")
                for r in violated:
                    f.write(f"{r['label']:<30} {r['target']:>8.1f} "
                            f"{r['actual']:>8.1f} {r['deviation']:>6.1f} "
                            f"{r['dphi']:>6.1f}\n")
        print(f"  Satisfaction report → {output_txt}")

    return stats


# ══════════════════════════════════════════════════════════════════════════════
# 8.  BACKBONE PDB GENERATION
# ══════════════════════════════════════════════════════════════════════════════

def _write_pdb_atom(serial, name, resname, chain, resseq,
                    x, y, z, occupancy, bfactor, element):
    """Format a single ATOM record in strict PDB fixed-width columns."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:<3s} {chain:1s}"
        f"{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{bfactor:6.2f}          "
        f"{element:>2s}\n"
    )


def generate_backbone_pdb(cg_gro, distributions, per_res_TDC,
                           per_res_weight, output_pdb, chain='A'):
    """
    Write a backbone-only PDB using CG bead positions as coordinates.

    B-factor  = TDC score (0-100, like pLDDT)
    Occupancy = dihedral weight (0-1)

    NOTE: Coordinates are CG interaction site positions, not exact
    atomistic coordinates. This file is for visualisation and
    trajectory-guided alignment only.
    """
    u        = mda.Universe(cg_gro)
    protein  = u.select_atoms("protein")
    residues = protein.residues
    cg_resids = sorted(distributions.keys())

    backbone_atoms = ['N', 'CA', 'C', 'O']
    serial = 1
    lines  = [
        'REMARK  CG-backbone-derived reference for CG2AT\n',
        'REMARK  Coordinates are CG interaction site positions (not exact AA)\n',
        'REMARK  B-factor  = TDC score (0-100, analogous to pLDDT)\n',
        'REMARK  Occupancy = dihedral weight (0-1)\n',
        'REMARK  Generated by traj_dih_guide.py\n',
    ]

    for i, res in enumerate(residues):
        if i >= len(cg_resids):
            break
        cg_resid = cg_resids[i]
        tdc      = per_res_TDC.get(cg_resid, 0.0)
        weight   = per_res_weight.get(cg_resid, 0.0)

        for aname in backbone_atoms:
            sel = res.atoms.select_atoms(f"name {aname}")
            if len(sel) == 0:
                continue
            pos     = sel[0].position
            element = aname[0]
            lines.append(
                _write_pdb_atom(
                    serial, f' {aname}', res.resname, chain, res.resid,
                    pos[0], pos[1], pos[2],
                    weight, tdc, element
                )
            )
            serial += 1

    lines.append('END\n')
    with open(output_pdb, 'w') as f:
        f.writelines(lines)

    print(f"Backbone PDB written → {output_pdb}  ({serial-1} atoms)\n")


# ══════════════════════════════════════════════════════════════════════════════
# 9.  PLOTS
# ══════════════════════════════════════════════════════════════════════════════

def plot_TDC_profile(per_res_TDC, distributions, prefix='protein'):
    """Bar plot of TDC score per residue — analogous to pLDDT profile."""
    resids = sorted(per_res_TDC.keys())
    tdc    = [per_res_TDC[r] for r in resids]
    colors = ['#2196F3' if t > TDC_HIGH else
              '#4CAF50' if t > TDC_LOW  else '#FF9800'
              for t in tdc]

    fig, ax = plt.subplots(figsize=(max(12, len(resids) * 0.15), 4))
    ax.bar(range(len(resids)), tdc, color=colors, alpha=0.85, width=0.9)
    ax.axhline(TDC_HIGH, color='#2196F3', ls='--', lw=1,
               label=f'High (>{TDC_HIGH:.0f})')
    ax.axhline(TDC_LOW,  color='#FF9800', ls='--', lw=1,
               label=f'Low (<{TDC_LOW:.0f})')
    ax.set_ylim(0, 105)
    ax.set_xlim(-0.5, len(resids) - 0.5)
    ax.set_xlabel('Residue index', fontsize=11)
    ax.set_ylabel('TDC Score', fontsize=11)
    ax.set_title(f'{prefix} — Trajectory Dihedral Confidence (TDC)', fontsize=12)
    ax.legend(fontsize=9)
    for i, (resid, t) in enumerate(zip(resids, tdc)):
        if t < TDC_LOW:
            rn = distributions[resid]['resname']
            ax.text(i, t + 2, f'{rn}{resid}',
                    ha='center', va='bottom', fontsize=6, rotation=90)
    plt.tight_layout()
    out = f'{prefix}_TDC_profile.png'
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved TDC profile  → {out}")


def plot_weight_profile(per_res_weight, distributions, prefix='protein'):
    """Bar plot of restraint weight per residue."""
    resids  = sorted(per_res_weight.keys())
    weights = [per_res_weight[r] for r in resids]
    colors  = ['#1565C0' if w > 0.7 else
               '#43A047' if w > 0.3 else '#E65100'
               for w in weights]

    fig, ax = plt.subplots(figsize=(max(12, len(resids) * 0.15), 4))
    ax.bar(range(len(resids)), weights, color=colors, alpha=0.85, width=0.9)
    ax.axhline(0.7, color='#1565C0', ls='--', lw=1, label='High (>0.7)')
    ax.axhline(0.3, color='#E65100', ls='--', lw=1, label='Low (<0.3)')
    ax.set_ylim(0, 1.05)
    ax.set_xlabel('Residue index', fontsize=11)
    ax.set_ylabel('Restraint weight', fontsize=11)
    ax.set_title(f'{prefix} — Per-residue restraint weight', fontsize=12)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = f'{prefix}_weight_profile.png'
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved weight profile → {out}")


def plot_phi_psi_distributions(distributions, per_res_TDC,
                                prefix='protein', n_per_row=6):
    """Per-residue φ/ψ scatter plots coloured by TDC."""
    resids = sorted(distributions.keys())
    n_res  = len(resids)
    n_cols = n_per_row
    n_rows = max(1, (n_res + n_cols - 1) // n_cols)

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 3, n_rows * 3))
    axes = np.array(axes).flatten()

    for idx, resid in enumerate(resids):
        ax      = axes[idx]
        data    = distributions[resid]
        phi_arr = data['phi']
        psi_arr = data['psi']
        tdc     = per_res_TDC.get(resid, 0.0)
        n       = min(len(phi_arr), len(psi_arr))

        if n == 0:
            ax.set_visible(False)
            continue

        color = ('#2196F3' if tdc > TDC_HIGH else
                 '#4CAF50' if tdc > TDC_LOW  else '#FF9800')
        ax.scatter(phi_arr[:n], psi_arr[:n], alpha=0.2, s=3, color=color)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.axhline(0, color='gray', lw=0.4)
        ax.axvline(0, color='gray', lw=0.4)
        ax.set_title(f"{data['resname']}{resid}  TDC={tdc:.0f}", fontsize=7)
        ax.tick_params(labelsize=6)

    for idx in range(n_res, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle(f'{prefix} — CG-backbone-derived φ/ψ distributions',
                 fontsize=11)
    plt.tight_layout()
    out = f'{prefix}_phi_psi_distributions.png'
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Saved φ/ψ distributions → {out}")


def plot_overall_ramachandran(distributions, per_res_TDC, prefix='protein'):
    """Overall Ramachandran — all residues, coloured by TDC."""
    all_phi, all_psi, all_tdc = [], [], []
    for resid, data in distributions.items():
        phi_arr = data['phi']
        psi_arr = data['psi']
        tdc     = per_res_TDC.get(resid, 0.0)
        n = min(len(phi_arr), len(psi_arr))
        if n == 0:
            continue
        all_phi.extend(phi_arr[:n])
        all_psi.extend(psi_arr[:n])
        all_tdc.extend([tdc] * n)

    all_phi = np.array(all_phi)
    all_psi = np.array(all_psi)
    all_tdc = np.array(all_tdc)

    colors = np.where(all_tdc > TDC_HIGH, '#2196F3',
             np.where(all_tdc > TDC_LOW,  '#4CAF50', '#FF9800'))

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(all_phi, all_psi, c=colors, alpha=0.15, s=2)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.axhline(0, color='k', lw=0.5, ls='--')
    ax.axvline(0, color='k', lw=0.5, ls='--')
    ax.set_xlabel('Φ (degrees)', fontsize=12)
    ax.set_ylabel('Ψ (degrees)', fontsize=12)
    ax.set_title(f'{prefix} — CG-backbone-derived Ramachandran\n'
                 'Blue=HIGH TDC  Green=MEDIUM  Orange=LOW', fontsize=11)
    plt.tight_layout()
    out = f'{prefix}_overall_ramachandran.png'
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved overall Ramachandran → {out}")


# ══════════════════════════════════════════════════════════════════════════════
# 10.  INSTRUCTIONS PRINTER
# ══════════════════════════════════════════════════════════════════════════════

def print_restraints_instructions(itp_file):
    """Print step-by-step GROMACS instructions after generating ITP."""
    sep = '=' * 70
    print(f"\n{sep}")
    print("  HOW TO APPLY DIHEDRAL RESTRAINTS")
    print(sep)
    print(f"""
The dihedral restraints ITP has been written to:
  {itp_file}

STEP 1 — Add restraints to the molecule topology
  Open PROTEIN_0.itp and add this line at the very end:

    #include "{os.path.basename(itp_file)}"

  The [ dihedral_restraints ] section MUST be inside the molecule
  .itp file (not in topol.top) for GROMACS to recognise it.

  Quick command:
    echo '#include "{os.path.basename(itp_file)}"' >> PROTEIN_0.itp

STEP 2 — Run energy minimisation (remove clashes first)

  gmx grompp -f minim.mdp -c final_cg2at_de_novo.pdb \\
             -p topol_final.top -o minim.tpr
  gmx mdrun -v -deffnm minim -ntmpi 1 -ntomp 4

STEP 3 — Run short NVT MD with restraints
  NOTE: Minimisation alone CANNOT rotate backbone dihedrals over
  energy barriers. Thermal energy at 300K is required.

  gmx grompp -f nvt_restrained.mdp -c minim.gro \\
             -p topol_final.top -o nvt_restrained.tpr
  gmx mdrun -v -deffnm nvt_restrained -ntmpi 1 -ntomp 4

STEP 4 — Extract final structure
  gmx editconf -f nvt_restrained.gro -o final_restrained.pdb

Or use cg2at_refine.py to automate steps 1-4:
  python cg2at_refine.py -cg2at CG2AT_* -tpr protein.tpr -xtc traj.xtc
""")
    print(sep + '\n')


# ══════════════════════════════════════════════════════════════════════════════
# 11.  ARGUMENT PARSER
# ══════════════════════════════════════════════════════════════════════════════

def build_parser():
    p = argparse.ArgumentParser(
        description=(
            'traj_dih_guide.py — CG-backbone-derived dihedral restraints '
            'for CG2AT.\n'
            'Requires new Martini 3 four-bead backbone (N, CA, C, O).'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MODES
-----
Standard mode (backbone PDB + plots):
  python traj_dih_guide.py -tpr protein.tpr -xtc traj.xtc
                           -c last_frame.gro -o backbone_ref.pdb
                           -step 10 -prefix myprotein

Restraints mode (GROMACS ITP):
  python traj_dih_guide.py --restraints
                           -tpr protein.tpr -xtc traj.xtc
                           -aa final_cg2at_de_novo.pdb
                           -o dihedral_restraints.itp
                           -step 10

Fixed dphi override (recommended for flexible/soluble proteins):
  python traj_dih_guide.py --restraints ... --dphi 10.0

Conservative kmax (default):
  k_max = 500 kJ/mol/rad²  (avoids over-restraining many simultaneous dihedrals)
  Use --kmax 1000 for stronger restraints.
""")

    p.add_argument('--restraints', action='store_true',
                   help='Restraints mode: generate GROMACS ITP file')
    p.add_argument('-tpr',  required=True,
                   help='CG topology file (.tpr or .gro)')
    p.add_argument('-xtc',  required=True,
                   help='CG trajectory file (.xtc)')
    p.add_argument('-c',    default=None,
                   help='[standard] Last frame .gro for backbone PDB')
    p.add_argument('-aa',   default=None,
                   help='[restraints] AA de novo PDB from CG2AT')
    p.add_argument('-o',    required=True,
                   help='Output file (.pdb standard / .itp restraints)')
    p.add_argument('-prefix', default='protein',
                   help='Prefix for plot filenames (default: protein)')
    p.add_argument('-start', type=int, default=0,
                   help='First frame index (default: 0)')
    p.add_argument('-stop',  type=int, default=None,
                   help='Last frame index (default: all)')
    p.add_argument('-step',  type=int, default=1,
                   help='Frame stride (default: 1)')
    p.add_argument('--kmax',  type=float, default=K_MAX,
                   help=f'Max force constant kJ/mol/rad² (default: {K_MAX})')
    p.add_argument('--sigma', type=float, default=SIG_INFLECTION,
                   help=f'Force constant decay scale degrees (default: {SIG_INFLECTION})')
    p.add_argument('--dphi',  type=float, default=None,
                   help='Fixed dphi tolerance degrees (default: adaptive)')
    p.add_argument('--min-samples', type=int, default=MIN_SAMPLES_UNIMODAL,
                   help=f'Min frames for multimodality check (default: {MIN_SAMPLES_UNIMODAL})')
    p.add_argument('--no-plots', action='store_true',
                   help='Skip generating plot files')
    p.add_argument('--validate-only', action='store_true',
                   help='Only validate backbone — do not extract phi/psi')
    return p


# ══════════════════════════════════════════════════════════════════════════════
# 12.  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    sep = '=' * 70
    print(f"\n{sep}")
    print("  traj_dih_guide.py — CG-backbone-derived dihedral restraints")
    print(f"{sep}\n")

    parser = build_parser()
    args   = parser.parse_args()

    # ── Backbone validation (point 7) ─────────────────────────────────────────
    print("Validating backbone representation ...")
    validate_backbone(args.tpr, args.xtc)
    print()

    if args.validate_only:
        print("Backbone validation complete — exiting (--validate-only).")
        return

    # ── Extract φ/ψ distributions ─────────────────────────────────────────────
    distributions = extract_phi_psi(
        args.tpr, args.xtc,
        start=args.start, stop=args.stop, step=args.step
    )

    # ── Compute per-dihedral TDC ──────────────────────────────────────────────
    global_TDC, global_weight, per_res_TDC, per_res_weight, per_res_stats = \
        compute_TDC(distributions)

    w_str = f"{global_weight:.2f}" if not np.isnan(global_weight) else "nan"
    print(f"Global TDC score: {global_TDC:.1f}/100")
    print(f"Global weight:    {w_str}\n")

    print_TDC_table(distributions, per_res_TDC, per_res_weight,
                    per_res_stats, global_TDC, global_weight)

    # ── Restraints mode ───────────────────────────────────────────────────────
    if args.restraints:
        if args.aa is None:
            parser.error("--restraints mode requires -aa <de_novo.pdb>")

        print("Generating dihedral restraints ITP ...")
        print(f"  Using AA atom indices from: {args.aa}\n")

        write_dihedral_restraints(
            distributions, per_res_TDC, per_res_stats,
            aa_pdb      = args.aa,
            output_itp  = args.o,
            k_max       = args.kmax,
            sigma       = args.sigma,
            tdc_low     = TDC_LOW,
            fixed_dphi  = args.dphi,
            min_samples = args.min_samples,
        )
        print_restraints_instructions(args.o)
        return

    # ── Standard mode: backbone PDB + plots ──────────────────────────────────
    if args.c is None:
        parser.error("Standard mode requires -c <last_frame.gro>")

    print("Generating backbone reference PDB ...")
    generate_backbone_pdb(
        args.c, distributions, per_res_TDC, per_res_weight, args.o
    )

    if not args.no_plots:
        print("Generating plots ...")
        plot_TDC_profile(per_res_TDC, distributions, prefix=args.prefix)
        plot_weight_profile(per_res_weight, distributions, prefix=args.prefix)
        plot_phi_psi_distributions(distributions, per_res_TDC,
                                   prefix=args.prefix)
        plot_overall_ramachandran(distributions, per_res_TDC,
                                  prefix=args.prefix)

    print(f"\nDone. Backbone reference: {args.o}")
    print("NOTE: Backbone PDB is for visualisation — "
          "use --restraints for GROMACS ITP.\n")


if __name__ == '__main__':
    main()
