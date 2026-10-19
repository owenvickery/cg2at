"""
cg2at_refine.py
===============
Post-processing refinement for CG2AT.

Applies trajectory-informed dihedral restraints (TDC method) to improve
Ramachandran quality of CG→AA backmapped structures — without biasing
toward any reference structure.

Workflow
--------
    # Step 1 — run CG2AT as normal
    cg2at.py -c last_frame.gro -fg martini_3-2_nextgen_protein_charmm36 \\
             -ff charmm36-cg2at -w tip3p

    # Step 2 — refine with trajectory dihedral restraints
    python cg2at_refine.py \\
        -cg2at CG2AT_2026-08-01_17-07-15 \\
        -tpr   protein_CG.tpr \\
        -xtc   trajectory.xtc \\
        -step  10

    # Output: CG2AT_2026-08-01_17-07-15/FINAL/final_cg2at_restrained.pdb

Scientific basis
----------------
TDC (Trajectory Dihedral Confidence)
    Analogous to pLDDT but measuring dihedral consistency.
    Narrow φ/ψ distributions → high TDC → strong restraints.
    Broad φ/ψ distributions  → low TDC  → weak/no restraints.

    TDC = 100 × (1 - σ_eff / σ_max)
    k   = k_max × exp(-σ_eff / σ_inflection)

Flags
-----
    Required:
      -cg2at    Path to CG2AT output folder (e.g. CG2AT_2026-08-01_17-07-15)
      -tpr      CG topology file (.tpr or .gro)
      -xtc      CG trajectory file (.xtc)

    Trajectory options:
      -step     Frame stride (default: 10)
      -b        Start time in ps (default: 0 = beginning)
      -e        End time in ps   (default: -1 = end)

    Restraint options:
      --dphi    Fixed dphi tolerance in degrees (default: adaptive σ×0.3)
      --kmax    Maximum force constant kJ/mol/rad² (default: 500)
      --sigma   Sigmoid inflection in degrees (default: 30)
      --tdc-low Minimum TDC to apply restraint (default: 40)

    GROMACS options:
      --gmx         GROMACS executable (default: auto-detect from log)
      --nsteps-nvt  NVT steps (default: 50000 = 100 ps at 2 fs)
      --temp        Temperature in K (default: 300)
      --no-nvt      Skip NVT — only run minimisation
      --ntmpi       MPI threads for mdrun (default: 1)
      --ntomp       OpenMP threads for mdrun (default: 4)

    Input options:
      --aa-input    AA structure for restraint template and GROMACS input.
                    'de_novo'  : use final_cg2at_de_novo.pdb  (default)
                    'aligned'  : use final_cg2at_aligned.pdb
                    Use 'aligned' when an experimental or AF3 reference
                    was provided via CG2AT -a flag.

    Output options:
      --keep-tmp    Keep intermediate GROMACS files in FINAL/
      --no-plots    Skip Ramachandran plot generation
      --prefix      Output file prefix (default: final_cg2at_restrained)

Requirements
    MDAnalysis >= 2.0
    numpy, scipy, matplotlib
    GROMACS 2020+ (auto-detected from cg2at output logs)

Notes
-----
    The Ramachandran background uses MDAnalysis built-in reference data
    (Richardson lab Top8000 dataset) stored locally in the MDAnalysis
    package — no internet access required.

    MDAnalysis only provides a general Ramachandran reference (Rama_ref).
    There is NO separate GLY or PRO reference in MDAnalysis.
    GLY and PRO classification uses analytical boundaries from
    Lovell et al. 2003 (the same source MDAnalysis cites).

    GLY and PRO residues are excluded from the general Ramachandran plot
    (as in Swiss-Model) because they have distinct allowed regions.
    Their count is reported in the footnote below each panel.
    A separate figure (ramachandran_gly_pro_restrained.png) shows
    GLY and PRO panels with residue-specific analytical backgrounds.

    If rama500 data files (from Lovell et al. 2003) are placed in a
    'rama_data/' folder next to this script, they will be used for
    GLY/PRO background contours — matching RAMPAGE/MolProbity exactly:

        mkdir -p rama_data
        wget -P rama_data https://github.com/donaldlab/BWM/raw/refs/heads/master/data/rama500-general.data
        wget -P rama_data https://github.com/donaldlab/BWM/raw/refs/heads/master/data/rama500-gly-sym.data
        wget -P rama_data https://github.com/donaldlab/BWM/raw/refs/heads/master/data/rama500-pro.data
        wget -P rama_data https://github.com/donaldlab/BWM/raw/refs/heads/master/data/rama500-prepro.data

    Falls back to analytical boundaries if files are not found.

Author : Hafez Razmazma
Contact: hafez.razmazma@warwick.ac.uk
"""

import os
import sys
import shutil
import subprocess
import warnings
import argparse
import textwrap
import time
import numpy as np
from pathlib import Path
from datetime import datetime

# ── Add script directory to path so traj_dih_guide.py is always found ────────
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

# ── Matplotlib backend — must be before pyplot ────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
except Exception:
    pass
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
from scipy.stats import circmean, circstd

# ── Suppress expected MDAnalysis warnings ─────────────────────────────────────
warnings.filterwarnings('ignore', message='Cannot determine phi and psi angles')
warnings.filterwarnings('ignore', message='Reader has no dt information')
warnings.filterwarnings('ignore', message='Element information is missing')
warnings.filterwarnings('ignore', message='Found no information for attr')
warnings.filterwarnings('ignore', message='Found missing chainIDs')

# ── Import core functions from traj_dih_guide ────────────────────────────────
try:
    try:
        from traj_dih_guide import (        # database/bin on sys.path (clone/conda, or via cg2at)
            extract_phi_psi,
            compute_TDC,
            write_dihedral_restraints,
            print_TDC_table,
        )
    except ImportError:
        from cg2at.bin.traj_dih_guide import (   # pip-installed console script
            extract_phi_psi,
            compute_TDC,
            write_dihedral_restraints,
            print_TDC_table,
        )
    _TDG_AVAILABLE = True
except ImportError as e:
    print(f"  WARNING: traj_dih_guide could not be imported ({e}); "
          f"dihedral restraints will be skipped.")
    _TDG_AVAILABLE = False


# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════
K_MAX            = 500.0
SIG_INFLECTION   = 30.0
MAX_CIRCULAR_STD = 82.7
TDC_LOW          = 40.0
DPHI_MIN         = 5.0
DPHI_MAX         = 20.0
DPHI_SCALE       = 0.3

# MDAnalysis Ramachandran contour levels (from MDAnalysis source)
# Z >= MDA_LEVEL_FAV  → favoured (dark blue)
# Z >= MDA_LEVEL_ALL  → allowed  (light blue)
# Z <  MDA_LEVEL_ALL  → outlier  (white)
MDA_LEVEL_FAV    = 17
MDA_LEVEL_ALL    = 1

# Dot colours — visible on MDAnalysis blue background
COLOR_FAV  = '#222222'   # near-black  — favoured
COLOR_ALL  = '#E07B00'   # amber       — allowed
COLOR_OUT  = '#B71C1C'   # dark red    — outlier
DOT_SIZE   = 35
DOT_ALPHA  = 0.90

# Residues excluded from general Ramachandran (same as Swiss-Model)
SPECIAL_RES = frozenset({'GLY', 'PRO'})

OUT_RESTRAINED_PDB   = 'final_cg2at_restrained.pdb'
OUT_TDC_PDB          = 'final_cg2at_restrained_TDC.pdb'
OUT_RESTRAINTS_ITP   = 'dihedral_restraints.itp'
OUT_QUALITY_DAT      = 'structure_quality_restrained.dat'
OUT_RAMA_PNG         = 'ramachandran_restrained.png'
OUT_RAMA_GLYRO_PNG   = 'ramachandran_gly_pro_restrained.png'
OUT_LOG              = 'cg2at_refine.log'
GROMACS_SUBDIR       = 'gromacs_outputs_refine'
MDP_MINIM_NAME       = 'minim_dih_res.mdp'
MDP_NVT_NAME         = 'nvt_dih_res.mdp'

# ── Shared plot style constants (applied identically to ALL panels) ───────────
PLOT_FIGSIZE_PER_PANEL = (5.5, 6.0)
PLOT_TITLE_FONTSIZE    = 11
PLOT_LABEL_FONTSIZE    = 10
PLOT_TICK_FONTSIZE     = 9
PLOT_TABLE_FONTSIZE    = 8
PLOT_TABLE_VAL_FONTSIZE= 9
PLOT_SUPTITLE_FONTSIZE = 13
PLOT_TICKS             = [-180, -90, 0, 90, 180]
PLOT_TICK_LABELS       = ['-180°', '-90°', '0°', '90°', '180°']


# ══════════════════════════════════════════════════════════════════════════════
# Timing — per-step stopwatch
# ══════════════════════════════════════════════════════════════════════════════

class StepTimer:
    """Record wall-clock time for each named step."""

    def __init__(self):
        self._start  = time.time()
        self._steps  = []
        self._t_step = time.time()

    def tick(self, label):
        elapsed = time.time() - self._t_step
        self._steps.append((label, elapsed))
        self._t_step = time.time()

    def total(self):
        return time.time() - self._start

    def _fmt(self, seconds):
        h = int(seconds // 3600)
        m = int((seconds % 3600) // 60)
        s = int(seconds % 60)
        return f"{h} hours {m:2d} min {s:2d} sec"

    def report(self):
        W   = 100
        sep = '-' * W
        lines = [
            sep, '',
            f"{'Job':<50} {'Time':>30}",
            f"{'---':<50} {'----':>30}", '',
        ]
        for label, elapsed in self._steps:
            lines.append(f"  {label:<48} {self._fmt(elapsed):>30}")
        lines += [
            sep,
            f"  {'Total run time:':<48} {self._fmt(self.total()):>30}",
            sep, '',
        ]
        return '\n'.join(lines)


# ══════════════════════════════════════════════════════════════════════════════
# Logging — tee to terminal + log file
# ══════════════════════════════════════════════════════════════════════════════

class Tee:
    def __init__(self, log_path):
        self._terminal = sys.stdout
        self._log      = open(log_path, 'w')
        self._log_only = False

    def write(self, msg):
        if not self._log_only:
            self._terminal.write(msg)
        self._log.write(msg)

    def flush(self):
        self._terminal.flush()
        self._log.flush()

    def close(self):
        self._log.close()

    def set_log_only(self, flag):
        self._log_only = flag


# ══════════════════════════════════════════════════════════════════════════════
# 1.  AUTO-DETECT GROMACS
# ══════════════════════════════════════════════════════════════════════════════

def find_gromacs_executable(cg2at_folder):
    for root, dirs, files in os.walk(cg2at_folder):
        for fname in ['gromacs_outputs', 'gromacs_output']:
            fpath = os.path.join(root, fname)
            if os.path.isfile(fpath):
                try:
                    with open(fpath, 'r', errors='ignore') as f:
                        for line in f:
                            if line.strip().startswith('Executable:'):
                                gmx = line.split('Executable:')[1].strip()
                                if os.path.isfile(gmx):
                                    return gmx
                except Exception:
                    pass
    for candidate in ['gmx', 'gmx_mpi', '/usr/bin/gmx']:
        if shutil.which(candidate):
            return candidate
    return 'gmx'


# ══════════════════════════════════════════════════════════════════════════════
# 2.  LOCATE CG2AT OUTPUT FILES
# ══════════════════════════════════════════════════════════════════════════════

def find_cg2at_outputs(cg2at_folder):
    cg2at_folder = Path(cg2at_folder).resolve()
    if not cg2at_folder.exists():
        sys.exit(f"ERROR: CG2AT folder not found: {cg2at_folder}")

    final_dir  = cg2at_folder / 'FINAL'
    merged_dir = cg2at_folder / 'MERGED'
    if not final_dir.exists():
        sys.exit(f"ERROR: FINAL/ subfolder not found in {cg2at_folder}")

    outputs = {
        'cg2at_folder': str(cg2at_folder),
        'final_dir'   : str(final_dir),
        'merged_dir'  : str(merged_dir) if merged_dir.exists() else None,
    }

    de_novo = final_dir / 'final_cg2at_de_novo.pdb'
    if not de_novo.exists():
        sys.exit(f"ERROR: final_cg2at_de_novo.pdb not found in {final_dir}")
    outputs['de_novo_pdb'] = str(de_novo)

    protein_itp = final_dir / 'PROTEIN_0.itp'
    if not protein_itp.exists():
        sys.exit(f"ERROR: PROTEIN_0.itp not found in {final_dir}")
    outputs['protein_itp'] = str(protein_itp)

    topol_top = final_dir / 'topol_final.top'
    if not topol_top.exists():
        sys.exit(f"ERROR: topol_final.top not found in {final_dir}")
    outputs['topol_top'] = str(topol_top)

    aligned = final_dir / 'final_cg2at_aligned.pdb'
    outputs['aligned_pdb'] = str(aligned) if aligned.exists() else None

    # Auto-detect CG reference: CG2AT_*/INPUT/CG_INPUT.pdb
    cg_input = cg2at_folder / 'INPUT' / 'CG_INPUT.pdb'
    outputs['cg_input_pdb'] = str(cg_input) if cg_input.exists() else None

    return outputs


# ══════════════════════════════════════════════════════════════════════════════
# 3.  TRAJECTORY TIME RANGE → FRAME INDICES
# ══════════════════════════════════════════════════════════════════════════════

def parse_time_range(tpr_file, xtc_file, b_ps, e_ps, step):
    u     = mda.Universe(tpr_file, xtc_file)
    times = np.array([ts.time for ts in u.trajectory])
    start = 0 if b_ps <= 0 else min(np.searchsorted(times, b_ps), len(times) - 1)
    stop  = None if e_ps < 0 else min(np.searchsorted(times, e_ps, side='right'), len(times))
    n_frames = len(u.trajectory[start:stop:step])
    return start, stop, n_frames


# ══════════════════════════════════════════════════════════════════════════════
# 4.  INLINE PHI/PSI EXTRACTION (fallback)
# ══════════════════════════════════════════════════════════════════════════════

def _extract_phi_psi_inline(tpr_file, xtc_file, start=0, stop=None, step=1):
    u        = mda.Universe(tpr_file, xtc_file)
    protein  = u.select_atoms("protein")
    residues = protein.residues
    n_res    = len(residues)

    distributions = {
        res.resid: {'resname': res.resname, 'phi': [], 'psi': []}
        for res in residues
    }

    n_proc = len(u.trajectory[start:stop:step])
    print(f"Found {n_res} residues | {len(u.trajectory)} total frames")
    print(f"Using every {step} frame(s) → processing {n_proc} frames")

    frame_count = 0
    for ts in u.trajectory[start:stop:step]:
        for i, res in enumerate(residues):
            if i > 0:
                prev = residues[i - 1]
                try:
                    C_prev  = prev.atoms.select_atoms("name C")[0].position
                    N_curr  = res.atoms.select_atoms("name N")[0].position
                    CA_curr = res.atoms.select_atoms("name CA")[0].position
                    C_curr  = res.atoms.select_atoms("name C")[0].position
                    phi = calc_dihedrals(
                        C_prev[np.newaxis], N_curr[np.newaxis],
                        CA_curr[np.newaxis], C_curr[np.newaxis])[0]
                    distributions[res.resid]['phi'].append(np.degrees(phi))
                except Exception:
                    pass
            if i < n_res - 1:
                nxt = residues[i + 1]
                try:
                    N_curr  = res.atoms.select_atoms("name N")[0].position
                    CA_curr = res.atoms.select_atoms("name CA")[0].position
                    C_curr  = res.atoms.select_atoms("name C")[0].position
                    N_next  = nxt.atoms.select_atoms("name N")[0].position
                    psi = calc_dihedrals(
                        N_curr[np.newaxis], CA_curr[np.newaxis],
                        C_curr[np.newaxis], N_next[np.newaxis])[0]
                    distributions[res.resid]['psi'].append(np.degrees(psi))
                except Exception:
                    pass
        frame_count += 1
        if frame_count % 200 == 0:
            print(f"  Processed {frame_count}/{n_proc} frames ...")

    for resid in distributions:
        distributions[resid]['phi'] = np.array(distributions[resid]['phi'])
        distributions[resid]['psi'] = np.array(distributions[resid]['psi'])

    print(f"Done — {frame_count} frames processed.\n")
    return distributions


# ══════════════════════════════════════════════════════════════════════════════
# 5.  INLINE TDC SCORING (fallback)
# ══════════════════════════════════════════════════════════════════════════════

def _compute_TDC_inline(distributions, min_samples=10):
    per_res_TDC    = {}
    per_res_weight = {}

    for resid, data in distributions.items():
        phi_arr = np.array(data['phi'])
        psi_arr = np.array(data['psi'])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            phi_std = np.degrees(circstd(np.radians(phi_arr))) \
                      if len(phi_arr) >= min_samples else np.nan
            psi_std = np.degrees(circstd(np.radians(psi_arr))) \
                      if len(psi_arr) >= min_samples else np.nan

        valid = [s for s in [phi_std, psi_std] if not np.isnan(s)]
        if not valid:
            per_res_TDC[resid]    = 0.0
            per_res_weight[resid] = 0.0
            continue

        uncertainty    = float(np.mean(valid))
        tdc            = 100.0 * max(0.0, 1.0 - uncertainty / MAX_CIRCULAR_STD)
        weight         = 1.0 / (1.0 + np.exp((uncertainty - SIG_INFLECTION) / 10.0))
        is_terminal    = np.isnan(phi_std) or np.isnan(psi_std)
        per_res_TDC[resid]    = tdc
        per_res_weight[resid] = 0.0 if is_terminal else weight

    global_TDC    = float(np.nanmean(list(per_res_TDC.values())))
    global_weight = float(np.nanmean(
        [w for w in per_res_weight.values() if w > 0]
    )) if any(w > 0 for w in per_res_weight.values()) else 0.0

    return global_TDC, global_weight, per_res_TDC, per_res_weight


# ══════════════════════════════════════════════════════════════════════════════
# 6.  INLINE DIHEDRAL RESTRAINT WRITER (fallback)
# ══════════════════════════════════════════════════════════════════════════════

def _write_restraints_inline(distributions, per_res_TDC, aa_pdb, output_itp,
                              k_max=K_MAX, sigma=SIG_INFLECTION,
                              tdc_low=TDC_LOW, dphi_fixed=None, min_samples=10):
    u        = mda.Universe(aa_pdb)
    protein  = u.select_atoms("protein")
    residues = protein.residues
    n_res    = len(residues)

    lines = [
        '; Dihedral restraints from CG trajectory φ/ψ distributions\n',
        '; Generated by cg2at_refine.py (TDC method)\n',
        ';\n',
        '[ dihedral_restraints ]\n',
        '; ai    aj    ak    al   type    phi(deg)  dphi(deg)  kfac(kJ/mol/rad2)\n',
    ]

    n_written = 0
    for i, res in enumerate(residues):
        cg_resid = i + 1
        if cg_resid not in distributions:
            continue
        tdc = per_res_TDC.get(cg_resid, 0.0)
        if tdc < tdc_low:
            continue

        phi_arr = np.array(distributions[cg_resid]['phi'])
        psi_arr = np.array(distributions[cg_resid]['psi'])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if i > 0 and len(phi_arr) >= min_samples:
                prev = residues[i - 1]
                try:
                    ai = prev.atoms.select_atoms("name C")[0].index + 1
                    aj = res.atoms.select_atoms("name N")[0].index  + 1
                    ak = res.atoms.select_atoms("name CA")[0].index + 1
                    al = res.atoms.select_atoms("name C")[0].index  + 1
                    phi_mean = ((np.degrees(circmean(np.radians(phi_arr)))+180)%360)-180
                    phi_std  = np.degrees(circstd(np.radians(phi_arr)))
                    k        = k_max * np.exp(-phi_std / sigma)
                    dphi     = float(dphi_fixed) if dphi_fixed is not None else \
                               float(np.clip(phi_std * DPHI_SCALE, DPHI_MIN, DPHI_MAX))
                    lines.append(f'; PHI {res.resname}{res.resid}  mean={phi_mean:.1f}°\n')
                    lines.append(f'  {ai:<6} {aj:<6} {ak:<6} {al:<6}  1  '
                                 f'{phi_mean:>8.2f}  {dphi:>5.2f}  {k:>8.2f}\n')
                    n_written += 1
                except Exception:
                    pass
            if i < n_res - 1 and len(psi_arr) >= min_samples:
                nxt = residues[i + 1]
                try:
                    ai = res.atoms.select_atoms("name N")[0].index  + 1
                    aj = res.atoms.select_atoms("name CA")[0].index + 1
                    ak = res.atoms.select_atoms("name C")[0].index  + 1
                    al = nxt.atoms.select_atoms("name N")[0].index  + 1
                    psi_mean = ((np.degrees(circmean(np.radians(psi_arr)))+180)%360)-180
                    psi_std  = np.degrees(circstd(np.radians(psi_arr)))
                    k        = k_max * np.exp(-psi_std / sigma)
                    dphi     = float(dphi_fixed) if dphi_fixed is not None else \
                               float(np.clip(psi_std * DPHI_SCALE, DPHI_MIN, DPHI_MAX))
                    lines.append(f'; PSI {res.resname}{res.resid}  mean={psi_mean:.1f}°\n')
                    lines.append(f'  {ai:<6} {aj:<6} {ak:<6} {al:<6}  1  '
                                 f'{psi_mean:>8.2f}  {dphi:>5.2f}  {k:>8.2f}\n')
                    n_written += 1
                except Exception:
                    pass

    with open(output_itp, 'w') as f:
        f.writelines(lines)

    data_lines = [l for l in lines if l.strip() and l[0] not in (';', '[')]
    max_idx    = max(int(l.split()[3]) for l in data_lines) if data_lines else 0
    n_atoms    = len(mda.Universe(aa_pdb).atoms)
    print(f"  Restraints written : {n_written}")
    print(f"  Max atom index     : {max_idx}  (AA structure has {n_atoms} atoms)")
    print(f"  Atom indices OK ✅" if max_idx <= n_atoms else
          f"  WARNING: max atom index {max_idx} > {n_atoms} ❌")
    return n_written


# ══════════════════════════════════════════════════════════════════════════════
# 7.  GROMACS-STYLE BACKUP
# ══════════════════════════════════════════════════════════════════════════════

def gromacs_style_backup(filepath):
    p = Path(filepath)
    n = 1
    while True:
        backup = p.parent / f'#{p.name}.{n}#'
        if not backup.exists():
            shutil.copy2(filepath, backup)
            return str(backup)
        n += 1


# ══════════════════════════════════════════════════════════════════════════════
# 8.  PATCH PROTEIN_0.ITP
# ══════════════════════════════════════════════════════════════════════════════

def patch_protein_itp(protein_itp_path, itp_filename='dihedral_restraints.itp'):
    with open(protein_itp_path, 'r') as f:
        content = f.read()

    include_line = f'#include "{itp_filename}"'
    if include_line in content:
        print(f"  PROTEIN_0.itp already includes {itp_filename} — skipping patch")
        return False

    backup = gromacs_style_backup(protein_itp_path)
    print(f"  Backup → {backup}")

    with open(protein_itp_path, 'a') as f:
        f.write(f'\n; Trajectory dihedral restraints (cg2at_refine.py)\n')
        f.write('#ifdef DIHEDRALPOSRES\n')
        f.write(f'{include_line}\n')
        f.write('#endif\n')

    print(f"  Added to PROTEIN_0.itp: {include_line}")
    return True


# ══════════════════════════════════════════════════════════════════════════════
# 9.  MDP FILE WRITERS
# ══════════════════════════════════════════════════════════════════════════════

def write_mdp_files(final_dir, nsteps_nvt=50000, temp=300):
    minim_mdp = os.path.join(final_dir, MDP_MINIM_NAME)
    nvt_mdp   = os.path.join(final_dir, MDP_NVT_NAME)

    with open(minim_mdp, 'w') as f:
        f.write(textwrap.dedent(f"""\
            ; Steepest descent minimisation — dihedral restraints active via topology
            ; Generated by cg2at_refine.py
            define        = -DDIHEDRALPOSRES
            integrator    = steep
            nsteps        = 100000
            emtol         = 10.0
            emstep        = 0.001
            cutoff-scheme = Verlet
            nstlist       = 10
            rcoulomb      = 1.2
            rvdw          = 1.2
            pbc           = xyz
        """))

    with open(nvt_mdp, 'w') as f:
        f.write(textwrap.dedent(f"""\
            ; Short NVT MD — dihedral restraints active via topology
            ; Thermal energy at {temp} K allows backbone to rotate toward target φ/ψ
            ; Generated by cg2at_refine.py
            define        = -DDIHEDRALPOSRES
            integrator    = md
            nsteps        = {nsteps_nvt}
            dt            = 0.001
            tcoupl        = V-rescale
            tc-grps       = System
            tau-t         = 0.1
            ref-t         = {temp}
            constraints   = h-bonds
            constraint-algorithm = LINCS
            coulombtype   = PME
            pme_order     = 4
            fourierspacing = 0.135
            pcoupl        = no
            cutoff-scheme = Verlet
            nstlist       = 25
            rcoulomb      = 1.2
            rvdw          = 1.2
            pbc           = xyz
            nstxout       = 5000
            nstenergy     = 1000
            nstlog        = 1000
            nstxout-compressed = 1000
        """))

    print(f"  MDP files written: {os.path.basename(minim_mdp)}, {os.path.basename(nvt_mdp)}")
    return minim_mdp, nvt_mdp


# ══════════════════════════════════════════════════════════════════════════════
# 10.  GROMACS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run_gromacs_pipeline(gmx, final_dir, aa_input_pdb, topol_top,
                          minim_mdp, nvt_mdp,
                          skip_nvt=False, ntmpi=1, ntomp=4,
                          out_prefix='final_cg2at_restrained'):
    gmx_dir = os.path.join(final_dir, GROMACS_SUBDIR)
    os.makedirs(gmx_dir, exist_ok=True)

    def run_cmd(cmd, log_name, label):
        log_path = os.path.join(gmx_dir, log_name)
        result   = subprocess.run(
            cmd, cwd=final_dir,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        with open(log_path, 'w') as lf:
            lf.write(' '.join(cmd) + '\n\n')
            lf.write(result.stdout)
        if result.returncode != 0:
            print(f"\n  ERROR in {label}")
            print(f"  See log: {log_path}")
            for line in result.stdout.splitlines()[-20:]:
                print(f"    {line}")
            sys.exit(1)
        print(f"  {label} ✅")

    def rel(path):
        return os.path.relpath(path, final_dir)

    minim_tpr = os.path.join(gmx_dir, 'minim_dih_res.tpr')
    minim_gro = os.path.join(gmx_dir, 'minim_dih_res.gro')
    nvt_tpr   = os.path.join(gmx_dir, 'nvt_dih_res.tpr')
    nvt_gro   = os.path.join(gmx_dir, 'nvt_dih_res.gro')
    out_pdb   = os.path.join(final_dir, f'{out_prefix}.pdb')

    run_cmd([gmx, 'grompp', '-f', rel(minim_mdp), '-c', rel(aa_input_pdb),
             '-p', rel(topol_top), '-o', rel(minim_tpr), '-maxwarn', '5'],
            'grompp_minim.log', 'Minimisation grompp')
    run_cmd([gmx, 'mdrun', '-v', '-s', rel(minim_tpr),
             '-deffnm', rel(minim_gro).replace('.gro', ''),
             '-ntmpi', str(ntmpi), '-ntomp', str(ntomp)],
            'mdrun_minim.log', 'Minimisation mdrun ')

    input_for_editconf = minim_gro

    if not skip_nvt:
        run_cmd([gmx, 'grompp', '-f', rel(nvt_mdp), '-c', rel(minim_gro),
                 '-p', rel(topol_top), '-o', rel(nvt_tpr), '-maxwarn', '5'],
                'grompp_nvt.log', 'NVT grompp        ')
        run_cmd([gmx, 'mdrun', '-v', '-s', rel(nvt_tpr),
                 '-deffnm', rel(nvt_gro).replace('.gro', ''),
                 '-ntmpi', str(ntmpi), '-ntomp', str(ntomp)],
                'mdrun_nvt.log', 'NVT mdrun         ')
        input_for_editconf = nvt_gro

    run_cmd([gmx, 'editconf', '-f', rel(input_for_editconf), '-o', rel(out_pdb)],
            'editconf.log', 'editconf          ')

    for fname in [MDP_MINIM_NAME, MDP_NVT_NAME, 'mdout.mdp']:
        src = os.path.join(final_dir, fname)
        dst = os.path.join(gmx_dir, fname)
        if os.path.exists(src):
            shutil.move(src, dst)

    return out_pdb


# ══════════════════════════════════════════════════════════════════════════════
# 11.  BACKBONE RMSD — Kabsch alignment (matches CG2AT method)
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# 11.  BACKBONE RMSD  (CG2AT methodology, standard ÷N formula)
# ══════════════════════════════════════════════════════════════════════════════

def _detect_cg_backbone_type(cg_file):
    """
    Returns 'explicit' if CG file has CA beads (Martini 3 explicit backbone
    N/CA/C/O), or 'bb' if it has BB beads (classic Martini 3 single bead).
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(cg_file)
        if len(u.select_atoms("name CA")) > 0:
            return 'explicit'
        if len(u.select_atoms("name BB")) > 0:
            return 'bb'
    except Exception:
        pass
    return 'explicit'   # safe default for nextgen Martini 3


def _backbone_com_per_residue(pdb_file):
    """
    Mass-weighted COM of backbone heavy atoms (N, CA, C, O) per residue.
    Matches CG2AT RMSD_measure_de_novo() exactly.
    Returns (N,3) array in Å.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(pdb_file)
    protein = u.select_atoms("protein")
    coords  = []
    for res in protein.residues:
        bb = res.atoms.select_atoms("name N CA C O")
        if len(bb) == 0:
            bb = res.atoms.select_atoms("name CA")
        if len(bb) == 0:
            continue
        coords.append(bb.center_of_mass())
    return np.array(coords)


def _cg_backbone_positions(cg_file, bb_type):
    """
    Returns (N,3) array of CG backbone positions in Å.
    explicit: CA bead positions
    bb      : BB bead positions
    MDAnalysis auto-converts nm→Å, so no manual scaling needed.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(cg_file)
    sel_name = "CA" if bb_type == 'explicit' else "BB"
    beads    = u.select_atoms(f"name {sel_name}")
    if len(beads) == 0:
        # Fallback: try protein selection
        beads = u.select_atoms(f"protein and name {sel_name}")
    return beads.positions.copy()


def _kabsch_rmsd(pos_mobile, pos_ref):
    """
    Standard backbone RMSD (÷N) after Kabsch optimal alignment.
    Matches the corrected CG2AT Calculate_RMSD with sum over axis=1.
    pos_mobile, pos_ref: (N,3) float arrays, already centred or not.
    """
    c1 = pos_mobile - pos_mobile.mean(axis=0)
    c2 = pos_ref    - pos_ref.mean(axis=0)
    H         = c1.T @ c2
    U, S, Vt  = np.linalg.svd(H)
    d         = np.linalg.det(Vt.T @ U.T)
    D         = np.diag([1.0, 1.0, d])
    R         = Vt.T @ D @ U.T
    rotated   = c1 @ R.T
    diff      = rotated - c2
    # Standard RMSD: sum squared distances per atom, mean over N atoms
    return float(np.round(np.sqrt(np.mean(np.sum(diff**2, axis=1))), 3))


def calculate_rmsd_ca(pdb1, pdb2):
    """
    Pairwise Cα RMSD between two AA structures after Kabsch alignment.
    Standard ÷N formula.
    """
    if not pdb1 or not pdb2:
        return np.nan
    if not os.path.exists(pdb1) or not os.path.exists(pdb2):
        return np.nan
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u1 = mda.Universe(pdb1)
            u2 = mda.Universe(pdb2)
        ca1 = u1.select_atoms("protein and name CA")
        ca2 = u2.select_atoms("protein and name CA")
        if len(ca1) == 0 or len(ca1) != len(ca2):
            return np.nan
        return _kabsch_rmsd(ca2.positions.copy(), ca1.positions.copy())
    except Exception:
        return np.nan


def calculate_rmsd_vs_cg(aa_pdb, cg_file):
    """
    Backbone COM (N/CA/C/O, mass-weighted) RMSD of AA structure vs CG backbone.
    Standard ÷N formula — directly comparable to corrected CG2AT output.
    MDAnalysis reads GRO/PDB in Å automatically (no manual nm→Å conversion).
    """
    if not aa_pdb or not cg_file:
        return np.nan
    if not os.path.exists(aa_pdb) or not os.path.exists(cg_file):
        return np.nan
    try:
        bb_type = _detect_cg_backbone_type(cg_file)
        pos_aa  = _backbone_com_per_residue(aa_pdb)
        pos_cg  = _cg_backbone_positions(cg_file, bb_type)
        n = min(len(pos_aa), len(pos_cg))
        if n == 0:
            return np.nan
        return _kabsch_rmsd(pos_aa[:n], pos_cg[:n])
    except Exception:
        return np.nan


def print_rmsd_report(de_novo_pdb, aligned_pdb, restrained_pdb, cg_ref=None):
    """
    Print two-section backbone RMSD report.

    Section 1 — Backbone COM vs CG reference (only when cg_ref is provided)
    Section 2 — Pairwise Cα RMSD between De novo / Aligned / Restrained

    Returns rmsd_data dict for write_quality_dat().
    """
    W = 100
    print(f"\n{'─'*W}")
    print("  Backbone RMSD report  (CG2AT methodology)")
    print(f"{'─'*W}")

    has_aligned = aligned_pdb and os.path.exists(aligned_pdb)

    rmsd_data = {
        'rmsd_dn_cg'  : np.nan,
        'rmsd_al_cg'  : np.nan,
        'rmsd_rest_cg': np.nan,
        'rmsd_dn_rest': np.nan,
        'rmsd_al_rest': np.nan,
        'rmsd_dn_al'  : np.nan,
    }

    # ── Section 1: vs CG reference ────────────────────────────────────────────
    if cg_ref and os.path.exists(cg_ref):
        bb_type = _detect_cg_backbone_type(cg_ref)
        bb_label = 'explicit backbone N/CA/C/O' if bb_type == 'explicit' else 'classic BB bead'
        print(f"\n  CG reference : {os.path.basename(cg_ref)}  ({bb_label})")
        print( "  AA method    : backbone COM (N/CA/C/O, mass-weighted) per residue")
        print( "  Alignment    : Kabsch optimal rotation")

        rmsd_data['rmsd_dn_cg']   = calculate_rmsd_vs_cg(de_novo_pdb,    cg_ref)
        rmsd_data['rmsd_rest_cg'] = calculate_rmsd_vs_cg(restrained_pdb, cg_ref)
        if has_aligned:
            rmsd_data['rmsd_al_cg'] = calculate_rmsd_vs_cg(aligned_pdb,  cg_ref)

        col_w = 18
        print(f"\n  Backbone COM RMSD vs CG:\n")
        hdr  = f"   {'chain':^7}  {'De novo (Å)':^{col_w}}"
        sep2 = f"   {'-----':^7}  {'------------------':^{col_w}}"
        if has_aligned:
            hdr  += f"  {'Aligned (Å)':^{col_w}}"
            sep2 += f"  {'------------------':^{col_w}}"
        hdr  += f"  {'Restrained (Å)':^{col_w}}"
        sep2 += f"  {'------------------':^{col_w}}"
        print(hdr)
        print(sep2)

        row = f"   {'0':^7}  {rmsd_data['rmsd_dn_cg']:^{col_w}.3f}"
        if has_aligned:
            row += f"  {rmsd_data['rmsd_al_cg']:^{col_w}.3f}"
        row += f"  {rmsd_data['rmsd_rest_cg']:^{col_w}.3f}"
        print(row)
    else:
        if cg_ref:
            print(f"\n  WARNING: CG reference not found: {cg_ref}")
        else:
            print("\n  No CG reference provided — skipping Section 1.")
        print("  (Pass --cg-ref <file.gro> or place CG_INPUT.pdb in INPUT/ folder)")

    # ── Section 2: pairwise Cα RMSD ───────────────────────────────────────────
    print(f"\n  Pairwise Cα RMSD between converted structures:\n")

    rmsd_data['rmsd_dn_rest'] = calculate_rmsd_ca(de_novo_pdb,   restrained_pdb)
    if has_aligned:
        rmsd_data['rmsd_dn_al']   = calculate_rmsd_ca(de_novo_pdb,  aligned_pdb)
        rmsd_data['rmsd_al_rest'] = calculate_rmsd_ca(aligned_pdb, restrained_pdb)

    col_w = 18
    if has_aligned:
        print(f"  {'':22}  {'De novo':^{col_w}}  {'Aligned':^{col_w}}  {'Restrained':^{col_w}}")
        print(f"  {'De novo':22}  {'0.000':^{col_w}}  "
              f"{rmsd_data['rmsd_dn_al']:^{col_w}.3f}  "
              f"{rmsd_data['rmsd_dn_rest']:^{col_w}.3f}")
        print(f"  {'Aligned':22}  {rmsd_data['rmsd_dn_al']:^{col_w}.3f}  "
              f"{'0.000':^{col_w}}  "
              f"{rmsd_data['rmsd_al_rest']:^{col_w}.3f}")
        print(f"  {'Restrained':22}  {rmsd_data['rmsd_dn_rest']:^{col_w}.3f}  "
              f"{rmsd_data['rmsd_al_rest']:^{col_w}.3f}  "
              f"{'0.000':^{col_w}}")
    else:
        print(f"  {'':22}  {'De novo':^{col_w}}  {'Restrained':^{col_w}}")
        print(f"  {'De novo':22}  {'0.000':^{col_w}}  "
              f"{rmsd_data['rmsd_dn_rest']:^{col_w}.3f}")
        print(f"  {'Restrained':22}  {rmsd_data['rmsd_dn_rest']:^{col_w}.3f}  "
              f"{'0.000':^{col_w}}")

    print(f"\n{'─'*W}")
    return rmsd_data


# ══════════════════════════════════════════════════════════════════════════════
# 12.  TDC B-FACTOR PDB
# ══════════════════════════════════════════════════════════════════════════════

def write_tdc_bfactor_pdb(restrained_pdb, per_res_TDC, output_pdb):
    """
    Write a copy of the restrained structure with TDC scores in B-factor column.
    Analogous to AlphaFold pLDDT B-factor colouring.

    Visualise in PyMOL:
        load final_cg2at_restrained_TDC.pdb
        spectrum b, blue_white_red, minimum=0, maximum=100

    Visualise in VMD:
        mol new final_cg2at_restrained_TDC.pdb
        mol modcolor 0 0 Beta
        color scale method BWR
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(restrained_pdb)
        protein  = u.select_atoms("protein")
        residues = protein.residues
        cg_resids = sorted(per_res_TDC.keys())

        tdc_per_atom = {}
        for i, res in enumerate(residues):
            cg_resid = cg_resids[i] if i < len(cg_resids) else None
            tdc      = per_res_TDC.get(cg_resid, 0.0) if cg_resid else 0.0
            for atom in res.atoms:
                tdc_per_atom[atom.index] = tdc

        lines_out = [
            'REMARK  TDC (Trajectory Dihedral Confidence) in B-factor column\n',
            'REMARK  Analogous to pLDDT: 0=low confidence, 100=high confidence\n',
            'REMARK  TDC > 70  HIGH    strong dihedral restraint applied\n',
            'REMARK  TDC 40-70 MEDIUM  moderate restraint applied\n',
            'REMARK  TDC < 40  LOW     no restraint applied\n',
            'REMARK\n',
            'REMARK  Visualise in PyMOL:\n',
            'REMARK    spectrum b, blue_white_red, minimum=0, maximum=100\n',
            'REMARK\n',
            'REMARK  Visualise in VMD:\n',
            'REMARK    mol modcolor 0 0 Beta  |  color scale method BWR\n',
        ]

        with open(restrained_pdb, 'r') as f:
            atom_serial = 0
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        atom_idx = protein.atoms[atom_serial].index
                        tdc      = tdc_per_atom.get(atom_idx, 0.0)
                        new_line = line[:60] + f'{tdc:6.2f}' + line[66:]
                        lines_out.append(new_line)
                        atom_serial += 1
                    except (IndexError, Exception):
                        lines_out.append(line)
                else:
                    lines_out.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(lines_out)

        print(f"  TDC B-factor PDB   → {os.path.basename(output_pdb)}")
    except Exception as e:
        print(f"  WARNING: could not write TDC PDB: {e}")


# ══════════════════════════════════════════════════════════════════════════════
# 13.  RAMACHANDRAN REFERENCE DATA
# ══════════════════════════════════════════════════════════════════════════════

_RAMA500_SETTINGS = {
    'general': {
        'filename': 'rama500-general.data',
        'levels'  : [0.0, 0.0005, 0.02],
        'colors'  : ['#FFFFFF', '#B3E8FF', '#7FD9FF'],
    },
    'gly': {
        'filename': 'rama500-gly-sym.data',
        'levels'  : [0.0, 0.002, 0.02],
        'colors'  : ['#FFFFFF', '#FFE8C5', '#FFCC7F'],
    },
    'pro': {
        'filename': 'rama500-pro.data',
        'levels'  : [0.0, 0.002, 0.02],
        'colors'  : ['#FFFFFF', '#D0FFC5', '#7FFF8C'],
    },
}


def _find_rama500_dir():
    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir / 'rama_data',
        script_dir.parent / 'rama_data',
        Path.home() / '.cg2at_refine' / 'rama_data',
        Path.cwd() / 'rama_data',
    ]
    for d in candidates:
        if d.is_dir() and (d / 'rama500-gly-sym.data').exists():
            return str(d)
    return None


def _load_rama500_file(filepath):
    mid_points = np.arange(-179, 180, 2, dtype=float)
    n = len(mid_points)
    Z = np.zeros((n, n), dtype=float)
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                phi   = float(parts[0])
                psi   = float(parts[1])
                value = float(parts[2])
                i_phi = int((phi + 179) / 2)
                i_psi = int((psi + 179) / 2)
                if 0 <= i_phi < n and 0 <= i_psi < n:
                    Z[i_psi, i_phi] = value
            except (ValueError, IndexError):
                continue
    return Z, mid_points


def _load_mda_reference():
    try:
        import MDAnalysis.analysis.data as mda_data
        data_dir = Path(mda_data.__file__).parent
        for pattern in ['*.npy', '*.npz']:
            for fname in data_dir.rglob(pattern):
                if 'rama' in fname.name.lower():
                    try:
                        data = np.load(str(fname))
                        Z = data if isinstance(data, np.ndarray) \
                            else data[list(data.keys())[0]]
                        if Z.ndim == 2:
                            n    = Z.shape[0]
                            bins = np.linspace(-180, 180, n + 1)
                            ctrs = (bins[:-1] + bins[1:]) / 2
                            return Z, ctrs, ctrs
                    except Exception:
                        continue
    except Exception:
        pass
    return None, None, None


def _classify_mda(phi, psi, Z, phi_bins, psi_bins):
    phi_idx = max(0, min(int(np.searchsorted(phi_bins, phi) - 1), Z.shape[1] - 1))
    psi_idx = max(0, min(int(np.searchsorted(psi_bins, psi) - 1), Z.shape[0] - 1))
    z_val   = Z[psi_idx, phi_idx]
    if z_val >= MDA_LEVEL_FAV:
        return 'favoured'
    elif z_val >= MDA_LEVEL_ALL:
        return 'allowed'
    return 'outlier'


def _classify_rama500(phi, psi, Z, mid_points, levels):
    lower_thresh = levels[1]
    upper_thresh = levels[2]
    i_phi = max(0, min(int(np.searchsorted(mid_points, phi) - 1), len(mid_points) - 1))
    i_psi = max(0, min(int(np.searchsorted(mid_points, psi) - 1), len(mid_points) - 1))
    z_val = Z[i_psi, i_phi]
    if z_val >= upper_thresh:
        return 'favoured'
    elif z_val >= lower_thresh:
        return 'allowed'
    return 'outlier'


def _classify_analytical_general(phi, psi):
    if -160 < phi < -20 and -120 < psi < 50:
        return 'favoured'
    if -180 < phi < -40 and (60 < psi < 180 or -180 < psi < -160):
        return 'favoured'
    if 20 < phi < 90 and 20 < psi < 90:
        return 'favoured'
    if phi < 0:
        return 'allowed'
    if 0 < phi < 120 and -20 < psi < 120:
        return 'allowed'
    return 'outlier'


def _classify_analytical_gly(phi, psi):
    phi_abs = abs(phi)
    if phi_abs < 160 and -120 < psi < 50:
        return 'favoured'
    if phi_abs < 180 and (60 < psi < 180 or -180 < psi < -160):
        return 'favoured'
    return 'allowed'


def _classify_analytical_pro(phi, psi):
    if -90 < phi < -40:
        if -80 < psi < 180 or -180 < psi < -120:
            return 'favoured'
        return 'allowed'
    if -120 < phi < -40:
        return 'allowed'
    return 'outlier'


# ══════════════════════════════════════════════════════════════════════════════
# 14.  RAMACHANDRAN QUALITY ASSESSMENT
# ══════════════════════════════════════════════════════════════════════════════

def _get_phi_psi_all(pdb_file):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u        = mda.Universe(pdb_file)
    protein  = u.select_atoms("protein")
    residues = protein.residues
    n_res    = len(residues)
    result   = []
    for i, res in enumerate(residues):
        if i == 0 or i == n_res - 1:
            continue
        try:
            prev = residues[i - 1]
            nxt  = residues[i + 1]
            C_prev  = prev.atoms.select_atoms("name C")[0].position
            N_curr  = res.atoms.select_atoms("name N")[0].position
            CA_curr = res.atoms.select_atoms("name CA")[0].position
            C_curr  = res.atoms.select_atoms("name C")[0].position
            N_next  = nxt.atoms.select_atoms("name N")[0].position
            phi = np.degrees(calc_dihedrals(
                C_prev[np.newaxis], N_curr[np.newaxis],
                CA_curr[np.newaxis], C_curr[np.newaxis])[0])
            psi = np.degrees(calc_dihedrals(
                N_curr[np.newaxis], CA_curr[np.newaxis],
                C_curr[np.newaxis], N_next[np.newaxis])[0])
            result.append((f"{res.resname}{res.resid}", phi, psi, res.resname))
        except Exception:
            pass
    return result


def _make_stats(favoured, allowed, outliers, gly_pro=None):
    total = len(favoured) + len(allowed) + len(outliers)
    return {
        'favoured' : favoured,
        'allowed'  : allowed,
        'outliers' : outliers,
        'gly_pro'  : gly_pro or [],
        'total'    : total,
        'pct_fav'  : 100.0 * len(favoured) / total if total > 0 else 0.0,
        'pct_all'  : 100.0 * len(allowed)  / total if total > 0 else 0.0,
        'pct_out'  : 100.0 * len(outliers) / total if total > 0 else 0.0,
        'n_gly_pro': len(gly_pro) if gly_pro else 0,
    }


def assess_ramachandran(pdb_file, Z=None, phi_bins=None, psi_bins=None,
                         rama500_general=None):
    if not pdb_file or not os.path.exists(pdb_file):
        return None
    all_res  = _get_phi_psi_all(pdb_file)
    favoured, allowed, outliers, gly_pro = [], [], [], []
    for label, phi, psi, resname in all_res:
        if resname in SPECIAL_RES:
            gly_pro.append((label, phi, psi, resname))
            continue
        if rama500_general is not None:
            Z_g, mp_g, lv_g = rama500_general
            region = _classify_rama500(phi, psi, Z_g, mp_g, lv_g)
        elif Z is not None:
            region = _classify_mda(phi, psi, Z, phi_bins, psi_bins)
        else:
            region = _classify_analytical_general(phi, psi)
        if region == 'favoured':
            favoured.append((label, phi, psi))
        elif region == 'allowed':
            allowed.append((label, phi, psi))
        else:
            outliers.append((label, phi, psi))
    return _make_stats(favoured, allowed, outliers, gly_pro)


def assess_ramachandran_gly(pdb_file, rama500_gly=None):
    if not pdb_file or not os.path.exists(pdb_file):
        return None
    all_res  = _get_phi_psi_all(pdb_file)
    favoured, allowed, outliers = [], [], []
    for label, phi, psi, resname in all_res:
        if resname != 'GLY':
            continue
        if rama500_gly is not None:
            Z_g, mp_g, lv_g = rama500_gly
            region = _classify_rama500(phi, psi, Z_g, mp_g, lv_g)
        else:
            region = _classify_analytical_gly(phi, psi)
        if region == 'favoured':
            favoured.append((label, phi, psi))
        elif region == 'allowed':
            allowed.append((label, phi, psi))
        else:
            outliers.append((label, phi, psi))
    return _make_stats(favoured, allowed, outliers)


def assess_ramachandran_pro(pdb_file, rama500_pro=None):
    if not pdb_file or not os.path.exists(pdb_file):
        return None
    all_res  = _get_phi_psi_all(pdb_file)
    favoured, allowed, outliers = [], [], []
    for label, phi, psi, resname in all_res:
        if resname != 'PRO':
            continue
        if rama500_pro is not None:
            Z_g, mp_g, lv_g = rama500_pro
            region = _classify_rama500(phi, psi, Z_g, mp_g, lv_g)
        else:
            region = _classify_analytical_pro(phi, psi)
        if region == 'favoured':
            favoured.append((label, phi, psi))
        elif region == 'allowed':
            allowed.append((label, phi, psi))
        else:
            outliers.append((label, phi, psi))
    return _make_stats(favoured, allowed, outliers)


# ══════════════════════════════════════════════════════════════════════════════
# 15.  PRINT / WRITE QUALITY SUMMARIES
# ══════════════════════════════════════════════════════════════════════════════

def print_rama_summary(label, stats):
    if stats is None:
        return
    w = 60
    print(f"\n  {'─' * w}")
    print(f"  Ramachandran — {label}")
    print(f"  {'─' * w}")
    print(f"  Favoured  : {stats['pct_fav']:>6.2f}%  "
          f"({len(stats['favoured'])}/{stats['total']})")
    print(f"  Allowed   : {stats['pct_all']:>6.2f}%  "
          f"({len(stats['allowed'])}/{stats['total']})")
    print(f"  Outliers  : {stats['pct_out']:>6.2f}%  "
          f"({len(stats['outliers'])}/{stats['total']})")
    if stats['outliers']:
        names = ', '.join(o[0] for o in stats['outliers'])
        print(f"  Outlier residues: {names}")
    if stats.get('n_gly_pro', 0) > 0:
        print(f"  (GLY/PRO excluded: {stats['n_gly_pro']} residues)")
    print(f"  {'─' * w}")


def write_quality_dat(output_path, stats_dn, stats_rest, stats_al=None,
                       stats_gly_dn=None, stats_gly_rest=None, stats_gly_al=None,
                       stats_pro_dn=None, stats_pro_rest=None, stats_pro_al=None,
                       rmsd_data=None):
    """
    Write structure_quality_restrained.dat.
    Includes Ramachandran statistics (General, GLY, PRO) and
    backbone RMSD report matching CG2AT format.
    """
    has_al = stats_al is not None

    def make_table(s_dn, s_rest, s_al):
        if s_dn is None or s_rest is None:
            return ['  (no data)\n']
        hdr = f'{"Metric":<35} {"De novo":>12} {"Restrained":>12}'
        sep = f'{"─"*35} {"─"*12} {"─"*12}'
        if s_al:
            hdr += f' {"Aligned":>12}'
            sep += f' {"─"*12}'
        lines = [hdr + '\n', sep + '\n']

        def row(label, key, fmt='.2f'):
            dn  = format(s_dn[key],   fmt)
            rs  = format(s_rest[key], fmt)
            s   = f"{label:<35} {dn:>12} {rs:>12}"
            if s_al:
                al = format(s_al[key], fmt)
                s += f" {al:>12}"
            return s + '\n'

        lines.append(row('Ramachandran favoured (%)', 'pct_fav'))
        lines.append(row('Ramachandran allowed  (%)', 'pct_all'))
        lines.append(row('Ramachandran outliers (%)', 'pct_out'))
        lines.append(row('N residues assessed',       'total',   'd'))
        lines.append('\nOutlier residues:\n')
        lines.append(f'  De novo    : '
                     f'{", ".join(o[0] for o in s_dn["outliers"]) or "none"}\n')
        lines.append(f'  Restrained : '
                     f'{", ".join(o[0] for o in s_rest["outliers"]) or "none"}\n')
        if s_al:
            lines.append(f'  Aligned    : '
                         f'{", ".join(o[0] for o in s_al["outliers"]) or "none"}\n')
        return lines

    sep_line = f'\n{"─"*60}\n'
    all_lines = [
        'Structure quality assessment — cg2at_refine.py (TDC method)\n',
        f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n',
    ]

    # ── Ramachandran sections ─────────────────────────────────────────────────
    all_lines += [sep_line, 'GENERAL RESIDUES (non-GLY, non-PRO)\n', sep_line]
    all_lines += make_table(stats_dn, stats_rest, stats_al if has_al else None)
    if stats_dn:
        all_lines.append(f'\n  GLY/PRO excluded: '
                         f'{stats_dn.get("n_gly_pro", 0)} residues\n')

    if stats_gly_dn is not None:
        all_lines += [sep_line, 'GLYCINE\n', sep_line]
        all_lines += make_table(stats_gly_dn, stats_gly_rest,
                                 stats_gly_al if has_al else None)

    if stats_pro_dn is not None:
        all_lines += [sep_line, 'PROLINE\n', sep_line]
        all_lines += make_table(stats_pro_dn, stats_pro_rest,
                                 stats_pro_al if has_al else None)

    # ── RMSD section ─────────────────────────────────────────────────────────
    if rmsd_data:
        has_cg  = not np.isnan(rmsd_data.get('rmsd_dn_cg',   np.nan))
        has_al  = not np.isnan(rmsd_data.get('rmsd_dn_al',   np.nan))
        has_alr = not np.isnan(rmsd_data.get('rmsd_al_rest', np.nan))

        all_lines += [sep_line, 'BACKBONE RMSD (standard ÷N formula, Kabsch alignment)\n', sep_line]
        all_lines.append(f'{"Metric":<42} {"Value (Å)":>10}\n')
        all_lines.append(f'{"─"*42} {"─"*10}\n')

        if has_cg:
            all_lines.append(f'{"Backbone COM RMSD: De novo vs CG":<42} {rmsd_data["rmsd_dn_cg"]:>10.3f}\n')
            all_lines.append(f'{"Backbone COM RMSD: Restrained vs CG":<42} {rmsd_data["rmsd_rest_cg"]:>10.3f}\n')
            if has_al:
                all_lines.append(f'{"Backbone COM RMSD: Aligned vs CG":<42} {rmsd_data["rmsd_al_cg"]:>10.3f}\n')
            all_lines.append('\n')

        all_lines.append(f'{"Cα RMSD: De novo vs Restrained":<42} {rmsd_data["rmsd_dn_rest"]:>10.3f}\n')
        if has_al:
            all_lines.append(f'{"Cα RMSD: De novo vs Aligned":<42} {rmsd_data["rmsd_dn_al"]:>10.3f}\n')
        if has_alr:
            all_lines.append(f'{"Cα RMSD: Aligned vs Restrained":<42} {rmsd_data["rmsd_al_rest"]:>10.3f}\n')

        all_lines.append('\nNote: standard RMSD (÷N, one value per residue/atom).\n'
                         '      Kabsch optimal rotation applied before RMSD.\n')
    with open(output_path, 'w') as f:
        f.writelines(all_lines)


# ══════════════════════════════════════════════════════════════════════════════
# 16.  SHARED PLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _style_ax(ax, title, is_first_col=False, row_label=None):
    ax.set_xticks(PLOT_TICKS)
    ax.set_xticklabels(PLOT_TICK_LABELS, fontsize=PLOT_TICK_FONTSIZE)
    ax.set_xlabel('Φ (°)', fontsize=PLOT_LABEL_FONTSIZE)
    ax.set_yticks(PLOT_TICKS)
    ax.set_yticklabels(PLOT_TICK_LABELS, fontsize=PLOT_TICK_FONTSIZE)
    ax.tick_params(labelleft=True)
    if is_first_col:
        if row_label:
            ax.set_ylabel(f'{row_label}\nΨ (°)', fontsize=PLOT_LABEL_FONTSIZE)
        else:
            ax.set_ylabel('Ψ (°)', fontsize=PLOT_LABEL_FONTSIZE)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_title(title, fontsize=PLOT_TITLE_FONTSIZE, fontweight='bold', pad=6)


def _add_stats_table(ax, stats):
    if stats is None:
        return
    table = ax.table(
        cellText=[[f"{stats['pct_fav']:.1f}%",
                   f"{stats['pct_all']:.1f}%",
                   f"{stats['pct_out']:.1f}%"]],
        colLabels=['Favoured', 'Allowed', 'Outlier'],
        cellLoc='center', colLoc='center',
        bbox=[0.10, -0.26, 0.80, 0.11],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(PLOT_TABLE_FONTSIZE)
    hdr_colors = [COLOR_FAV, COLOR_ALL, COLOR_OUT]
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor('0.70')
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor(hdr_colors[c])
            cell.get_text().set_color('white')
            cell.get_text().set_weight('bold')
            cell.get_text().set_fontsize(PLOT_TABLE_FONTSIZE)
        else:
            cell.set_facecolor('white')
            cell.get_text().set_weight('bold')
            cell.get_text().set_fontsize(PLOT_TABLE_VAL_FONTSIZE)


def _scatter_classified(ax, favoured, allowed, outliers):
    kw = dict(alpha=DOT_ALPHA, zorder=5, edgecolors='white', linewidths=0.3)
    if favoured:
        ax.scatter([p for _, p, _ in favoured], [s for _, _, s in favoured],
                   c=COLOR_FAV, s=DOT_SIZE, **kw)
    if allowed:
        ax.scatter([p for _, p, _ in allowed],  [s for _, _, s in allowed],
                   c=COLOR_ALL, s=DOT_SIZE, **kw)
    if outliers:
        ax.scatter([p for _, p, _ in outliers], [s for _, _, s in outliers],
                   c=COLOR_OUT, s=DOT_SIZE + 10,
                   alpha=DOT_ALPHA, zorder=6, edgecolors='white', linewidths=0.5)


# ══════════════════════════════════════════════════════════════════════════════
# 17.  RAMACHANDRAN PLOT — GENERAL
# ══════════════════════════════════════════════════════════════════════════════

def plot_ramachandran(de_novo_pdb, restrained_pdb, aligned_pdb, output_png,
                       Z=None, phi_bins=None, psi_bins=None,
                       rama500_general_data=None):
    try:
        from MDAnalysis.analysis.dihedrals import Ramachandran as MDArama
    except ImportError:
        print("  WARNING: MDAnalysis Ramachandran not available — skipping plot")
        return

    pdbs = [(de_novo_pdb, 'De novo'), (restrained_pdb, 'TDC Restrained')]
    if aligned_pdb and os.path.exists(aligned_pdb):
        pdbs.append((aligned_pdb, 'Aligned (-a flag)'))

    n = len(pdbs)
    fig, axes = plt.subplots(1, n,
                              figsize=(PLOT_FIGSIZE_PER_PANEL[0] * n,
                                       PLOT_FIGSIZE_PER_PANEL[1]),
                              sharey=True)
    if n == 1:
        axes = [axes]

    for col_idx, (ax, (pdb, title)) in enumerate(zip(axes, pdbs)):
        is_first = (col_idx == 0)
        if not os.path.exists(pdb):
            ax.set_title(f'{title}\n(not found)', fontsize=PLOT_TITLE_FONTSIZE)
            continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u   = mda.Universe(pdb)
                sel = u.select_atoms("protein")
                R   = MDArama(sel).run()
                R.plot(ax=ax, ref=True, color='black', s=0)
            for line in ax.lines:
                x, y = line.get_xdata(), line.get_ydata()
                if ((len(y) == 2 and np.allclose(y, [0, 0])) or
                        (len(x) == 2 and np.allclose(x, [0, 0]))):
                    line.set_color('0.60')
                    line.set_linestyle('--')
                    line.set_linewidth(0.8)
        except Exception:
            ax.set_facecolor('#F0F0F0')
            ax.axhline(0, color='0.60', linestyle='--', linewidth=0.8)
            ax.axvline(0, color='0.60', linestyle='--', linewidth=0.8)

        all_res   = _get_phi_psi_all(pdb)
        n_gly_pro = sum(1 for _, _, _, rn in all_res if rn in SPECIAL_RES)
        favoured, allowed, outliers = [], [], []

        for label, phi, psi, resname in all_res:
            if resname in SPECIAL_RES:
                continue
            if rama500_general_data is not None:
                Z_g, mp_g, lv_g = rama500_general_data
                region = _classify_rama500(phi, psi, Z_g, mp_g, lv_g)
            elif Z is not None:
                region = _classify_mda(phi, psi, Z, phi_bins, psi_bins)
            else:
                region = _classify_analytical_general(phi, psi)
            if region == 'favoured':
                favoured.append((label, phi, psi))
            elif region == 'allowed':
                allowed.append((label, phi, psi))
            else:
                outliers.append((label, phi, psi))

        _scatter_classified(ax, favoured, allowed, outliers)
        _style_ax(ax, title, is_first_col=is_first)

        total = len(favoured) + len(allowed) + len(outliers)
        _add_stats_table(ax, {
            'pct_fav': 100.0 * len(favoured) / total if total > 0 else 0.0,
            'pct_all': 100.0 * len(allowed)  / total if total > 0 else 0.0,
            'pct_out': 100.0 * len(outliers) / total if total > 0 else 0.0,
            'favoured': favoured, 'allowed': allowed, 'outliers': outliers,
        })

        ax.text(0.5, -0.33, f'* {n_gly_pro} GLY/PRO excluded (separate plot)',
                transform=ax.transAxes, ha='center',
                fontsize=7, color='0.50', style='italic')

    fig.suptitle('Ramachandran Plot — General Residues (GLY/PRO excluded)',
                 fontsize=PLOT_SUPTITLE_FONTSIZE, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0.06, 1, 1.0])
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  General Ramachandran → {os.path.basename(output_png)}")


# ══════════════════════════════════════════════════════════════════════════════
# 18.  RAMACHANDRAN PLOT — GLY / PRO
# ══════════════════════════════════════════════════════════════════════════════

def plot_ramachandran_gly_pro(de_novo_pdb, restrained_pdb, aligned_pdb,
                               output_png,
                               rama500_gly_data=None, rama500_pro_data=None):
    pdbs = [(de_novo_pdb, 'De novo'), (restrained_pdb, 'TDC Restrained')]
    if aligned_pdb and os.path.exists(aligned_pdb):
        pdbs.append((aligned_pdb, 'Aligned (-a flag)'))

    n_cols = len(pdbs)

    row_configs = [
        {
            'resname'      : 'GLY',
            'row_label'    : 'Glycine',
            'rama500_data' : rama500_gly_data,
            'settings'     : _RAMA500_SETTINGS['gly'],
            'classify_fn'  : lambda phi, psi: (
                _classify_rama500(phi, psi, *rama500_gly_data)
                if rama500_gly_data else _classify_analytical_gly(phi, psi)
            ),
        },
        {
            'resname'      : 'PRO',
            'row_label'    : 'Proline',
            'rama500_data' : rama500_pro_data,
            'settings'     : _RAMA500_SETTINGS['pro'],
            'classify_fn'  : lambda phi, psi: (
                _classify_rama500(phi, psi, *rama500_pro_data)
                if rama500_pro_data else _classify_analytical_pro(phi, psi)
            ),
        },
    ]

    fig, axes = plt.subplots(
        2, n_cols,
        figsize=(PLOT_FIGSIZE_PER_PANEL[0] * n_cols,
                 PLOT_FIGSIZE_PER_PANEL[1] * 2 + 0.5),
    )
    if n_cols == 1:
        axes = [[axes[0]], [axes[1]]]

    for row_idx, cfg in enumerate(row_configs):
        for col_idx, (pdb, col_title) in enumerate(pdbs):
            ax       = axes[row_idx][col_idx]
            is_first = (col_idx == 0)

            if cfg['rama500_data'] is not None:
                Z_r, mp_r, lv_r = cfg['rama500_data']
                settings         = cfg['settings']
                levels           = settings['levels'] + [Z_r.max() + 1]
                ax.contourf(mp_r, mp_r, Z_r, levels=levels,
                            colors=settings['colors'])
                ax.contour(mp_r, mp_r, Z_r, levels=settings['levels'][1:],
                           colors=['0.5'], linewidths=0.5, alpha=0.5)
            else:
                ax.set_facecolor('#F8F8F8')

            ax.axhline(0, color='0.60', linestyle='--', linewidth=0.8, zorder=3)
            ax.axvline(0, color='0.60', linestyle='--', linewidth=0.8, zorder=3)

            if not os.path.exists(pdb):
                ax.set_title(f'{col_title}\n(not found)', fontsize=PLOT_TITLE_FONTSIZE)
                continue

            all_res  = _get_phi_psi_all(pdb)
            favoured, allowed, outliers = [], [], []

            for label, phi, psi, resname in all_res:
                if resname != cfg['resname']:
                    continue
                region = cfg['classify_fn'](phi, psi)
                if region == 'favoured':
                    favoured.append((label, phi, psi))
                elif region == 'allowed':
                    allowed.append((label, phi, psi))
                else:
                    outliers.append((label, phi, psi))

            _scatter_classified(ax, favoured, allowed, outliers)
            _style_ax(ax, col_title, is_first_col=is_first,
                      row_label=cfg['row_label'])

            total = len(favoured) + len(allowed) + len(outliers)
            _add_stats_table(ax, {
                'pct_fav' : 100.0 * len(favoured) / total if total > 0 else 0.0,
                'pct_all' : 100.0 * len(allowed)  / total if total > 0 else 0.0,
                'pct_out' : 100.0 * len(outliers) / total if total > 0 else 0.0,
                'favoured': favoured, 'allowed': allowed, 'outliers': outliers,
            })

            src = ('rama500 (Lovell 2003)' if cfg['rama500_data']
                   else 'analytical boundaries (Lovell 2003)')
            ax.text(0.5, -0.33, f'Background: {src}',
                    transform=ax.transAxes, ha='center',
                    fontsize=7, color='0.50', style='italic')

    fig.suptitle('Ramachandran Plot — GLY and PRO residues',
                 fontsize=PLOT_SUPTITLE_FONTSIZE, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0.04, 1, 1.0])
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  GLY/PRO Ramachandran → {os.path.basename(output_png)}")


# ══════════════════════════════════════════════════════════════════════════════
# 19.  ARGUMENT PARSER
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        prog='cg2at_refine.py',
        description=textwrap.dedent("""\
            Post-processing refinement for CG2AT.
            Applies trajectory-informed dihedral restraints (TDC method)
            to improve Ramachandran quality of CG→AA backmapped structures.
            Only possible with the new Martini 3 explicit backbone (N/CA/C/O).
        """),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              # Basic usage
              python cg2at_refine.py -cg2at CG2AT_2026-08-01_17-07-15 \\
                  -tpr protein_CG.tpr -xtc trajectory.xtc -step 10

              # Use aligned structure as starting point
              python cg2at_refine.py -cg2at CG2AT_* \\
                  -tpr protein_CG.tpr -xtc trajectory.xtc --aa-input aligned

              # Fixed dphi (good for flexible proteins)
              python cg2at_refine.py -cg2at CG2AT_* \\
                  -tpr protein_CG.tpr -xtc trajectory.xtc --dphi 10.0

              # Use portion of trajectory
              python cg2at_refine.py -cg2at CG2AT_* \\
                  -tpr protein_CG.tpr -xtc trajectory.xtc -b 5000 -e 10000
        """),
    )

    req = p.add_argument_group('Required arguments')
    req.add_argument('-cg2at', required=True, metavar='DIR',
        help='Path to CG2AT output folder (e.g. CG2AT_2026-08-01_17-07-15). '
             'Must contain FINAL/final_cg2at_de_novo.pdb and FINAL/PROTEIN_0.itp.')
    req.add_argument('-tpr', required=True, metavar='FILE',
        help='CG topology file (.tpr or .gro). Used to load the CG trajectory.')
    req.add_argument('-xtc', required=True, metavar='FILE',
        help='CG trajectory file (.xtc). φ/ψ distributions are extracted from this.')

    traj = p.add_argument_group('Trajectory options')
    traj.add_argument('-step', type=int, default=10, metavar='N',
        help='Frame stride for φ/ψ extraction (default: 10). '
             'Use 1 for maximum statistical quality, 10 for speed.')
    traj.add_argument('-b', type=float, default=0.0, metavar='ps',
        help='Start time in ps (default: 0 = beginning). '
             'Use to skip equilibration, e.g. -b 5000.')
    traj.add_argument('-e', type=float, default=-1.0, metavar='ps',
        help='End time in ps (default: -1 = full trajectory).')

    rest = p.add_argument_group('Restraint options')
    rest.add_argument('--dphi', type=float, default=None, metavar='DEG',
        help='Fixed dphi tolerance in degrees (default: adaptive σ×0.3, '
             'clipped 5°–20°). Use --dphi 10.0 for flexible proteins.')
    rest.add_argument('--kmax', type=float, default=K_MAX, metavar='kJ',
        help=f'Maximum force constant kJ/mol/rad² (default: {K_MAX}). '
             'Use 1000 for stronger restraints if improvement is insufficient.')
    rest.add_argument('--sigma', type=float, default=SIG_INFLECTION, metavar='DEG',
        help=f'Sigmoid inflection point in degrees (default: {SIG_INFLECTION}°). '
             'k = kmax × exp(-σ/sigma).')
    rest.add_argument('--tdc-low', type=float, default=TDC_LOW, metavar='TDC',
        help=f'Minimum TDC to apply a restraint (default: {TDC_LOW}). '
             'Residues below this are left to energy minimisation.')

    inp = p.add_argument_group('Input structure options')
    inp.add_argument('--aa-input',
        choices=['de_novo', 'aligned'],
        default='de_novo',
        metavar='STRUCT',
        help='AA structure for restraint template and GROMACS input. '
             '"de_novo" (default): use final_cg2at_de_novo.pdb. '
             '"aligned": use final_cg2at_aligned.pdb — recommended when an '
             'experimental or AF3 reference was provided via CG2AT -a flag.')

    inp.add_argument('--cg-ref', default=None, metavar='FILE',
        help='CG reference GRO/PDB for backbone RMSD (same file used with -c in '
             'CG2AT). If omitted, auto-detected from CG2AT_*/INPUT/CG_INPUT.pdb.')

    gmx_grp = p.add_argument_group('GROMACS options')
    gmx_grp.add_argument('--gmx', default=None, metavar='EXE',
        help='GROMACS executable (default: auto-detected from CG2AT logs).')
    gmx_grp.add_argument('--nsteps-nvt', type=int, default=50000, metavar='N',
        help='NVT MD steps (default: 50000 = 100 ps at dt=0.001 ps).')
    gmx_grp.add_argument('--temp', type=float, default=300.0, metavar='K',
        help='NVT temperature in K (default: 300).')
    gmx_grp.add_argument('--no-nvt', action='store_true',
        help='Skip NVT — minimisation only. NOT recommended for Ramachandran improvement.')
    gmx_grp.add_argument('--ntmpi', type=int, default=1,
        help='MPI threads for mdrun (default: 1).')
    gmx_grp.add_argument('--ntomp', type=int, default=4,
        help='OpenMP threads per MPI rank (default: 4).')

    out_grp = p.add_argument_group('Output options')
    out_grp.add_argument('--keep-tmp', action='store_true',
        help='Keep MDP files in FINAL/ instead of gromacs_outputs_refine/.')
    out_grp.add_argument('--no-plots', action='store_true',
        help='Skip Ramachandran plot generation.')
    out_grp.add_argument('--prefix', default='final_cg2at_restrained',
        help='Output PDB prefix (default: final_cg2at_restrained).')

    return p.parse_args()


# ══════════════════════════════════════════════════════════════════════════════
# 20.  BANNER
# ══════════════════════════════════════════════════════════════════════════════

def print_banner():
    W   = 100
    sep = '-' * W
    lines = [
        sep, '',
        f"{'CG2AT-REFINE':^{W}}",
        f"{'Trajectory-Informed Dihedral Restraint Refinement for CG2AT':^{W}}",
        '',
        f"{'Written by Hafez Razmazma  |  Supervised by Phillip J. Stansfeld':^{W}}",
        f"{'University of Warwick, Coventry, UK':^{W}}",
        f"{'Contact: hafez.razmazma@warwick.ac.uk':^{W}}",
        '',
        f"{'Please cite:':^{W}}",
        f"{'CG2AT2: Vickery & Stansfeld, JCTC 2021, DOI: 10.1021/acs.jctc.1c00295':^{W}}",
        '',
        sep, '',
    ]
    print('\n'.join(lines))


# ══════════════════════════════════════════════════════════════════════════════
# 21.  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    args  = parse_args()
    timer = StepTimer()

    outputs   = find_cg2at_outputs(args.cg2at)
    final_dir = outputs['final_dir']

    log_path   = os.path.join(final_dir, OUT_LOG)
    tee        = Tee(log_path)
    sys.stdout = tee

    print_banner()
    print(f"  Started : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Log     : {log_path}")
    print()
    timer.tick('Initialisation and setup')

    # ── 1 — outputs ────────────────────────────────────────────────────────────
    print("[1/8] Locating CG2AT outputs ...")
    print(f"  CG2AT folder : {outputs['cg2at_folder']}")
    print(f"  FINAL dir    : {os.path.basename(final_dir)}")
    print(f"  De novo PDB  : {os.path.basename(outputs['de_novo_pdb'])}")
    if outputs['aligned_pdb']:
        print(f"  Aligned PDB  : {os.path.basename(outputs['aligned_pdb'])}")
    else:
        print(f"  Aligned PDB  : not found")

    # ── Select AA input structure ──────────────────────────────────────────────
    if args.aa_input == 'aligned':
        if outputs['aligned_pdb'] and os.path.exists(outputs['aligned_pdb']):
            aa_input_pdb = outputs['aligned_pdb']
            print(f"  AA input     : final_cg2at_aligned.pdb  (--aa-input aligned)")
        else:
            print(f"  WARNING: aligned PDB not found — falling back to de novo")
            aa_input_pdb = outputs['de_novo_pdb']
    else:
        aa_input_pdb = outputs['de_novo_pdb']
        print(f"  AA input     : final_cg2at_de_novo.pdb  (default)")

    # ── 2 — GROMACS ────────────────────────────────────────────────────────────
    print("\n[2/8] Detecting GROMACS ...")
    gmx = args.gmx if args.gmx else find_gromacs_executable(outputs['cg2at_folder'])
    print(f"  GMX executable: {gmx}")

    # ── 3 — time range ─────────────────────────────────────────────────────────
    print("\n[3/8] Parsing trajectory time range ...")
    start_frame, stop_frame, n_frames = parse_time_range(
        args.tpr, args.xtc, args.b, args.e, args.step)
    if args.b > 0 or args.e > 0:
        print(f"  Time range : {args.b:.0f} – "
              f"{'end' if args.e < 0 else f'{args.e:.0f}'} ps  |  "
              f"{n_frames} frames (stride {args.step})")
    else:
        print(f"  Using full trajectory ({n_frames} frames, stride {args.step})")

    # ── 4 — phi/psi ────────────────────────────────────────────────────────────
    print("\n[4/8] Extracting CG-backbone-derived φ/ψ from trajectory ...")
    tee.set_log_only(True)
    if _TDG_AVAILABLE:
        distributions = extract_phi_psi(
            args.tpr, args.xtc,
            start=start_frame, stop=stop_frame, step=args.step)
    else:
        distributions = _extract_phi_psi_inline(
            args.tpr, args.xtc,
            start=start_frame, stop=stop_frame, step=args.step)
    tee.set_log_only(False)
    print(f"  Done — {n_frames} frames processed for {len(distributions)} residues.")
    timer.tick('Extract φ/ψ from CG trajectory')

    # ── 5 — TDC ────────────────────────────────────────────────────────────────
    print("\n[5/8] Computing TDC scores ...")
    if _TDG_AVAILABLE:
        result         = compute_TDC(distributions)
        global_TDC     = result[0]
        global_weight  = result[1]
        per_res_TDC    = result[2]
        per_res_weight = result[3]
        per_res_stats  = result[4]
    else:
        global_TDC, global_weight, per_res_TDC, per_res_weight = \
            _compute_TDC_inline(distributions)
        per_res_stats = {}

    n_high   = sum(1 for t in per_res_TDC.values() if t >= 70)
    n_medium = sum(1 for t in per_res_TDC.values() if 40 <= t < 70)
    n_low    = sum(1 for t in per_res_TDC.values() if t <  40)
    print(f"  Global TDC score : {global_TDC:.1f}/100")
    print(f"  Global weight    : {global_weight:.2f}")
    print(f"  HIGH   (TDC≥70)  : {n_high} residues → strong restraints")
    print(f"  MEDIUM (40-70)   : {n_medium} residues → moderate restraints")
    print(f"  LOW    (TDC<40)  : {n_low} residues → not restrained")

    # TDC table → log only
    if _TDG_AVAILABLE and per_res_stats:
        tee.set_log_only(True)
        print_TDC_table(distributions, per_res_TDC, per_res_weight,
                        per_res_stats, global_TDC, global_weight)
        tee.set_log_only(False)

    timer.tick('Compute TDC scores')

    # ── 6 — restraints ─────────────────────────────────────────────────────────
    print("\n[6/8] Generating dihedral restraints ...")
    itp_path = os.path.join(final_dir, OUT_RESTRAINTS_ITP)
    if _TDG_AVAILABLE:
        write_dihedral_restraints(
            distributions, per_res_TDC, per_res_stats,
            aa_pdb     = aa_input_pdb,
            output_itp = itp_path,
            k_max      = args.kmax,
            sigma      = args.sigma,
            tdc_low    = args.tdc_low,
            fixed_dphi = args.dphi)
    else:
        _write_restraints_inline(
            distributions, per_res_TDC,
            aa_pdb=aa_input_pdb, output_itp=itp_path,
            k_max=args.kmax, sigma=args.sigma,
            tdc_low=args.tdc_low, dphi_fixed=args.dphi)
    patch_protein_itp(outputs['protein_itp'], OUT_RESTRAINTS_ITP)
    timer.tick('Generate dihedral restraints ITP')

    # ── 7 — MDP files ──────────────────────────────────────────────────────────
    print("\n[7/8] Writing MDP files ...")
    minim_mdp, nvt_mdp = write_mdp_files(
        final_dir, nsteps_nvt=args.nsteps_nvt, temp=args.temp)

    # ── 8 — GROMACS pipeline ───────────────────────────────────────────────────
    print("\n[8/8] Running GROMACS pipeline ...")
    out_pdb = run_gromacs_pipeline(
        gmx=gmx, final_dir=final_dir,
        aa_input_pdb=aa_input_pdb, topol_top=outputs['topol_top'],
        minim_mdp=minim_mdp, nvt_mdp=nvt_mdp,
        skip_nvt=args.no_nvt, ntmpi=args.ntmpi, ntomp=args.ntomp,
        out_prefix=args.prefix)
    timer.tick('GROMACS minimisation and NVT MD')

    # ── TDC B-factor PDB ───────────────────────────────────────────────────────
    tdc_pdb = os.path.join(final_dir, OUT_TDC_PDB)
    write_tdc_bfactor_pdb(out_pdb, per_res_TDC, tdc_pdb)

    # ── RMSD report ────────────────────────────────────────────────────────────
    # Resolve CG reference:
    #   1. --cg-ref flag (explicit user choice)
    #   2. Auto-detected CG2AT_*/INPUT/CG_INPUT.pdb
    #   3. None — Section 1 skipped, pairwise only
    if args.cg_ref:
        cg_ref = args.cg_ref
        print(f"  CG reference : {os.path.basename(cg_ref)}  (--cg-ref flag)")
    elif outputs.get('cg_input_pdb'):
        cg_ref = outputs['cg_input_pdb']
        print(f"  CG reference : auto-detected INPUT/CG_INPUT.pdb")
    else:
        cg_ref = None
        print("  CG reference : not found — pairwise RMSD only")

    rmsd_data = print_rmsd_report(
        de_novo_pdb    = outputs['de_novo_pdb'],
        aligned_pdb    = outputs['aligned_pdb'],
        restrained_pdb = out_pdb,
        cg_ref         = cg_ref,
    )

    # ── Load reference data ────────────────────────────────────────────────────
    print('\n── Ramachandran quality assessment ' + '─' * 35)

    Z, phi_bins, psi_bins = _load_mda_reference()
    if Z is not None:
        print(f"  MDAnalysis reference loaded ✅")
    else:
        print("  MDAnalysis reference not found — using analytical boundaries")

    rama500_dir          = _find_rama500_dir()
    rama500_general_data = None
    rama500_gly_data     = None
    rama500_pro_data     = None

    if rama500_dir:
        print(f"  rama500 data: {rama500_dir}")
        for key, cfg in _RAMA500_SETTINGS.items():
            fpath = os.path.join(rama500_dir, cfg['filename'])
            if os.path.exists(fpath):
                try:
                    Z_r, mp_r  = _load_rama500_file(fpath)
                    data_tuple = (Z_r, mp_r, cfg['levels'])
                    if key == 'general':
                        rama500_general_data = data_tuple
                    elif key == 'gly':
                        rama500_gly_data = data_tuple
                    elif key == 'pro':
                        rama500_pro_data = data_tuple
                except Exception:
                    pass
    else:
        print("  rama500 not found — GLY/PRO will use analytical boundaries")

    # ── Assessment ─────────────────────────────────────────────────────────────
    def assess(pdb):
        if not pdb or not os.path.exists(pdb):
            return None, None, None
        g   = assess_ramachandran(pdb, Z, phi_bins, psi_bins, rama500_general_data)
        gly = assess_ramachandran_gly(pdb, rama500_gly_data)
        pro = assess_ramachandran_pro(pdb, rama500_pro_data)
        return g, gly, pro

    stats_dn_g,   stats_dn_gly,   stats_dn_pro   = assess(outputs['de_novo_pdb'])
    stats_rest_g, stats_rest_gly, stats_rest_pro = assess(out_pdb)
    stats_al_g,   stats_al_gly,   stats_al_pro   = assess(outputs['aligned_pdb'])

    print_rama_summary('De novo — General',        stats_dn_g)
    print_rama_summary('TDC Restrained — General', stats_rest_g)
    if stats_al_g:
        print_rama_summary('Aligned — General',    stats_al_g)

    print_rama_summary('De novo — GLY',            stats_dn_gly)
    print_rama_summary('TDC Restrained — GLY',     stats_rest_gly)
    if stats_al_gly:
        print_rama_summary('Aligned — GLY',        stats_al_gly)

    print_rama_summary('De novo — PRO',            stats_dn_pro)
    print_rama_summary('TDC Restrained — PRO',     stats_rest_pro)
    if stats_al_pro:
        print_rama_summary('Aligned — PRO',        stats_al_pro)

    def red(s_dn, s_rest):
        if s_dn and s_rest:
            delta = s_dn['pct_out'] - s_rest['pct_out']
            return (f"{s_dn['pct_out']:.1f}% → {s_rest['pct_out']:.1f}%  "
                    f"({delta:+.1f}%)")
        return 'n/a'

    print(f"\n  Outlier reduction (General): {red(stats_dn_g,   stats_rest_g)}")
    print(f"  Outlier reduction (GLY)    : {red(stats_dn_gly, stats_rest_gly)}")
    print(f"  Outlier reduction (PRO)    : {red(stats_dn_pro, stats_rest_pro)}")

    quality_dat = os.path.join(final_dir, OUT_QUALITY_DAT)
    write_quality_dat(
        quality_dat,
        stats_dn_g,   stats_rest_g,   stats_al_g,
        stats_dn_gly, stats_rest_gly, stats_al_gly,
        stats_pro_dn  = stats_dn_pro,
        stats_pro_rest= stats_rest_pro,
        stats_pro_al  = stats_al_pro,
        rmsd_data     = rmsd_data,
    )
    print(f"  Structure quality → {os.path.basename(quality_dat)}")
    timer.tick('Ramachandran quality assessment')

    # ── Plots ──────────────────────────────────────────────────────────────────
    if not args.no_plots:
        rama_png = os.path.join(final_dir, OUT_RAMA_PNG)
        plot_ramachandran(
            de_novo_pdb          = outputs['de_novo_pdb'],
            restrained_pdb       = out_pdb,
            aligned_pdb          = outputs['aligned_pdb'],
            output_png           = rama_png,
            Z                    = Z,
            phi_bins             = phi_bins,
            psi_bins             = psi_bins,
            rama500_general_data = rama500_general_data,
        )
        rama_gp_png = os.path.join(final_dir, OUT_RAMA_GLYRO_PNG)
        plot_ramachandran_gly_pro(
            de_novo_pdb      = outputs['de_novo_pdb'],
            restrained_pdb   = out_pdb,
            aligned_pdb      = outputs['aligned_pdb'],
            output_png       = rama_gp_png,
            rama500_gly_data = rama500_gly_data,
            rama500_pro_data = rama500_pro_data,
        )
    timer.tick('Generate Ramachandran plots')

    # ── Timing report ──────────────────────────────────────────────────────────
    print()
    print(timer.report())

    # ── Final summary ──────────────────────────────────────────────────────────
    W   = 100
    sep = '-' * W
    print(sep)
    print(f"{'Refinement complete!':^{W}}")
    print(sep)
    print()
    print(f"  Final restrained structure : {os.path.basename(out_pdb)}")
    print(f"  TDC B-factor structure     : {os.path.basename(tdc_pdb)}")
    print(f"  Ramachandran quality       : {os.path.basename(quality_dat)}")
    if not args.no_plots:
        print(f"  General Ramachandran plot  : {os.path.basename(rama_png)}")
        print(f"  GLY/PRO Ramachandran plot  : {os.path.basename(rama_gp_png)}")
    print(f"  GROMACS logs               : {GROMACS_SUBDIR}/")
    print(f"  Full run log               : {os.path.basename(log_path)}")
    print()
    print(f"  TDC visualisation in PyMOL:")
    print(f"    spectrum b, blue_white_red, minimum=0, maximum=100")
    print()
    print(f"  TDC visualisation in VMD:")
    print(f"    mol modcolor 0 0 Beta  |  color scale method BWR")
    print()

    sys.stdout = tee._terminal
    tee.close()


if __name__ == '__main__':
    main()
