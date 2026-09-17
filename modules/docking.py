# =============================================================================
# modules/docking.py
# =============================================================================
# Runs AutoDock Vina for every ligand–target combination and parses results.
#
# How AutoDock Vina works (simplified):
#   - Vina treats the ligand as flexible (it can rotate around single bonds)
#   - The receptor is treated as rigid (fixed protein structure)
#   - Vina searches the grid box for the orientation and conformation of the
#     ligand that gives the lowest binding energy (most negative kcal/mol)
#   - Lower (more negative) = better binding
#   - Threshold used in this pipeline: ≤ −7.0 kcal/mol (from Section 3.7)
#
# Output for each run:
#   - A .pdbqt file with the docked poses (3D coordinates of each binding mode)
#   - A .log file with the binding energy table
#   - Parsed energies stored as a list of dicts for the scorer module
# =============================================================================

import logging
import subprocess
import re
from pathlib import Path
from typing import List, Dict, Optional

log = logging.getLogger(__name__)


def build_vina_config(
    receptor_pdbqt : Path,
    center         : tuple,
    box_size       : tuple,
    exhaustiveness : int,
    num_modes      : int,
) -> dict:
    """
    Build the parameter dictionary for an AutoDock Vina run.

    Args:
        receptor_pdbqt : path to the prepared receptor .pdbqt file
        center         : (x, y, z) coordinates of the binding site center
        box_size       : (x, y, z) dimensions of the search grid box in Å
        exhaustiveness : search thoroughness (16 = paper default)
        num_modes      : number of binding poses to return

    Returns:
        Dict of Vina command-line arguments
    """
    cx, cy, cz = center
    sx, sy, sz = box_size

    return {
        "--receptor"      : str(receptor_pdbqt),
        "--center_x"      : cx,
        "--center_y"      : cy,
        "--center_z"      : cz,
        "--size_x"        : sx,
        "--size_y"        : sy,
        "--size_z"        : sz,
        "--exhaustiveness": exhaustiveness,
        "--num_modes"     : num_modes,
    }


def run_vina(
    vina_exe       : str,
    ligand_pdbqt   : Path,
    config         : dict,
    out_pdbqt      : Path,
    log_path       : Path,
) -> Optional[str]:
    """
    Execute a single AutoDock Vina docking run.

    Args:
        vina_exe     : path or name of the Vina executable (e.g. "vina")
        ligand_pdbqt : prepared ligand .pdbqt file
        config       : dict from build_vina_config()
        out_pdbqt    : where to save the output poses
        log_path     : where to save the Vina log (contains energy table)

    Returns:
        Contents of the Vina log as a string, or None on failure
    """
    # Build the command list
    # e.g.: vina --receptor rec.pdbqt --ligand lig.pdbqt --center_x 20.5 ...
    cmd = [vina_exe, "--ligand", str(ligand_pdbqt), "--out", str(out_pdbqt)]
    for flag, value in config.items():
        cmd += [str(flag), str(value)]

    log.info(f"    [docking] Running Vina: {ligand_pdbqt.stem} → {out_pdbqt.parent.name}")
    log.debug(f"    [docking] Full command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,   # 10-minute timeout per ligand-target pair
        )
    except subprocess.TimeoutExpired:
        log.error(f"    [docking] ✗ Vina timed out for {ligand_pdbqt.stem} — try reducing exhaustiveness")
        return None
    except FileNotFoundError:
        raise FileNotFoundError(
            f"AutoDock Vina not found at '{vina_exe}'. "
            f"Run: bash scripts/install_tools.sh"
        )

    # Vina writes its output to stdout — save it as the log
    vina_log = result.stdout + result.stderr

    if log_path:
        log_path.write_text(vina_log)

    if result.returncode != 0:
        log.error(f"    [docking] ✗ Vina returned error code {result.returncode}")
        log.error(f"    [docking]   Output: {vina_log[:500]}")
        return None

    if not out_pdbqt.exists():
        log.error(f"    [docking] ✗ No output PDBQT generated — Vina may have failed silently")
        return None

    log.info(f"    [docking]   ✓ Docking complete, output: {out_pdbqt.name}")
    return vina_log


def parse_vina_log(vina_log: str, ligand_name: str, target_name: str) -> List[Dict]:
    """
    Parse the binding energy table from a Vina log.

    AutoDock Vina log format (the section we care about):
        -----+------------+----------+----------
        mode |   affinity  | dist from best mode
             | (kcal/mol)  | rmsd l.b.| rmsd u.b.
        -----+------------+----------+----------
           1         -8.5      0.000      0.000
           2         -7.9      1.432      3.201
           ...

    Args:
        vina_log    : raw string output from Vina
        ligand_name : name label for this ligand
        target_name : name label for this target

    Returns:
        List of dicts, one per binding mode, with keys:
          ligand, target, mode, affinity_kcal_mol, rmsd_lb, rmsd_ub
    """
    results = []

    # Regex that matches each result row: mode number, affinity, two RMSD values
    # Example line: "   1         -8.5      0.000      0.000"
    pattern = re.compile(
        r"^\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)",
        re.MULTILINE
    )

    matches = pattern.findall(vina_log)

    if not matches:
        log.warning(
            f"    [docking] ⚠ Could not parse binding energies for "
            f"{ligand_name} vs {target_name}. "
            f"The log may be malformed — check data/results/*.log"
        )
        return results

    for match in matches:
        mode, affinity, rmsd_lb, rmsd_ub = match
        results.append({
            "ligand"           : ligand_name,
            "target"           : target_name,
            "mode"             : int(mode),
            "affinity_kcal_mol": float(affinity),
            "rmsd_lb"          : float(rmsd_lb),
            "rmsd_ub"          : float(rmsd_ub),
        })

    # The best (lowest) affinity is always mode 1
    best = results[0]["affinity_kcal_mol"]
    log.info(f"    [docking]   Best binding energy: {best} kcal/mol ({len(results)} modes)")

    return results


def dock_all(
    targets         : List[Dict],
    prepared_ligands: List[Dict],
    receptor_pdbqts : Dict[str, Path],
    results_dir     : Path,
    vina_exe        : str,
    exhaustiveness  : int,
    num_modes       : int,
) -> List[Dict]:
    """
    Run docking for ALL ligand–target combinations.

    This is the main function called by run_pipeline.py.
    It loops over every (target, ligand) pair and collects all results.

    Args:
        targets          : list of target dicts from settings.py
        prepared_ligands : list of ligand dicts with "pdbqt_path" key added
        receptor_pdbqts  : dict mapping pdb_id → Path of receptor .pdbqt
        results_dir      : folder to store docking output files
        vina_exe         : name/path of AutoDock Vina executable
        exhaustiveness   : Vina exhaustiveness parameter
        num_modes        : number of binding modes to generate

    Returns:
        Flat list of result dicts (one per mode, per ligand-target pair)
    """
    all_results = []
    total_runs  = len(targets) * len(prepared_ligands)
    run_n       = 0

    log.info(f"\n{'═'*50}")
    log.info(f"  DOCKING  ({len(prepared_ligands)} ligands × {len(targets)} targets = {total_runs} runs)")
    log.info(f"{'═'*50}")

    for target in targets:
        pdb_id      = target["pdb_id"]
        target_name = target["name"]
        receptor    = receptor_pdbqts.get(pdb_id)

        if receptor is None or not receptor.exists():
            log.error(f"  [docking] Receptor PDBQT missing for {pdb_id} — skipping target")
            continue

        # Make a subfolder per target for clean organisation
        target_dir = results_dir / pdb_id
        target_dir.mkdir(parents=True, exist_ok=True)

        config = build_vina_config(
            receptor_pdbqt = receptor,
            center         = target["center"],
            box_size       = target["box_size"],
            exhaustiveness = exhaustiveness,
            num_modes      = num_modes,
        )

        for ligand in prepared_ligands:
            run_n += 1
            ligand_name  = ligand["name"]
            ligand_pdbqt = ligand["pdbqt_path"]

            log.info(f"\n  Run {run_n}/{total_runs}: {ligand_name} vs {target_name}")

            out_pdbqt = target_dir / f"{ligand_name}_docked.pdbqt"
            log_path  = target_dir / f"{ligand_name}_vina.log"

            # Skip if already done (allows resuming interrupted runs)
            if out_pdbqt.exists() and log_path.exists():
                log.info(f"    [docking] Already done — loading existing log")
                vina_log = log_path.read_text()
            else:
                vina_log = run_vina(
                    vina_exe     = vina_exe,
                    ligand_pdbqt = ligand_pdbqt,
                    config       = config,
                    out_pdbqt    = out_pdbqt,
                    log_path     = log_path,
                )

            if vina_log is None:
                continue

            # Parse the energy table from the log
            run_results = parse_vina_log(vina_log, ligand_name, target_name)
            all_results.extend(run_results)

    log.info(f"\n  Docking complete — {len(all_results)} total pose records collected")
    return all_results