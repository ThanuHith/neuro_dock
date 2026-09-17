# =============================================================================
# modules/ligand_prep.py
# =============================================================================
# Converts SMILES strings into 3D ligand files ready for AutoDock Vina.
#
# What is SMILES?
#   SMILES (Simplified Molecular Input Line Entry System) is a text notation
#   for molecules. Example: "CCO" = ethanol. You can get SMILES for any
#   compound from PubChem (https://pubchem.ncbi.nlm.nih.gov).
#
# What this module does, step by step:
#   1. Parses the SMILES string using RDKit
#   2. Adds missing hydrogen atoms (needed for 3D geometry)
#   3. Generates a 3D conformer using the ETKDG algorithm
#      (ETKDG = Experimental-Torsion Distance Geometry — gives realistic
#       bond angles and distances for drug-like molecules)
#   4. Minimizes the geometry using MMFF94 force field
#      (MMFF94 = Merck Molecular Force Field — industry standard for
#       small-molecule geometry optimization)
#   5. Saves as a .sdf file, then converts to .pdbqt via Open Babel
#      (Open Babel adds Gasteiger charges and rotatable bond markers
#       that Vina uses during flexible docking)
# =============================================================================

import logging
import subprocess
from pathlib import Path
from typing import List, Dict

log = logging.getLogger(__name__)


def smiles_to_pdbqt(ligand: Dict, ligand_dir: Path) -> Path | None:
    """
    Convert one ligand (from SMILES) to a Vina-ready .pdbqt file.

    Args:
        ligand     : dict with keys "name", "smiles", "note"
        ligand_dir : folder to save prepared ligand files

    Returns:
        Path to the .pdbqt file, or None if preparation failed
    """
    # Import RDKit here so errors are clear if it's not installed
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError:
        raise ImportError(
            "RDKit is not installed. Run: pip install rdkit"
        )

    ligand_dir.mkdir(parents=True, exist_ok=True)

    name   = ligand["name"]
    smiles = ligand["smiles"]

    # Output paths
    sdf_path   = ligand_dir / f"{name}.sdf"
    pdbqt_path = ligand_dir / f"{name}.pdbqt"

    if pdbqt_path.exists():
        log.info(f"  [ligand_prep] {name}.pdbqt already exists — skipping.")
        return pdbqt_path

    log.info(f"  [ligand_prep] Preparing ligand: {name}")
    log.info(f"  [ligand_prep]   SMILES: {smiles}")

    # ── Step 1: Parse SMILES ──────────────────────────────────────────────────
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.error(f"  [ligand_prep] ✗ Invalid SMILES for {name}: '{smiles}'")
        log.error(f"  [ligand_prep]   Tip: validate your SMILES at https://pubchem.ncbi.nlm.nih.gov")
        return None

    # ── Step 2: Add hydrogens ─────────────────────────────────────────────────
    # RDKit's mol objects omit implicit hydrogens by default.
    # We need explicit H atoms for proper 3D geometry.
    mol = Chem.AddHs(mol)

    # ── Quick drug-likeness check (Lipinski Rule of Five) ─────────────────────
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    log.info(f"  [ligand_prep]   MW={mw:.1f}, logP={logp:.2f}, HBD={hbd}, HBA={hba} — Lipinski violations: {violations}")

    if violations > 1:
        log.warning(f"  [ligand_prep]   ⚠ {name} has {violations} Lipinski violations — may have poor oral bioavailability")

    # ── Step 3: Generate 3D conformer (ETKDG algorithm) ──────────────────────
    # ETKDG uses experimental torsion angle preferences from the Cambridge
    # Structural Database to generate realistic 3D starting geometries.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42       # Fixed seed for reproducibility
    params.numThreads = 1        # Single thread keeps output deterministic

    result = AllChem.EmbedMolecule(mol, params)

    if result == -1:
        log.error(
            f"  [ligand_prep] ✗ 3D conformer generation failed for {name}. "
            f"The SMILES may describe a structure too complex for ETKDG. "
            f"Try simplifying the molecule or use a pre-built SDF from PubChem."
        )
        return None

    # ── Step 4: Geometry optimization with MMFF94 ────────────────────────────
    # MMFF94 minimizes bond lengths, angles, and torsions to a local energy minimum.
    # This gives a more realistic 3D structure before docking.
    ff_result = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=2000)

    if ff_result == 1:
        log.warning(f"  [ligand_prep]   ⚠ MMFF94 optimization did not fully converge for {name} — using best geometry found")
    elif ff_result == -1:
        log.warning(f"  [ligand_prep]   ⚠ MMFF94 not applicable for {name} — using unoptimized geometry")

    # ── Step 5: Save as SDF ──────────────────────────────────────────────────
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    log.info(f"  [ligand_prep]   3D SDF saved: {sdf_path.name}")

    # ── Step 6: Convert SDF → PDBQT via Open Babel ───────────────────────────
    # Open Babel adds:
    #   - Gasteiger partial charges (same method used for the receptor)
    #   - Rotatable bond definitions (tells Vina which bonds can flex during docking)
    #   - AutoDock atom type codes
    cmd = [
        "obabel",
        "-isdf", str(sdf_path),
        "-opdbqt",
        "-O", str(pdbqt_path),
        "--partialcharge", "gasteiger",
        "-h",            # keep all hydrogens
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not pdbqt_path.exists():
        log.error(
            f"  [ligand_prep] ✗ Open Babel conversion failed for {name}.\n"
            f"  STDERR: {result.stderr}"
        )
        return None

    log.info(f"  [ligand_prep]   ✓ PDBQT ready: {pdbqt_path.name}")
    return pdbqt_path


def prepare_all_ligands(ligands: List[Dict], ligand_dir: Path) -> List[Dict]:
    """
    Prepare all ligands in the LIGANDS list from settings.py.

    Args:
        ligands    : list of dicts from config/settings.py LIGANDS
        ligand_dir : folder to save prepared ligand files

    Returns:
        List of dicts with added "pdbqt_path" key for successfully prepared ligands.
        Failed ligands are skipped with a warning.
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  LIGAND PREPARATION  ({len(ligands)} compounds)")
    log.info(f"{'═'*50}")

    prepared = []
    failed   = []

    for ligand in ligands:
        pdbqt = smiles_to_pdbqt(ligand, ligand_dir)

        if pdbqt is not None:
            prepared.append({**ligand, "pdbqt_path": pdbqt})
        else:
            failed.append(ligand["name"])

    log.info(f"\n  Ligand prep complete: {len(prepared)} ready, {len(failed)} failed")
    if failed:
        log.warning(f"  Failed ligands: {', '.join(failed)}")

    return prepared