# =============================================================================
# modules/protein_prep.py
# =============================================================================
# Downloads protein structures from RCSB and prepares them for AutoDock Vina.
#
# What this module does, step by step:
#   1. Downloads the raw .pdb file from https://rcsb.org using the PDB ID
#   2. Keeps only the chain you specified (usually chain A)
#   3. Removes crystallographic water molecules (HETATM HOH lines)
#   4. Removes any co-crystallized ligand (non-protein HETATM records)
#   5. Calls Open Babel to add polar hydrogens and convert to .pdbqt format
#      (PDBQT is the format AutoDock Vina requires — it adds partial charges
#       and atom types on top of normal PDB coordinates)
#
# Why we need PDBQT:
#   AutoDock Vina scores binding using partial atomic charges (Gasteiger method).
#   PDBQT files embed those charges and the atom-type codes Vina understands.
# =============================================================================

import logging
import requests
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)


def download_pdb(pdb_id: str, out_dir: Path) -> Path:
    """
    Download a PDB file from RCSB by its 4-letter ID.

    Args:
        pdb_id  : e.g. "4AT3"
        out_dir : folder to save the downloaded file

    Returns:
        Path to the downloaded .pdb file
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = out_dir / f"{pdb_id}.pdb"

    if pdb_path.exists():
        log.info(f"  [protein_prep] {pdb_id}.pdb already exists — skipping download.")
        return pdb_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log.info(f"  [protein_prep] Downloading {pdb_id} from {url} ...")

    response = requests.get(url, timeout=30)

    if response.status_code != 200:
        raise FileNotFoundError(
            f"Could not download {pdb_id} from RCSB "
            f"(HTTP {response.status_code}). "
            f"Check that the PDB ID is correct: {url}"
        )

    pdb_path.write_text(response.text)
    log.info(f"  [protein_prep] Saved to {pdb_path}")
    return pdb_path


def clean_pdb(pdb_path: Path, chain: str = "A") -> Path:
    """
    Clean a raw PDB file:
      - Keep only the specified chain
      - Remove water molecules (HOH)
      - Remove non-protein HETATM records (ligands, ions)
      - Keep only ATOM records (protein backbone + sidechains)

    Args:
        pdb_path : path to the downloaded .pdb file
        chain    : chain ID to keep (default "A")

    Returns:
        Path to the cleaned .pdb file (saved alongside original with _clean suffix)
    """
    clean_path = pdb_path.parent / f"{pdb_path.stem}_clean.pdb"

    log.info(f"  [protein_prep] Cleaning {pdb_path.name} — keeping chain {chain} ...")

    kept_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()

            # Keep only ATOM records (protein atoms)
            # HETATM records are heteroatoms: water, ligands, cofactors — skip them
            if record != "ATOM":
                continue

            # Check the chain ID (column 22 in PDB format, 0-indexed = column 21)
            line_chain = line[21].strip()
            if line_chain != chain:
                continue

            kept_lines.append(line)

    if not kept_lines:
        raise ValueError(
            f"No ATOM records found for chain '{chain}' in {pdb_path.name}. "
            f"Try a different chain letter — check the PDB file header for available chains."
        )

    # Always end with an END record
    kept_lines.append("END\n")

    clean_path.write_text("".join(kept_lines))
    log.info(f"  [protein_prep] Cleaned file saved: {clean_path.name} ({len(kept_lines)-1} ATOM lines)")
    return clean_path


def convert_to_pdbqt(pdb_path: Path) -> Path:
    """
    Convert a cleaned .pdb protein file to .pdbqt format using Open Babel.

    Open Babel adds:
      - Polar hydrogens (hydrogens on O, N, S — needed for H-bond scoring)
      - Gasteiger partial atomic charges
      - AutoDock atom-type codes

    Args:
        pdb_path : path to the cleaned .pdb file

    Returns:
        Path to the output .pdbqt file
    """
    pdbqt_path = pdb_path.with_suffix(".pdbqt")

    if pdbqt_path.exists():
        log.info(f"  [protein_prep] {pdbqt_path.name} already exists — skipping conversion.")
        return pdbqt_path

    log.info(f"  [protein_prep] Converting {pdb_path.name} → PDBQT via Open Babel ...")

    # obabel command explanation:
    #   -ipdb        : input format is PDB
    #   -opdbqt      : output format is PDBQT
    #   -O <file>    : output file path
    #   -r           : keep only the largest fragment (removes stray atoms)
    #   --partialcharge gasteiger : add Gasteiger charges (same method Vina uses)
    cmd = [
        "obabel",
        "-ipdb", str(pdb_path),
        "-opdbqt",
        "-O", str(pdbqt_path),
        "-r",
        "--partialcharge", "gasteiger",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not pdbqt_path.exists():
        raise RuntimeError(
            f"Open Babel failed for {pdb_path.name}.\n"
            f"STDERR: {result.stderr}\n"
            f"Make sure Open Babel is installed: bash scripts/install_tools.sh"
        )

    log.info(f"  [protein_prep] PDBQT saved: {pdbqt_path.name}")
    return pdbqt_path


def prepare_receptor(target: dict, protein_dir: Path) -> Path:
    """
    Full receptor preparation pipeline for one target.
    Calls download_pdb → clean_pdb → convert_to_pdbqt in sequence.

    Args:
        target      : one entry from config/settings.py TARGETS list
        protein_dir : folder to store protein files

    Returns:
        Path to the final ready-to-dock .pdbqt file
    """
    pdb_id = target["pdb_id"]
    chain  = target.get("chain", "A")

    log.info(f"\n{'─'*50}")
    log.info(f"  Preparing receptor: {target['name']} ({pdb_id})")
    log.info(f"{'─'*50}")

    # Step 1 — Download
    pdb_raw = download_pdb(pdb_id, protein_dir)

    # Step 2 — Clean
    pdb_clean = clean_pdb(pdb_raw, chain=chain)

    # Step 3 — Convert to PDBQT
    pdbqt = convert_to_pdbqt(pdb_clean)

    log.info(f"  [protein_prep] ✓ Receptor ready: {pdbqt.name}\n")
    return pdbqt