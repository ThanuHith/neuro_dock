# =============================================================================
# modules/interaction.py
# =============================================================================
# Protein-Ligand Interaction Analysis
#
# Analyses the docked binding poses to identify:
#   1. Hydrogen bonds         — donor/acceptor pairs within 3.5 Å
#   2. Hydrophobic contacts   — C-C contacts within 4.0 Å
#   3. Pi-stacking            — aromatic ring pairs within 5.5 Å
#   4. Salt bridges           — charged group pairs within 4.0 Å
#   5. Key residue contacts   — any atom within 4.5 Å of ligand
#
# Uses only RDKit + standard Python — no extra installs needed.
#
# Input:  docked .pdbqt files from data/results/<PDB_ID>/<ligand>_docked.pdbqt
# Output: data/results/interactions/
#           interactions_summary.csv   — all contacts, one row per interaction
#           interactions_report.txt    — human-readable summary
# =============================================================================

import logging
import math
import csv
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict

log = logging.getLogger(__name__)


# ── Interaction distance thresholds (Å) ──────────────────────────────────────
THRESHOLDS = {
    "hbond"       : 3.5,    # Hydrogen bond (donor-acceptor distance)
    "hydrophobic" : 4.0,    # Hydrophobic C-C contact
    "pistacking"  : 5.5,    # Pi-stacking (ring centroid distance)
    "saltbridge"  : 4.0,    # Salt bridge (charged groups)
    "contact"     : 4.5,    # General close contact (any atoms)
}

# ── Atom type definitions ─────────────────────────────────────────────────────
HBOND_DONORS    = {"N", "O", "S"}       # Atoms that can donate H-bonds
HBOND_ACCEPTORS = {"N", "O", "S", "F"} # Atoms that can accept H-bonds
HYDROPHOBIC     = {"C"}                 # Carbon atoms for hydrophobic contacts
CHARGED_POS     = {"N"}                 # Positively charged (protonated N)
CHARGED_NEG     = {"O"}                 # Negatively charged (carboxylate O)

# ── Amino acid hydrophobicity classification ──────────────────────────────────
HYDROPHOBIC_RESIDUES = {
    "ALA", "VAL", "ILE", "LEU", "MET",
    "PHE", "TRP", "PRO", "TYR"
}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
CHARGED_POS_RESIDUES = {"ARG", "LYS", "HIS"}
CHARGED_NEG_RESIDUES = {"ASP", "GLU"}
HBOND_RESIDUES = {
    "SER", "THR", "CYS", "TYR", "ASN", "GLN",
    "HIS", "ARG", "LYS", "ASP", "GLU", "TRP"
}


# =============================================================================
# PDB/PDBQT Parsing
# =============================================================================

def parse_pdbqt(pdbqt_path: Path) -> List[Dict]:
    """
    Parse a PDBQT file and return a list of atom dicts.

    Each atom dict has:
        name, resname, chain, resnum, x, y, z, element, record
    """
    atoms = []
    with open(pdbqt_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            try:
                atom = {
                    "record" : record,
                    "name"   : line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain"  : line[21].strip(),
                    "resnum" : int(line[22:26].strip()),
                    "x"      : float(line[30:38]),
                    "y"      : float(line[38:46]),
                    "z"      : float(line[46:54]),
                    "element": line[76:78].strip() if len(line) > 76 else line[12:14].strip()[0],
                }
                atoms.append(atom)
            except (ValueError, IndexError):
                continue
    return atoms


def parse_docked_ligand(pdbqt_path: Path) -> List[Dict]:
    """
    Parse only the best pose (MODEL 1) from a docked ligand PDBQT file.
    Vina outputs multiple poses — we take only the first (best energy).
    """
    atoms = []
    in_model1 = False

    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                model_num = line.split()[1] if len(line.split()) > 1 else "1"
                in_model1 = (model_num == "1")
                continue
            if line.startswith("ENDMDL"):
                if in_model1:
                    break   # Stop after first model
                continue

            # If file has no MODEL records, treat everything as model 1
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue

            try:
                atom = {
                    "record" : record,
                    "name"   : line[12:16].strip(),
                    "resname": "LIG",
                    "chain"  : "L",
                    "resnum" : 1,
                    "x"      : float(line[30:38]),
                    "y"      : float(line[38:46]),
                    "z"      : float(line[46:54]),
                    "element": line[76:78].strip() if len(line) > 76 else "C",
                }
                atoms.append(atom)
            except (ValueError, IndexError):
                continue

    return atoms


# =============================================================================
# Distance & Geometry Utilities
# =============================================================================

def distance(a1: Dict, a2: Dict) -> float:
    """Euclidean distance between two atom dicts."""
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )


def centroid(atoms: List[Dict]) -> Tuple[float, float, float]:
    """Calculate the centroid (centre of mass) of a list of atoms."""
    n = len(atoms)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (
        sum(a["x"] for a in atoms) / n,
        sum(a["y"] for a in atoms) / n,
        sum(a["z"] for a in atoms) / n,
    )


def centroid_distance(c1: Tuple, c2: Tuple) -> float:
    """Distance between two centroids."""
    return math.sqrt(sum((a-b)**2 for a, b in zip(c1, c2)))


# =============================================================================
# Interaction Detection
# =============================================================================

def find_contacts(
    receptor_atoms : List[Dict],
    ligand_atoms   : List[Dict],
    cutoff         : float,
) -> List[Dict]:
    """
    Find all receptor-ligand atom pairs within cutoff distance.
    Returns list of contact dicts.
    """
    contacts = []
    for r_atom in receptor_atoms:
        for l_atom in ligand_atoms:
            d = distance(r_atom, l_atom)
            if d <= cutoff:
                contacts.append({
                    "res_name"  : r_atom["resname"],
                    "res_num"   : r_atom["resnum"],
                    "res_atom"  : r_atom["name"],
                    "res_elem"  : r_atom["element"],
                    "lig_atom"  : l_atom["name"],
                    "lig_elem"  : l_atom["element"],
                    "distance"  : round(d, 3),
                })
    return contacts


def classify_interactions(contacts: List[Dict]) -> List[Dict]:
    """
    Classify each contact by interaction type based on atom elements
    and residue identity.

    Returns list of interaction dicts with added "type" field.
    """
    interactions = []

    for c in contacts:
        r_elem   = c["res_elem"].upper()
        l_elem   = c["lig_elem"].upper()
        res_name = c["res_name"].upper()
        d        = c["distance"]
        itype    = None

        # ── Hydrogen bond ─────────────────────────────────────────────────────
        # Both atoms must be potential H-bond donors or acceptors (N, O, S)
        # Distance ≤ 3.5 Å
        if (r_elem in HBOND_DONORS or r_elem in HBOND_ACCEPTORS) and \
           (l_elem in HBOND_DONORS or l_elem in HBOND_ACCEPTORS) and \
           d <= THRESHOLDS["hbond"]:
            itype = "Hydrogen Bond"

        # ── Salt bridge ───────────────────────────────────────────────────────
        # Charged residue (Arg, Lys, Asp, Glu) + charged ligand atom
        elif res_name in CHARGED_POS_RESIDUES and l_elem == "O" and \
             d <= THRESHOLDS["saltbridge"]:
            itype = "Salt Bridge"
        elif res_name in CHARGED_NEG_RESIDUES and l_elem == "N" and \
             d <= THRESHOLDS["saltbridge"]:
            itype = "Salt Bridge"

        # ── Hydrophobic contact ───────────────────────────────────────────────
        # Carbon-carbon contact between hydrophobic residue and ligand
        elif r_elem == "C" and l_elem == "C" and \
             res_name in HYDROPHOBIC_RESIDUES and \
             d <= THRESHOLDS["hydrophobic"]:
            itype = "Hydrophobic"

        # ── General close contact ─────────────────────────────────────────────
        elif d <= THRESHOLDS["contact"]:
            itype = "Close Contact"

        if itype:
            interactions.append({**c, "type": itype})

    return interactions


def find_pistacking(
    receptor_atoms : List[Dict],
    ligand_atoms   : List[Dict],
) -> List[Dict]:
    """
    Detect pi-stacking between aromatic residues and ligand aromatic rings.

    Strategy:
    - Group receptor atoms by residue, find aromatic residues
    - Find groups of aromatic-looking ligand atoms (C/N in rings)
    - Check centroid distances
    """
    pistacking = []

    # Group receptor atoms by residue
    by_residue = defaultdict(list)
    for atom in receptor_atoms:
        key = (atom["resname"], atom["resnum"])
        by_residue[key].append(atom)

    # Find aromatic residues in receptor
    aromatic_receptor = {}
    for (resname, resnum), atoms in by_residue.items():
        if resname.upper() in AROMATIC_RESIDUES:
            # Take ring atoms (non-backbone C and N atoms)
            ring_atoms = [a for a in atoms if a["element"] in ("C", "N")
                          and a["name"] not in ("CA", "CB", "C")]
            if len(ring_atoms) >= 5:
                aromatic_receptor[(resname, resnum)] = ring_atoms

    # Simple ligand aromatic detection: cluster carbon/nitrogen atoms
    # that are close together (within 2 Å = aromatic ring neighbours)
    lig_aromatic = [a for a in ligand_atoms if a["element"] in ("C", "N")]
    lig_centroid = centroid(lig_aromatic) if lig_aromatic else None

    if not lig_centroid:
        return []

    for (resname, resnum), ring_atoms in aromatic_receptor.items():
        rec_centroid = centroid(ring_atoms)
        d = centroid_distance(rec_centroid, lig_centroid)
        if d <= THRESHOLDS["pistacking"]:
            pistacking.append({
                "res_name" : resname,
                "res_num"  : resnum,
                "res_atom" : "ring",
                "res_elem" : "C",
                "lig_atom" : "ring",
                "lig_elem" : "C",
                "distance" : round(d, 3),
                "type"     : "Pi-Stacking",
            })

    return pistacking


# =============================================================================
# Main Analysis Function
# =============================================================================

def analyse_pair(
    receptor_pdbqt : Path,
    docked_pdbqt   : Path,
    ligand_name    : str,
    target_name    : str,
    target_pdb_id  : str,
) -> List[Dict]:
    """
    Full interaction analysis for one ligand-target pair.

    Args:
        receptor_pdbqt : prepared receptor PDBQT file
        docked_pdbqt   : docked ligand PDBQT file (output from Vina)
        ligand_name    : name label
        target_name    : target label
        target_pdb_id  : PDB ID of target

    Returns:
        List of interaction dicts, one per unique interaction found
    """
    log.info(f"  [interaction] Analysing: {ligand_name} vs {target_name}")

    if not receptor_pdbqt.exists():
        log.error(f"  [interaction] Receptor not found: {receptor_pdbqt}")
        return []

    if not docked_pdbqt.exists():
        log.error(f"  [interaction] Docked pose not found: {docked_pdbqt}")
        return []

    # Parse files
    receptor_atoms = parse_pdbqt(receptor_pdbqt)
    ligand_atoms   = parse_docked_ligand(docked_pdbqt)

    if not receptor_atoms:
        log.error(f"  [interaction] No atoms parsed from receptor")
        return []

    if not ligand_atoms:
        log.error(f"  [interaction] No atoms parsed from docked pose")
        return []

    log.info(f"  [interaction]   Receptor: {len(receptor_atoms)} atoms, "
             f"Ligand: {len(ligand_atoms)} atoms")

    # Find all close contacts first
    contacts = find_contacts(receptor_atoms, ligand_atoms, THRESHOLDS["contact"])

    # Classify by interaction type
    interactions = classify_interactions(contacts)

    # Add pi-stacking separately
    pistacking = find_pistacking(receptor_atoms, ligand_atoms)
    interactions.extend(pistacking)

    # Add metadata to each interaction
    for ix in interactions:
        ix["ligand"]     = ligand_name
        ix["target"]     = target_name
        ix["pdb_id"]     = target_pdb_id

    # Remove duplicate residue-type combinations
    # (keep only the closest contact per residue per interaction type)
    seen = {}
    unique = []
    for ix in interactions:
        key = (ix["res_name"], ix["res_num"], ix["type"])
        if key not in seen or ix["distance"] < seen[key]["distance"]:
            seen[key] = ix
    unique = sorted(seen.values(), key=lambda x: x["distance"])

    # Count by type
    type_counts = defaultdict(int)
    for ix in unique:
        type_counts[ix["type"]] += 1

    log.info(f"  [interaction]   Found {len(unique)} unique interactions:")
    for itype, count in sorted(type_counts.items()):
        log.info(f"  [interaction]     {itype}: {count}")

    return unique


def run_interaction_analysis(
    results         : List[Dict],
    receptor_pdbqts : Dict,
    results_dir     : Path,
    targets         : List[Dict],
    top_n           : int = 5,
) -> List[Dict]:
    """
    Run interaction analysis on the top N docking results.

    Args:
        results         : scored results from scorer (ranked by affinity)
        receptor_pdbqts : dict mapping pdb_id → receptor PDBQT path
        results_dir     : base results directory
        targets         : target list from settings.py
        top_n           : how many top hits to analyse (default 5)

    Returns:
        Full list of all interaction records
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  INTERACTION ANALYSIS (top {top_n} docking results)")
    log.info(f"{'═'*50}")

    # Build PDB ID lookup
    pdb_lookup = {t["name"]: t["pdb_id"] for t in targets}

    # Take top N results by binding energy
    top_results = sorted(
        results,
        key=lambda x: x["affinity_kcal_mol"]
    )[:top_n]

    all_interactions = []
    interaction_dir = results_dir / "interactions"
    interaction_dir.mkdir(parents=True, exist_ok=True)

    for row in top_results:
        ligand_name = row["ligand"]
        target_name = row["target"]
        pdb_id      = pdb_lookup.get(target_name, "")
        affinity    = row["affinity_kcal_mol"]

        log.info(f"\n  --- {ligand_name} vs {target_name} "
                 f"(ΔG = {affinity} kcal/mol) ---")

        # Locate files
        receptor_pdbqt = receptor_pdbqts.get(pdb_id)
        docked_pdbqt   = results_dir / pdb_id / f"{ligand_name}_docked.pdbqt"

        if receptor_pdbqt is None:
            log.warning(f"  [interaction] No receptor path for {pdb_id}")
            continue

        interactions = analyse_pair(
            receptor_pdbqt = receptor_pdbqt,
            docked_pdbqt   = docked_pdbqt,
            ligand_name    = ligand_name,
            target_name    = target_name,
            target_pdb_id  = pdb_id,
        )

        # Add affinity to each interaction row
        for ix in interactions:
            ix["affinity_kcal_mol"] = affinity

        all_interactions.extend(interactions)

    # Save to CSV
    if all_interactions:
        csv_path = interaction_dir / "interactions_summary.csv"
        cols = [
            "ligand", "target", "pdb_id", "affinity_kcal_mol",
            "type", "res_name", "res_num", "res_atom",
            "lig_atom", "distance"
        ]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(all_interactions)
        log.info(f"\n  [interaction] CSV saved: {csv_path}")

        # Save human-readable report
        _write_text_report(all_interactions, top_results, interaction_dir)

    log.info(f"\n  Interaction analysis complete — "
             f"{len(all_interactions)} total interactions found")

    return all_interactions


def _write_text_report(
    interactions : List[Dict],
    top_results  : List[Dict],
    out_dir      : Path,
) -> None:
    """Write a human-readable interaction report."""
    report_path = out_dir / "interactions_report.txt"

    # Group by ligand-target pair
    pairs = defaultdict(list)
    for ix in interactions:
        key = (ix["ligand"], ix["target"], ix["affinity_kcal_mol"])
        pairs[key].append(ix)

    lines = []
    lines.append("=" * 65)
    lines.append("  NEURO DOCK — PROTEIN-LIGAND INTERACTION REPORT")
    lines.append("=" * 65)
    lines.append("")

    for (ligand, target, affinity), ixs in sorted(pairs.items(),
                                                   key=lambda x: x[0][2]):
        lines.append(f"┌─ {ligand}  vs  {target}")
        lines.append(f"│  Binding affinity: {affinity} kcal/mol")
        lines.append(f"│  Total interactions: {len(ixs)}")
        lines.append("│")

        # Group by interaction type
        by_type = defaultdict(list)
        for ix in ixs:
            by_type[ix["type"]].append(ix)

        type_order = [
            "Hydrogen Bond", "Salt Bridge",
            "Pi-Stacking", "Hydrophobic", "Close Contact"
        ]

        for itype in type_order:
            if itype not in by_type:
                continue
            type_ixs = sorted(by_type[itype], key=lambda x: x["distance"])
            lines.append(f"│  {itype} ({len(type_ixs)}):")
            for ix in type_ixs[:8]:  # Show up to 8 per type
                lines.append(
                    f"│    {ix['res_name']}{ix['res_num']:>4} "
                    f"{ix['res_atom']:<4} ··· {ix['lig_atom']:<4}  "
                    f"{ix['distance']:.2f} Å"
                )

        lines.append("└" + "─" * 63)
        lines.append("")

    # Key residues summary
    lines.append("=" * 65)
    lines.append("  KEY BINDING RESIDUES SUMMARY")
    lines.append("=" * 65)

    hbond_residues = [
        ix for ix in interactions if ix["type"] == "Hydrogen Bond"
    ]
    if hbond_residues:
        lines.append("\n  Hydrogen Bond Partners:")
        seen = set()
        for ix in sorted(hbond_residues, key=lambda x: x["distance"]):
            key = f"{ix['res_name']}{ix['res_num']}"
            if key not in seen:
                seen.add(key)
                lines.append(
                    f"    {ix['ligand']} ↔ {key} "
                    f"({ix['distance']:.2f} Å) in {ix['target']}"
                )

    hydrophobic = [
        ix for ix in interactions if ix["type"] == "Hydrophobic"
    ]
    if hydrophobic:
        lines.append("\n  Hydrophobic Contacts:")
        seen = set()
        for ix in sorted(hydrophobic, key=lambda x: x["distance"]):
            key = f"{ix['res_name']}{ix['res_num']}"
            if key not in seen:
                seen.add(key)
                lines.append(
                    f"    {ix['ligand']} ↔ {key} "
                    f"({ix['distance']:.2f} Å) in {ix['target']}"
                )

    lines.append("\n" + "=" * 65)
    lines.append("  These residues are the primary docking targets.")
    lines.append("  Use PyMOL to visualise: open the receptor PDBQT")
    lines.append("  and docked ligand PDBQT, then select these residues.")
    lines.append("=" * 65)

    report_path.write_text("\n".join(lines))
    log.info(f"  [interaction] Text report saved: {report_path}")