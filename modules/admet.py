# =============================================================================
# modules/admet.py
# =============================================================================
# ADMET = Absorption, Distribution, Metabolism, Excretion, Toxicity
#
# This module takes the docking hits and applies CNS drug filters.
# All calculations are done using RDKit — completely free, no API needed.
#
# Why ADMET matters:
#   A compound can bind beautifully to a target but still fail as a drug if:
#     - It can't be absorbed orally (Absorption)
#     - It can't cross the blood-brain barrier (Distribution) ← critical for CNS
#     - It gets broken down too fast by liver enzymes (Metabolism)
#     - It damages kidneys or accumulates toxically (Excretion/Toxicity)
#
# What this module calculates:
#   1. Lipinski Rule of Five     — basic oral drug-likeness
#   2. logBB estimate            — blood-brain barrier permeability
#   3. TPSA                      — topological polar surface area (CNS filter)
#   4. Rotatable bonds           — molecular flexibility
#   5. Synthetic Accessibility   — how hard is it to make this molecule?
#   6. PAINS alerts              — known problematic/toxic substructures
# =============================================================================

import logging
from typing import List, Dict, Tuple

log = logging.getLogger(__name__)


# ── CNS Drug Thresholds (from your paper Section 3.7) ────────────────────────
# These are well-established rules for CNS-targeted compounds
THRESHOLDS = {
    "mw"              : (0, 450),      # Molecular weight (Da) — CNS stricter than Ro5's 500
    "logp"            : (-0.5, 5.0),   # Lipophilicity — CNS needs logP 1–5 ideally
    "hbd"             : (0, 3),        # H-bond donors — CNS stricter than Ro5's 5
    "hba"             : (0, 7),        # H-bond acceptors
    "tpsa"            : (0, 90),       # TPSA < 90 Å² for CNS (vs 140 for periphery)
    "rot_bonds"       : (0, 8),        # Rotatable bonds — CNS prefers rigid scaffolds
    "logbb_min"       : -1.0,          # logBB > -1.0 means good BBB penetration
    "sa_score_max"    : 6.0,           # Synthetic accessibility (1=easy, 10=hard)
    "lipinski_max_viol": 1,            # Max Lipinski violations allowed
}


def calculate_properties(smiles: str) -> Dict:
    """
    Calculate all ADMET-relevant molecular properties from a SMILES string.

    Args:
        smiles : SMILES string for the compound

    Returns:
        Dict of calculated properties, or empty dict if SMILES is invalid
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
        from rdkit.Chem.FilterCatalog import FilterCatalogParams
    except ImportError:
        raise ImportError("RDKit required: pip install rdkit")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.error(f"  [admet] Invalid SMILES: {smiles}")
        return {}

    # ── Basic Lipinski properties ─────────────────────────────────────────────
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)     # H-bond donors
    hba  = rdMolDescriptors.CalcNumHBA(mol)     # H-bond acceptors
    tpsa = Descriptors.TPSA(mol)
    rot  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings= rdMolDescriptors.CalcNumRings(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    hac  = mol.GetNumHeavyAtoms()               # Heavy atom count

    # ── Lipinski violations ───────────────────────────────────────────────────
    # Standard Ro5: MW≤500, logP≤5, HBD≤5, HBA≤10
    lipinski_violations = sum([
        mw   > 500,
        logp > 5,
        hbd  > 5,
        hba  > 10,
    ])

    # ── logBB estimation (Clark 1999 model) ───────────────────────────────────
    # logBB ≈ -0.0148 × TPSA + 0.152 × logP + 0.139
    # This is a widely used empirical model for BBB permeability estimation.
    # logBB > 0    = readily crosses BBB
    # logBB -1 to 0 = moderate BBB penetration
    # logBB < -1   = poor BBB penetration (CNS drugs should be > -1.0)
    logbb = -0.0148 * tpsa + 0.152 * logp + 0.139

    # ── Synthetic Accessibility Score (Ertl & Schuffenhauer 2009) ────────────
    # Scores from 1 (trivially easy) to 10 (nearly impossible to synthesize)
    # Drug-like compounds typically score 1–5
    try:
        from rdkit.Chem import RDConfig
        import os, sys
        sa_path = os.path.join(RDConfig.RDContribDir, 'SA_Score')
        sys.path.append(sa_path)
        import sascorer
        sa_score = sascorer.calculateScore(mol)
    except Exception:
        # SA score contrib not always available — estimate from heavy atom count
        # Simple fallback: larger/more complex molecules score higher
        sa_score = min(1.0 + (hac / 10.0) + (rings * 0.3), 9.0)
        log.debug("  [admet] SA Score contrib not found — using estimate")

    # ── PAINS (Pan-Assay Interference Compounds) alerts ──────────────────────
    # PAINS are substructures known to cause false positives in drug screens.
    # Examples: quinones, catechols, Michael acceptors.
    # Any PAINS alert is a red flag.
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        pains_alerts = catalog.GetMatches(mol)
        pains_count = len(pains_alerts)
        pains_names = [m.GetDescription() for m in pains_alerts]
    except Exception:
        pains_count = 0
        pains_names = []

    return {
        "mw"                  : round(mw, 2),
        "logp"                : round(logp, 2),
        "hbd"                 : hbd,
        "hba"                 : hba,
        "tpsa"                : round(tpsa, 2),
        "rot_bonds"           : rot,
        "rings"               : rings,
        "arom_rings"          : arom,
        "heavy_atoms"         : hac,
        "logbb"               : round(logbb, 3),
        "sa_score"            : round(sa_score, 2),
        "lipinski_violations" : lipinski_violations,
        "pains_count"         : pains_count,
        "pains_names"         : pains_names,
    }


def apply_cns_filters(props: Dict, name: str) -> Tuple[str, List[str]]:
    """
    Apply CNS drug filters and return PASS/FAIL with reasons.

    Args:
        props : property dict from calculate_properties()
        name  : compound name (for logging)

    Returns:
        Tuple of (decision, list_of_failure_reasons)
        decision is "PASS" or "FAIL"
    """
    if not props:
        return "FAIL", ["Could not calculate properties — check SMILES"]

    failures = []
    t = THRESHOLDS

    # Check each property against threshold
    if not (t["mw"][0] <= props["mw"] <= t["mw"][1]):
        failures.append(f"MW={props['mw']} (need ≤{t['mw'][1]} Da for CNS)")

    if not (t["logp"][0] <= props["logp"] <= t["logp"][1]):
        failures.append(f"logP={props['logp']} (need {t['logp'][0]}–{t['logp'][1]})")

    if props["hbd"] > t["hbd"][1]:
        failures.append(f"HBD={props['hbd']} (need ≤{t['hbd'][1]} for CNS)")

    if props["hba"] > t["hba"][1]:
        failures.append(f"HBA={props['hba']} (need ≤{t['hba'][1]})")

    if props["tpsa"] > t["tpsa"][1]:
        failures.append(f"TPSA={props['tpsa']} Å² (need <{t['tpsa'][1]} for CNS)")

    if props["rot_bonds"] > t["rot_bonds"][1]:
        failures.append(f"RotBonds={props['rot_bonds']} (need ≤{t['rot_bonds'][1]})")

    if props["logbb"] < t["logbb_min"]:
        failures.append(f"logBB={props['logbb']} (need >{t['logbb_min']} for CNS)")

    if props["sa_score"] > t["sa_score_max"]:
        failures.append(f"SA Score={props['sa_score']} (need ≤{t['sa_score_max']})")

    if props["lipinski_violations"] > t["lipinski_max_viol"]:
        failures.append(f"Lipinski violations={props['lipinski_violations']} (need ≤{t['lipinski_max_viol']})")

    if props["pains_count"] > 0:
        failures.append(f"PAINS alerts={props['pains_count']}: {', '.join(props['pains_names'])}")

    decision = "PASS" if not failures else "FAIL"

    if decision == "PASS":
        log.info(f"  [admet] ✓ {name} — PASS (all CNS filters met)")
    else:
        log.info(f"  [admet] ✗ {name} — FAIL ({len(failures)} issue(s))")
        for f in failures:
            log.info(f"  [admet]     • {f}")

    return decision, failures


def run_admet(scored_results: List[Dict], ligands: List[Dict]) -> List[Dict]:
    """
    Run full ADMET analysis on all docking hits.

    This is the main function called by run_pipeline.py.
    It processes ALL compounds (not just hits) so you can see
    why non-hits also failed ADMET.

    Args:
        scored_results : ranked results list from scorer.score_results()
        ligands        : original ligand list from config/settings.py
                         (needed to look up SMILES for each compound)

    Returns:
        scored_results list with ADMET properties and decision added to each row
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  ADMET & BBB ANALYSIS")
    log.info(f"{'═'*50}")

    # Build a name → SMILES lookup from the ligands list
    smiles_lookup = {lig["name"]: lig["smiles"] for lig in ligands}

    # Cache both properties AND decision so each compound is only
    # calculated and logged ONCE, even if it appears in multiple rows
    # (same ligand docked against multiple targets)
    cache = {}   # name -> {"props": ..., "decision": ..., "failures": ...}

    for row in scored_results:
        name   = row["ligand"]
        smiles = smiles_lookup.get(name)

        if smiles is None:
            log.warning(f"  [admet] No SMILES found for {name} — skipping ADMET")
            row["admet_decision"] = "UNKNOWN"
            row["admet_failures"] = ["SMILES not found"]
            continue

        # First time seeing this compound — calculate and log
        if name not in cache:
            log.info(f"\n  [admet] Analysing: {name}")
            props    = calculate_properties(smiles)
            decision, failures = apply_cns_filters(props, name)
            cache[name] = {
                "props"   : props,
                "decision": decision,
                "failures": failures,
            }

        # Reuse cached result silently — no duplicate logging
        props    = cache[name]["props"]
        decision = cache[name]["decision"]
        failures = cache[name]["failures"]

        # Add everything to the result row
        row.update(props)
        row["admet_decision"] = decision
        row["admet_failures"] = "; ".join(failures) if failures else "None"

        # Combined final decision: must pass BOTH docking AND ADMET
        row["final_decision"] = (
            "✓ LEAD CANDIDATE"
            if row["decision"] == "PASS" and decision == "PASS"
            else "✗ FAIL"
        )

    # Summary
    leads = [r for r in scored_results if r.get("final_decision") == "✓ LEAD CANDIDATE"]
    log.info(f"\n  ADMET complete — {len(leads)} lead candidate(s) passed both docking + ADMET")

    return scored_results