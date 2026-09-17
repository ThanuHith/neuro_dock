#!/usr/bin/env python3
# =============================================================================
# run_pipeline.py — NeuroDock Main Entry Point
# =============================================================================
#
# This is the ONLY file you need to run:
#   python run_pipeline.py
#
# It orchestrates the full docking pipeline in order:
#   1. Setup   — create folders, configure logging
#   2. Proteins — download & prepare receptor PDBQT files
#   3. Ligands  — convert SMILES → 3D PDBQT files
#   4. Docking  — run AutoDock Vina for all ligand-target pairs
#   5. Scoring  — filter results by binding energy threshold
#   6. Report   — save CSV results and print summary table
#
# To customise targets or ligands, edit config/settings.py.
# =============================================================================

import logging
import sys
from pathlib import Path

# ── Make sure Python can find our modules ────────────────────────────────────
# (Adds the project root to the Python path — needed when running as a script)
sys.path.insert(0, str(Path(__file__).resolve().parent))

# ── Import our modules ────────────────────────────────────────────────────────
from config   import settings
from modules import protein_prep, ligand_prep, docking, scorer, admet, reporter, pdf_report, interaction, pubchem_fetch, figures, nitrosamine


def setup_logging() -> None:
    """
    Configure the Python logging system.

    Logs go to BOTH:
      - The terminal (INFO level and above)
      - A log file at data/results/pipeline.log (DEBUG level and above)

    This means you see concise output in the terminal but have full
    debug information saved if something goes wrong.
    """
    settings.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    log_format = "%(asctime)s  %(levelname)-8s  %(message)s"
    date_format = "%H:%M:%S"

    # Root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)   # Capture everything internally

    # Terminal handler — INFO and above
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(getattr(logging, settings.LOG_LEVEL, logging.INFO))
    console.setFormatter(logging.Formatter(log_format, datefmt=date_format))
    root_logger.addHandler(console)

    # File handler — DEBUG and above (full trace)
    file_handler = logging.FileHandler(settings.LOG_FILE, mode="w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
    root_logger.addHandler(file_handler)


def check_dependencies() -> None:
    """
    Check that required tools are installed before starting.
    Gives clear error messages if something is missing.
    """
    import subprocess

    log = logging.getLogger(__name__)
    errors = []

    # Check AutoDock Vina
    try:
        result = subprocess.run(
            [settings.VINA_EXECUTABLE, "--version"],
            capture_output=True, text=True, timeout=5
        )
        log.info(f"  AutoDock Vina: {result.stdout.strip()[:60]}")
    except FileNotFoundError:
        errors.append(
            f"AutoDock Vina not found (tried '{settings.VINA_EXECUTABLE}'). "
            f"Run: bash scripts/install_tools.sh"
        )

    # Check Open Babel
    try:
        result = subprocess.run(
            ["obabel", "--version"],
            capture_output=True, text=True, timeout=5
        )
        log.info(f"  Open Babel: found")
    except FileNotFoundError:
        errors.append(
            "Open Babel (obabel) not found. "
            "Run: bash scripts/install_tools.sh"
        )

    # Check RDKit
    try:
        import rdkit
        log.info(f"  RDKit: {rdkit.__version__}")
    except ImportError:
        errors.append("RDKit not installed. Run: pip install rdkit")

    if errors:
        log.error("\n  ✗ Missing dependencies:\n")
        for e in errors:
            log.error(f"    • {e}")
        log.error("\n  Fix the above, then re-run: python run_pipeline.py")
        sys.exit(1)

    log.info("  All dependencies found ✓\n")


def main() -> None:
    # ── Setup ────────────────────────────────────────────────────────────────
    setup_logging()
    log = logging.getLogger(__name__)

    log.info("=" * 60)
    log.info("  NeuroDock — AI-Driven Docking Pipeline")
    log.info("  MDD & Alzheimer's Disease Drug Discovery")
    log.info("=" * 60)

    # Create all required directories
    for d in [settings.PROTEIN_DIR, settings.LIGAND_DIR, settings.RESULTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # ── Build the ligand library (original 5 + optional PubChem expansion) ──
    ligand_library = list(settings.ORIGINAL_LIGANDS)

    if getattr(settings, "ENABLE_PUBCHEM_EXPANSION", False):
        log.info(f"\n  Expanding compound library via PubChem...")
        expanded = pubchem_fetch.fetch_compound_library(
            names      = pubchem_fetch.NEUROPROTECTIVE_LIBRARY,
            cache_path = settings.PUBCHEM_CACHE_FILE,
        )
        # Merge, skipping duplicates by name (original 5 take priority)
        existing_names = {lig["name"].lower() for lig in ligand_library}
        added = 0
        for compound in expanded:
            if compound["name"].lower() not in existing_names:
                ligand_library.append(compound)
                existing_names.add(compound["name"].lower())
                added += 1
        log.info(f"  Library expansion: {len(settings.ORIGINAL_LIGANDS)} original "
                 f"+ {added} new = {len(ligand_library)} total compounds")

    log.info(f"\n  Config:")
    log.info(f"    Targets  : {len(settings.TARGETS)}")
    log.info(f"    Ligands  : {len(ligand_library)}")
    log.info(f"    Cutoff   : ≤ {settings.BINDING_ENERGY_CUTOFF} kcal/mol")
    log.info(f"    Exhaustiveness: {settings.VINA_EXHAUSTIVENESS}")

    # Check tools
    log.info(f"\n  Checking dependencies...")
    check_dependencies()

    # ── Step 1: Prepare receptors ─────────────────────────────────────────────
    receptor_pdbqts = {}   # Maps pdb_id → Path to its .pdbqt file

    for target in settings.TARGETS:
        try:
            pdbqt = protein_prep.prepare_receptor(target, settings.PROTEIN_DIR)
            receptor_pdbqts[target["pdb_id"]] = pdbqt
        except Exception as e:
            log.error(f"  Receptor prep failed for {target['pdb_id']}: {e}")
            log.error(f"  Skipping this target and continuing...")

    if not receptor_pdbqts:
        log.error("  No receptors prepared — cannot continue. Check your internet connection.")
        sys.exit(1)

    # ── Step 2: Prepare ligands ───────────────────────────────────────────────
    prepared_ligands = ligand_prep.prepare_all_ligands(
        ligand_library,
        settings.LIGAND_DIR,
    )

    if not prepared_ligands:
        log.error("  No ligands prepared — check your SMILES strings in config/settings.py")
        sys.exit(1)

    # ── Step 3: Run docking ───────────────────────────────────────────────────
    all_results = docking.dock_all(
        targets          = settings.TARGETS,
        prepared_ligands = prepared_ligands,
        receptor_pdbqts  = receptor_pdbqts,
        results_dir      = settings.RESULTS_DIR,
        vina_exe         = settings.VINA_EXECUTABLE,
        exhaustiveness   = settings.VINA_EXHAUSTIVENESS,
        num_modes        = settings.VINA_NUM_MODES,
    )

    # ── Step 4: Score & filter ────────────────────────────────────────────────
    scored_results = scorer.score_results(
        all_results = all_results,
        cutoff      = settings.BINDING_ENERGY_CUTOFF,
    )

    # ── Step 5: ADMET & BBB analysis ─────────────────────────────────────────
    scored_results = admet.run_admet(
        scored_results = scored_results,
        ligands        = ligand_library,
    )

    # ── Step 6: Interaction Analysis ─────────────────────────────────────────
    all_interactions = interaction.run_interaction_analysis(
        results         = scored_results,
        receptor_pdbqts = receptor_pdbqts,
        results_dir     = settings.RESULTS_DIR,
        targets         = settings.TARGETS,
        top_n           = 5,
    )

    # ── Step 6: Report ────────────────────────────────────────────────────────
    reporter.generate_report(
        results     = scored_results,
        results_dir = settings.RESULTS_DIR,
        cutoff      = settings.BINDING_ENERGY_CUTOFF,
    )

    # ── Step 7: PDF Report ───────────────────────────────────────────────────
    pdf_report.generate_pdf_report(
        results      = scored_results,
        results_dir  = settings.RESULTS_DIR,
        cutoff       = settings.BINDING_ENERGY_CUTOFF,
        targets      = settings.TARGETS,
        ligands      = ligand_library,
        interactions = all_interactions,
    )

    # ── Step 8: Publication Figures ──────────────────────────────────────────
    figures.generate_all_figures(
        results_dir = settings.RESULTS_DIR,
        cutoff      = settings.BINDING_ENERGY_CUTOFF,
    )

    # ── Nitrosamine Risk Assessment ──────────────────────────────────────────



    lead_smiles = []
    for r in scored_results:
        if r.get("final_decision") == "✓ LEAD CANDIDATE":
            smiles = r.get("smiles", "")
            if smiles:
                lead_smiles.append((r["ligand"], smiles))

    if lead_smiles:
        nitro_results = nitrosamine.run_nitrosamine_assessment(lead_smiles)
        report_text = nitrosamine.format_report(nitro_results)
        print(report_text)

        nitro_path = settings.RESULTS_DIR / "nitrosamine_risk_report.txt"
        with open(nitro_path, "w") as f:
            f.write(report_text)
        log.info(f"  Nitrosamine risk report saved: {nitro_path}")

    log.info(f"  Pipeline complete. Log saved to: {settings.LOG_FILE}")


if __name__ == "__main__":
    main()