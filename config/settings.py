# =============================================================================
# config/settings.py
# =============================================================================
# Central configuration for the NeuroDock pipeline.
# Everything you might want to change is here — no need to dig into modules.
# =============================================================================

from pathlib import Path

# ── Project root (auto-resolved — don't change this) ──────────────────────────
ROOT = Path(__file__).resolve().parent.parent

# ── Data directories ──────────────────────────────────────────────────────────
PROTEIN_DIR = ROOT / "data" / "proteins"
LIGAND_DIR  = ROOT / "data" / "ligands"
RESULTS_DIR = ROOT / "data" / "results"

# ── AutoDock Vina executable ───────────────────────────────────────────────────
# If vina is on your PATH, "vina" is fine.
# If you installed it to a custom location, put the full path here:
#   e.g.  "/home/user/tools/vina_1.2.5_linux_x86_64"
VINA_EXECUTABLE = "vina"

# ── Docking thresholds (from your paper — Section 3.7) ───────────────────────
BINDING_ENERGY_CUTOFF = -7.0   # kcal/mol  — only hits ≤ this advance
VINA_EXHAUSTIVENESS  = 16      # Higher = more thorough search (slower)
VINA_NUM_MODES       = 9       # Number of binding poses Vina returns

# ── Grid box buffer (Å added to each side of binding site) ───────────────────
GRID_BUFFER = 6.0              # From your paper: "≥ 6 Å buffer"

# ── Targets — add or remove entries here to change what gets screened ─────────
# Each entry needs:
#   pdb_id   : 4-letter RCSB code
#   name     : human-readable label
#   chain    : which chain to keep (usually "A")
#   center   : (x, y, z) of the binding site in Å
#              → find this in PyMOL: select ligand, run "get_position"
#   box_size : (x, y, z) dimensions of the search box in Å
#              → typically 20–30 Å covers most binding pockets
#   note     : short description (for the report)
TARGETS = [
    {
        "pdb_id"   : "2V61",
        "name"     : "MAO-B",
        "chain"    : "A",
        "center"   : (55.6, 163.3, 42.1),  # FAD-binding hydrophobic pocket
        "box_size" : (20.0, 20.0, 20.0),
        "note"     : "Oxidative neuroprotection — dopamine catabolism",
    },
]

# ── Other targets (disabled for expanded MAO-B screen) ───────────────────────
# Uncomment any of these and add back to TARGETS above to re-include them.
#
# TrkB receptor (NTRK2):
#   {"pdb_id": "4AT3", "name": "TrkB receptor (NTRK2)", "chain": "A",
#    "center": (61.0, 21.0, -10.2), "box_size": (140.0, 56.5, 49.3),
#    "note": "BDNF-dependent neurotrophic signaling — MDD & AD target"},
#
# Acetylcholinesterase (AChE):
#   {"pdb_id": "4EY7", "name": "Acetylcholinesterase (AChE)", "chain": "A",
#    "center": (2.3, -13.7, 51.2), "box_size": (24.0, 24.0, 24.0),
#    "note": "Dual-site inhibition: cholinergic + anti-amyloid (PAS)"},
#
# GSK-3β:
#   {"pdb_id": "1Q5K", "name": "GSK-3β", "chain": "A",
#    "center": (83.0, 23.4, 41.3), "box_size": (184.0, 68.7, 82.6),
#    "note": "Tau hyperphosphorylation driver — AD pathology"},

# ── Original 5 ligands (kept as the founding compound set) ──────────────────
ORIGINAL_LIGANDS = [
    {
        "name"  : "NP-alpha1",
        "smiles": "Cc1ccc(OC)c(CN2CCNCC2)c1",
        "note"  : "Top hit from paper — TrkB binder, Vina −8.5 kcal/mol",
    },
    {
        "name"  : "Donepezil",
        "smiles": "COc1cc2c(cc1OC)C(=O)C(CC2)CC1CCN(CC1)Cc1ccccc1",
        "note"  : "Approved AChE inhibitor — positive control",
    },
    {
        "name"  : "Selegiline",
        "smiles": "C#CCN(C)Cc1ccccc1",
        "note"  : "Approved MAO-B inhibitor — positive control",
    },
    {
        "name"  : "Quercetin",
        "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "note"  : "Natural neuroprotective flavonoid — reference compound",
    },
    {
        "name"  : "7,8-DHF",
        "smiles": "O=c1cc(-c2ccccc2)oc2cc(O)c(O)cc12",
        "note"  : "TrkB agonist proof-of-concept from paper (Section 5.1)",
    },
]

# ── Expanded compound library (auto-fetched from PubChem) ────────────────────
# See modules/pubchem_fetch.py for the NEUROPROTECTIVE_LIBRARY name list
# and the fetch logic. Results are cached to data/pubchem_cache.json so
# PubChem isn't re-queried on every pipeline run.
#
# LIGANDS is built dynamically in run_pipeline.py by merging
# ORIGINAL_LIGANDS with the PubChem-fetched expanded library, with
# duplicate names removed (original entries take priority).
LIGANDS = ORIGINAL_LIGANDS   # placeholder — overwritten at runtime if expansion is enabled

# ── Expansion toggle ──────────────────────────────────────────────────────────
# Set to True to fetch the 25-compound expanded library from PubChem
# and merge it with ORIGINAL_LIGANDS. Set to False to use only the
# original 5 hardcoded compounds (faster, no internet required).
ENABLE_PUBCHEM_EXPANSION = True
PUBCHEM_CACHE_FILE = ROOT / "data" / "pubchem_cache.json"

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_LEVEL = "INFO"    # Options: "DEBUG", "INFO", "WARNING", "ERROR"
LOG_FILE  = ROOT / "data" / "results" / "pipeline.log"