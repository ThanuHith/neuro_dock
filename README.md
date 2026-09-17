# NeuroDock — Computational Drug Discovery Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform: Linux/WSL](https://img.shields.io/badge/platform-Linux%20%7C%20WSL-lightgrey.svg)]()

**NeuroDock** is a fully open-source, end-to-end virtual screening pipeline for
computational drug discovery targeting neurodegeneration. One command takes you
from SMILES strings to publication-quality results.

## Pipeline

```
PubChem API -> 3D Ligand Prep -> AutoDock Vina Docking -> ADMET Filtering
           -> Interaction Analysis -> Nitrosamine Risk -> Figures + Report
```

## Key Results (MAO-B / PDB: 2V61)

| Compound      | dG (kcal/mol) | logBB  | TPSA  | Nitrosamine Risk |
|---------------|---------------|--------|-------|-----------------|
| Donepezil     | -8.84         | +0.287 | 38.77 | MODERATE-HIGH   |
| Chrysin       | -8.48         | -0.470 | 70.67 | LOW             |
| Resveratrol   | -8.35         | -0.307 | 60.69 | LOW             |
| Safinamide*   | -8.13         | -0.453 | 64.35 | MODERATE        |
| Galantamine*  | -7.74         | -0.200 | 41.93 | MODERATE        |
| Rivastigmine* | -7.12         | +0.073 | 32.78 | HIGH            |

*Clinically approved CNS drugs recovered as internal pipeline validation.

## Features

- Automated compound retrieval via PubChem PUG REST API with CID verification
- 3D conformer generation with RDKit ETKDGv3 + MMFF94 (fixed seed, reproducible)
- Molecular docking with AutoDock Vina 1.2.5
- CNS ADMET filtering: MW, logP, HBD/HBA, TPSA, logBB, Lipinski, PAINS
- Protein-ligand interaction analysis: H-bonds, hydrophobic, salt bridges
- Nitrosamine impurity risk assessment: ICH M7(R2) / EMA-FDA 2023 framework
- 5 publication-quality figures at 300 DPI
- Auto-generated PDF and CSV reports

## Quick Start

```bash
# 1. Clone
git clone https://github.com/YOUR_USERNAME/neuro_dock.git
cd neuro_dock

# 2. Create virtual environment
python3 -m venv venv
source venv/bin/activate

# 3. Install Python packages
pip install -r requirements.txt

# 4. Install AutoDock Vina + Open Babel
bash scripts/install_tools.sh

# 5. Run the pipeline
python3 run_pipeline.py
```

Results appear in `data/results/`.

## Project Structure

```
neuro_dock/
+-- run_pipeline.py          # Main orchestrator - run this
+-- config/
|   +-- settings.py          # Targets, ligands, thresholds
+-- modules/
|   +-- protein_prep.py      # PDB download + PDBQT conversion
|   +-- ligand_prep.py       # SMILES to 3D to PDBQT
|   +-- docking.py           # AutoDock Vina wrapper
|   +-- scorer.py            # Binding energy filtering
|   +-- admet.py             # CNS ADMET + PAINS screening
|   +-- interaction.py       # Protein-ligand interactions
|   +-- nitrosamine.py       # ICH M7 nitrosamine risk assessment
|   +-- figures.py           # Publication figure generation
|   +-- reporter.py          # CSV output
|   +-- pdf_report.py        # PDF report generation
|   +-- pubchem_fetch.py     # PubChem API retrieval
+-- data/
|   +-- proteins/            # PDB + PDBQT receptor files
|   +-- ligands/             # Prepared ligand PDBQT files
|   +-- results/             # All outputs (CSV, PDF, figures, reports)
+-- scripts/
|   +-- install_tools.sh     # System tool installer
+-- requirements.txt
+-- LICENSE
+-- README.md
```

## Customising

Edit `config/settings.py` to change:
- Target protein (PDB ID, grid box coordinates)
- Screening library (add/remove compounds)
- Docking threshold (default: -7.0 kcal/mol)
- ADMET thresholds

## Citation

```
Thanuhith J. (2025). In Silico Identification of Donepezil and Natural
Polyphenols as Monoamine Oxidase-B (MAO-B) Inhibitors: A Structure-Based
Virtual Screening, ADMET Profiling, and Protein-Ligand Interaction Study.
Dr. M.G.R. Educational and Research Institute, Chennai, India.
```

## Author

Thanuhith J. | Department of Pharmacy
Dr. M.G.R. Educational and Research Institute, Chennai, Tamil Nadu, India

## License

MIT License - free to use, modify, and distribute with attribution.
