#!/bin/bash
# NeuroDock system tool installer
# Installs AutoDock Vina 1.2.5 and Open Babel
# Usage: bash scripts/install_tools.sh

set -e

echo "NeuroDock: Installing system tools..."
echo ""

if [[ "$OSTYPE" != "linux-gnu"* ]]; then
    echo "ERROR: This script requires Linux or WSL (Ubuntu)."
    echo "Windows users: run this inside WSL2."
    exit 1
fi

echo "Step 1: Installing Open Babel..."
sudo apt-get update -qq
sudo apt-get install -y openbabel
echo "Open Babel installed."

echo ""
echo "Step 2: Installing AutoDock Vina 1.2.5..."
wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /tmp/vina
chmod +x /tmp/vina
sudo mv /tmp/vina /usr/local/bin/vina
echo "AutoDock Vina installed."

echo ""
echo "Verifying..."
vina --version 2>&1 | head -1
obabel --version 2>&1 | head -1

echo ""
echo "All tools ready. Run: python3 run_pipeline.py"
