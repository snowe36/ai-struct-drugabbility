#!/usr/bin/env bash
# GPU pod: Dyna-1 on VP35 3FKE, GNINA TEM-1/KRAS, ligand xTB strain.
# Expects /workspace/pocket-atlas (this repo) already copied.
set -euo pipefail
cd /workspace
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/mpl-cache
mkdir -p /workspace/bin /workspace/out /workspace/xtb

if [ ! -d /workspace/pocket-atlas ]; then
  echo "missing /workspace/pocket-atlas"
  exit 1
fi

python3 -m venv --system-site-packages /workspace/venv
# shellcheck disable=SC1091
source /workspace/venv/bin/activate
python -m pip install -U pip
python -m pip install -e "/workspace/pocket-atlas[dev,md]"

if [ ! -x /workspace/bin/gnina ]; then
  curl -L -o /workspace/bin/gnina \
    https://github.com/gnina/gnina/releases/download/v1.3/gnina
  chmod +x /workspace/bin/gnina
fi
export PATH="/workspace/bin:$PATH"

if [ ! -x /workspace/xtb/bin/xtb ]; then
  curl -L -o /tmp/xtb.tar.xz \
    https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz
  tar -xJf /tmp/xtb.tar.xz -C /workspace/xtb --strip-components=1
fi
export PATH="/workspace/xtb/bin:$PATH"

cd /workspace/pocket-atlas
pocket-fetch --case tem1_horn
pocket-fetch --case kras_switch2
pocket-fetch --case vp35_iid
pocket-prepare-crystal --case vp35_iid

# Dyna-1 on VP35 crystal
if [ ! -f /workspace/Dyna-1/dyna1.py ]; then
  git clone --depth 1 https://github.com/WaymentSteeleLab/Dyna-1.git /workspace/Dyna-1
fi
cd /workspace/Dyna-1
python -m pip install huggingface_hub MDAnalysis transformers pandas
if [ -f requirements.txt ]; then
  python -m pip install -r requirements.txt || true
fi
mkdir -p model/weights
python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download("gelnesr/Dyna-1", local_dir="model/weights", local_dir_use_symlinks=False)
print("downloaded gelnesr/Dyna-1")
PY
VP35_PDB=$(ls /workspace/pocket-atlas/data/processed/vp35_iid_crystal_3FKE.pdb)
python dyna1.py --pdb "$VP35_PDB" --chain A --name vp35_iid --use_pdb_seq --write_to_pdb --save_dir /workspace/out

cd /workspace/pocket-atlas
python - <<'PY'
from pathlib import Path
from pocket_atlas.dynamics.dyna1 import _scores_from_csv, save_scores
csv = Path("/workspace/out/vp35_iid-Dyna1.csv")
scores = _scores_from_csv(csv)
print("vp35", len(scores), "->", save_scores("vp35_iid", scores))
PY

echo "==== GNINA ===="
pocket-dock --cases tem1_horn kras_switch2 --prior literature --exhaustiveness 8

echo "==== xTB strain ===="
pocket-strain --cases tem1_horn kras_switch2 || true

ls -la data/processed /workspace/out data/processed/dock || true
echo DONE
