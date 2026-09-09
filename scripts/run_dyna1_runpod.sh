#!/usr/bin/env bash
# Runs on the Runpod box. Do not commit weights or scored PDBs from here.
set -euo pipefail
cd /workspace
mkdir -p pdbs out Dyna-1

if [ ! -f /workspace/Dyna-1/dyna1.py ]; then
  git clone --depth 1 https://github.com/WaymentSteeleLab/Dyna-1.git /workspace/Dyna-1
fi
cd /workspace/Dyna-1

python -m pip install -U pip
if [ -f requirements.txt ]; then
  python -m pip install -r requirements.txt
fi
python -m pip install huggingface_hub MDAnalysis transformers pandas

mkdir -p model/weights
python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download("gelnesr/Dyna-1", local_dir="model/weights", local_dir_use_symlinks=False)
print("downloaded gelnesr/Dyna-1")
PY

ls -la model/weights
python dyna1.py --help || true

for spec in \
  "tem1_horn /workspace/pdbs/tem1_horn_apo_1BTL.pdb A" \
  "kras_switch2 /workspace/pdbs/kras_switch2_apo_5V9U.pdb A"
do
  set -- $spec
  name=$1; pdb=$2; chain=$3
  echo "==== Dyna-1 $name ===="
  python dyna1.py --pdb "$pdb" --chain "$chain" --name "$name" --use_pdb_seq --write_to_pdb --save_dir /workspace/out
done
ls -la /workspace/out
echo DONE
