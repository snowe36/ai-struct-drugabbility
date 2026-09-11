#!/usr/bin/env bash
# CPU pod: download Zaire FAH, stride, CA 225–295 occupancy.
# Expects /workspace/pocket-atlas already copied. Do not copy trajectories back.
set -euo pipefail
cd /workspace
if [ ! -d /workspace/pocket-atlas ]; then
  echo "missing /workspace/pocket-atlas"
  exit 1
fi

python3 -m venv /workspace/venv
# shellcheck disable=SC1091
source /workspace/venv/bin/activate
python -m pip install -U pip
python -m pip install -e "/workspace/pocket-atlas[md]"

cd /workspace/pocket-atlas
pocket-fetch-md --yes --extract --analyze
# FAH tarball is multi-GB; keep occupancy JSON only.
rm -f data/raw/vp35/*.tar.gz
rm -rf data/processed/vp35/extract
ls -la data/processed/vp35 || true
echo DONE
