"""Branch protection and local workflow. The GitHub repo settings are yours.

- Protect `main`: no direct pushes, no force-push, PRs required.
- Require CI (`ruff` + `pytest`) to merge.
- Feature branches; linear history preferred.
- Do not commit `data/raw` PDBs, Zenodo tarballs, ESM weights, or `.venv`.
- Do commit `out/figures/*.png` used in the README.
- Dyna-1 is optional. The demo must stay green without weights or FAH downloads.
"""
