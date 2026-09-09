.PHONY: install lint test reproduce clean

install:
	python3 -m pip install -U pip
	pip install -e ".[dev]"

lint:
	ruff check .

test:
	MPLBACKEND=Agg pytest -q

reproduce:
	bash scripts/reproduce.sh

clean:
	rm -rf data/raw/* data/processed/*
	rm -rf out/*.json out/*.log demo/outputs/*
	touch data/raw/.gitkeep data/processed/.gitkeep demo/outputs/.gitkeep
