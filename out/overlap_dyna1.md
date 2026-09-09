# Pocket Atlas — NMR dynamics vs cryptic sites

Dyna-1 predicts *where* micro–millisecond motion occurs. It does not
generate molecular dynamics. 100 ns public atlases are the wrong
ensemble for cryptic opening; VP35 FAH/FAST is the long-MD arm.

## TEM-1 β-lactamase cryptic horn pocket

- Apo `1BTL` site recovered: **True** (clearance 3.45 Å; volume 1561 Å³ is secondary)
- Holo `1PZO` site recovered: **True** (clearance 4.02 Å; volume 1161 Å³ is secondary)
- NMR ∩ cryptic lining: 1/18 (enrichment 0.29)
- NMR ∩ control site: 0/6 (enrichment 0.00)
- Odds ratio (cryptic vs control): inf
- Prior: `dyna1` (Dyna-1 vs literature when both exist)

TEM-1 horn lining is mostly buried hydrophobic core, not the Ω-loop.
Savard/Gagné μs–ms exchange is reported at the Ω-loop and active-site
vicinity. If cryptic enrichment is low, that is a result: NMR dynamics
are a prior for *motion*, not a pocket oracle.

## KRAS switch-II cryptic pocket

- Apo `5V9U` site recovered: **True** (clearance 4.21 Å; volume 909 Å³ is secondary)
- Holo `6OIM` site recovered: **True** (clearance 4.40 Å; volume 602 Å³ is secondary)
- NMR ∩ cryptic lining: 7/26 (enrichment 1.32)
- NMR ∩ control site: 4/16 (enrichment 1.23)
- Odds ratio (cryptic vs control): 1.11
- Prior: `dyna1` (Dyna-1 vs literature when both exist)

KRAS switch-I/II carry most published μs–ms NMR signal and also line
the sotorasib site. Enrichment here is the expected positive control
for the same question.

## Limitations

- `--prior literature` is the curated YAML (CI/demo, offline).
- `--prior relaxdb` uses RelaxDB-CPMG. KRAS is in that set; TEM-1 is not (BLAC_CPMG is Mtb BlaC P9WKD3, not TEM-1).
- Dyna-1 weights are optional; captions must say Dyna-1 vs literature.
- Headline geometry is seed clearance, not matched-site volume.
- VP35 trajectories are not downloaded by the demo (multi-GB Zenodo archives).
