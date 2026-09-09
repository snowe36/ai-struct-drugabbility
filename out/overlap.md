# Pocket Atlas — NMR dynamics vs cryptic sites

Dyna-1 predicts *where* micro–millisecond motion occurs. It does not
generate molecular dynamics. 100 ns public atlases are the wrong
ensemble for cryptic opening; VP35 FAH/FAST is the long-MD arm.

## TEM-1 β-lactamase cryptic horn pocket

- Apo `1BTL` site recovered: **True** (1561 Å³)
- Holo `1PZO` site recovered: **True** (1161 Å³)
- NMR ∩ cryptic lining: 0/18 (enrichment 0.00)
- NMR ∩ control site: 6/6 (enrichment 8.77)
- Odds ratio (cryptic vs control): inf
- Prior: `nmr_literature`

TEM-1 horn lining is mostly buried hydrophobic core, not the Ω-loop.
Savard/Gagné μs–ms exchange is reported at the Ω-loop and active-site
vicinity. If cryptic enrichment is low, that is a result: NMR dynamics
are a prior for *motion*, not a pocket oracle.

## KRAS switch-II cryptic pocket

- Apo `5V9U` site recovered: **True** (909 Å³)
- Holo `6OIM` site recovered: **True** (602 Å³)
- NMR ∩ cryptic lining: 19/26 (enrichment 4.07)
- NMR ∩ control site: 3/16 (enrichment 1.04)
- Odds ratio (cryptic vs control): 11.76
- Prior: `nmr_literature`

KRAS switch-I/II carry most published μs–ms NMR signal and also line
the sotorasib site. Enrichment here is the expected positive control
for the same question.

## Limitations

- Literature NMR residue lists are curated subsets, not the full RelaxDB dump.
- Dyna-1 weights are optional; when absent the prior is those literature labels.
- Pocket volumes are a geometric analogue, not SiteMap Dscore.
- VP35 trajectories are not downloaded by the demo (multi-GB Zenodo archives).
