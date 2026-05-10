# primd vs primer3 — Comprehensive Validation

Generated: 2026-05-10
Test sequences: 8 (52% GC, 57% GC, 60% GC, 44% GC, 64% GC, 38% GC, 29% GC, 66% GC)
Windows per sequence: up to 3 × 400 bp non-overlapping
Total windows: 24

## A: Condition A: 50 mM NaCl, 0 mM Mg²⁺ (isolates NN accuracy)

**17/24** windows where both tools returned primers (70.8%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.931 | 0.879 |
| RMSE (°C) | 1.148 | 1.104 |
| Max Δ (°C) | 2.178 | 2.239 |
| Pearson r (all Tms) | 0.442 | — |
| Systematic bias (primd − p3, mean) | -0.623°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 20.7 bp | 21.5 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 0% |

## B: Condition B: 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs (typical PCR buffer)

**19/24** windows where both tools returned primers (79.2%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.684 | 0.568 |
| RMSE (°C) | 0.896 | 0.745 |
| Max Δ (°C) | 2.03 | 2.124 |
| Pearson r (all Tms) | 0.3233 | — |
| Systematic bias (primd − p3, mean) | -0.09°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 19.9 bp | 19.5 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 5.3% |

## Per-sequence Tm MAE
| Sequence | GC% | Cond A Fwd MAE | Cond A Rev MAE | Cond B Fwd MAE | Cond B Rev MAE |
|----------|-----|----------------|----------------|----------------|----------------|
| pUC19 | 52% | 1.873 | 1.562 | 1.231 | 0.527 |
| ecoli_rpoB | 57% | 0.557 | 0.186 | 0.701 | 0.302 |
| human_GAPDH | 60% | 1.255 | 1.272 | 1.097 | 1.143 |
| human_TP53_exon5 | 44% | 0.288 | 0.85 | 0.567 | 0.553 |
| MTB_IS6110 | 64% | 0.629 | 0.971 | 0.068 | 0.667 |
| SARS2_ORF1a | 38% | 1.025 | 0.34 | 0.687 | 0.68 |
| Pfalciparum_MSP1 | 29% | — | — | — | — |
| Streptomyces_actII | 66% | 1.068 | 1.076 | 0.371 | 0.328 |

## Interpretation notes
- **Tm MAE and Pearson r here compare Tm of *different* primer sequences** (each tool selects its own). This is NOT a thermodynamic accuracy metric. For that see validate-tm-direct.mjs: MAE=0.002°C, r=1.000 on identical sequences.
- **Positional agreement is expected to be low**: primer3 and primd use different penalty weightings. Position disagreement does not indicate error — both produce primers that amplify the target.
- **Plasmodium MSP1 (29% GC): 0 pairs from both tools** — correct behaviour. At 29% GC, even a 27 bp primer has Tm ≈ 52°C, below both tools' minimum Tm (57°C). Users targeting AT-rich organisms should lower tmTarget to 52–55°C.
- **Systematic bias**: the per-tool bias reflects different primer sequences selected (different length, GC%, position) — not a thermodynamic error.

