# primd vs primer3 — Comprehensive Validation

Generated: 2026-05-10
Test sequences: 8 (52% GC, 57% GC, 60% GC, 44% GC, 64% GC, 38% GC, 29% GC, 66% GC)
Windows per sequence: up to 3 × 400 bp non-overlapping
Total windows: 24

## A: Condition A: 50 mM NaCl, 0 mM Mg²⁺ (isolates NN accuracy)

**16/24** windows where both tools returned primers (66.7%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.588 | 0.811 |
| RMSE (°C) | 0.769 | 1.017 |
| Max Δ (°C) | 2.063 | 1.947 |
| Pearson r (all Tms) | 0.27 | — |
| Systematic bias (primd − p3, mean) | -0.297°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 20.8 bp | 24 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 12.5% |

## B: Condition B: 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs (typical PCR buffer)

**17/24** windows where both tools returned primers (70.8%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.556 | 0.528 |
| RMSE (°C) | 0.697 | 0.641 |
| Max Δ (°C) | 1.365 | 1.158 |
| Pearson r (all Tms) | 0.5262 | — |
| Systematic bias (primd − p3, mean) | 0.019°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 19.9 bp | 19.6 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 5.9% |

## Per-sequence Tm MAE
| Sequence | GC% | Cond A Fwd MAE | Cond A Rev MAE | Cond B Fwd MAE | Cond B Rev MAE |
|----------|-----|----------------|----------------|----------------|----------------|
| pUC19 | 52% | 0.178 | 0.339 | 0.905 | 0.62 |
| ecoli_rpoB | 57% | 0.57 | 0.322 | 0.701 | 0.302 |
| human_GAPDH | 60% | 0.604 | 0.641 | 0.206 | 0.344 |
| human_TP53_exon5 | 44% | 0.482 | 0.719 | 0.511 | 0.694 |
| MTB_IS6110 | 64% | 1.198 | 1.161 | 0.068 | 0.667 |
| SARS2_ORF1a | 38% | 0.331 | 1.79 | 0.687 | 0.352 |
| Pfalciparum_MSP1 | 29% | — | — | — | — |
| Streptomyces_actII | 66% | 0.603 | 0.847 | 0.701 | 0.838 |

## Interpretation notes
- **Positional agreement** is expected to be low — both tools find valid primers at similar Tm but use different penalty weightings. primer3 prefers 20bp; primd scores accessibility and allows up to 27bp. Position disagreement does not indicate error.
- **Tm MAE < 0.5°C** is within experimental measurement noise for standard oligonucleotides.
- **Systematic bias**: positive bias means primd predicts higher Tm than primer3. This is expected if longer primers are selected (more NN stacking).
- **Pearson r > 0.99** indicates the thermodynamic model is well-calibrated.

