# primd vs primer3 — Comprehensive Validation

Generated: 2026-05-10
Test sequences: 8 (52% GC, 57% GC, 60% GC, 44% GC, 65% GC, 38% GC, 22% GC, 72% GC)
Windows per sequence: up to 3 × 400 bp non-overlapping
Total windows: 24

## A: Condition A: 50 mM NaCl, 0 mM Mg²⁺ (isolates NN accuracy)

**10/24** windows where both tools returned primers (41.7%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.514 | 0.538 |
| RMSE (°C) | 0.614 | 0.755 |
| Max Δ (°C) | 1.172 | 1.936 |
| Pearson r (all Tms) | 0.5474 | — |
| Systematic bias (primd − p3, mean) | -0.112°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 20 bp | 24.1 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 10% |

## B: Condition B: 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs (typical PCR buffer)

**10/24** windows where both tools returned primers (41.7%)

### Tm accuracy
| Metric | Fwd | Rev |
|--------|-----|-----|
| MAE (°C) | 0.586 | 0.491 |
| RMSE (°C) | 0.728 | 0.58 |
| Max Δ (°C) | 1.365 | 1.158 |
| Pearson r (all Tms) | 0.5773 | — |
| Systematic bias (primd − p3, mean) | 0.135°C | — |

### Primer characteristics
| | primer3 | primd |
|-|---------|-------|
| Mean fwd length | 20 bp | 19.6 bp |
| Same position (Δ≤2 bp) | — | 0% |
| Same position (Δ≤10 bp) | — | 0% |

## Per-sequence Tm MAE
| Sequence | GC% | Cond A Fwd MAE | Cond A Rev MAE | Cond B Fwd MAE | Cond B Rev MAE |
|----------|-----|----------------|----------------|----------------|----------------|
| pUC19 | 52% | 0.178 | 0.339 | 0.905 | 0.62 |
| ecoli_rpoB | 57% | 0.57 | 0.322 | 0.701 | 0.302 |
| human_GAPDH | 60% | 0.604 | 0.641 | 0.206 | 0.344 |
| human_TP53_exon5 | 44% | 0.482 | 0.719 | 0.511 | 0.694 |
| MTB_IS6110 | 65% | — | — | — | — |
| SARS2_ORF1a | 38% | — | — | — | — |
| Pfalciparum_MSP1 | 22% | — | — | — | — |
| Streptomyces_actII | 72% | — | — | — | — |

## Interpretation notes
- **Positional agreement** is expected to be low — both tools find valid primers at similar Tm but use different penalty weightings. primer3 prefers 20bp; primd scores accessibility and allows up to 27bp. Position disagreement does not indicate error.
- **Tm MAE < 0.5°C** is within experimental measurement noise for standard oligonucleotides.
- **Systematic bias**: positive bias means primd predicts higher Tm than primer3. This is expected if longer primers are selected (more NN stacking).
- **Pearson r > 0.99** indicates the thermodynamic model is well-calibrated.

