# primd vs primer3 — Direct Tm Comparison

**Methodology**: Identical oligo sequences computed by both tools. Differences reflect only the thermodynamic implementation.
**Test oligos**: 150 (lengths 18–27 bp, GC 22–74%, 5 seeds per cell)

## Condition A: 50 mM NaCl, 0 mM Mg²⁺

**n = 81 oligos**

| Metric | Value |
|--------|-------|
| MAE (°C) | **0.0021** |
| RMSE (°C) | 0.0022 |
| 95th percentile Δ (°C) | 0.003 |
| Max Δ (°C) | 0.003 |
| Mean bias (primd − p3) | 0.0021°C |
| Pearson r | 1 |

### By primer length
| Length | n | MAE | Max Δ | Bias |
|--------|---|-----|-------|------|
| 18 bp | 18 | 0.002 | 0.003 | 0.002 |
| 20 bp | 19 | 0.002 | 0.003 | 0.002 |
| 22 bp | 19 | 0.002 | 0.003 | 0.002 |
| 24 bp | 16 | 0.002 | 0.002 | 0.002 |
| 27 bp | 9 | 0.002 | 0.002 | 0.002 |

### By GC content
| GC target | n | MAE | Max Δ | Bias |
|-----------|---|-----|-------|------|
| 22% | 19 | 0.002 | 0.003 | 0.002 |
| 35% | 22 | 0.002 | 0.003 | 0.002 |
| 45% | 18 | 0.002 | 0.003 | 0.002 |
| 55% | 11 | 0.002 | 0.003 | 0.002 |
| 65% | 10 | 0.002 | 0.003 | 0.002 |
| 74% | 1 | 0.003 | 0.003 | 0.003 |

## Condition B: 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs

**n = 52 oligos**

| Metric | Value |
|--------|-------|
| MAE (°C) | **0.1131** |
| RMSE (°C) | 0.1228 |
| 95th percentile Δ (°C) | 0.198 |
| Max Δ (°C) | 0.219 |
| Mean bias (primd − p3) | -0.1131°C |
| Pearson r | 0.99997 |

### By primer length
| Length | n | MAE | Max Δ | Bias |
|--------|---|-----|-------|------|
| 18 bp | 12 | 0.147 | 0.219 | -0.147 |
| 20 bp | 15 | 0.134 | 0.199 | -0.134 |
| 22 bp | 12 | 0.102 | 0.148 | -0.102 |
| 24 bp | 9 | 0.075 | 0.134 | -0.075 |
| 27 bp | 4 | 0.054 | 0.073 | -0.054 |

### By GC content
| GC target | n | MAE | Max Δ | Bias |
|-----------|---|-----|-------|------|
| 22% | 18 | 0.065 | 0.098 | -0.065 |
| 35% | 18 | 0.11 | 0.137 | -0.11 |
| 45% | 12 | 0.16 | 0.179 | -0.16 |
| 55% | 4 | 0.203 | 0.219 | -0.203 |

## Notes
- Condition A uses Owczarzy 2004 in both tools. Agreement here reflects NN parameter parity.
- Condition B uses Owczarzy 2004 divalent (primer3) vs Owczarzy 2008 (primd). Larger Δ here is expected and reflects the different Mg²⁺ correction formulas, not an error.
- Pearson r > 0.999 on identical oligos confirms the NN implementation is correct.
