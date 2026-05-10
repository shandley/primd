# primd vs primer3 — Direct Tm Comparison

**Methodology**: Identical oligo sequences computed by both tools. Differences reflect only the thermodynamic implementation.
**Test oligos**: 150 (lengths 18–27 bp, GC 22–74%, 5 seeds per cell)

## Condition A: 50 mM NaCl, 0 mM Mg²⁺

**n = 20 oligos**

| Metric | Value |
|--------|-------|
| MAE (°C) | **0.0022** |
| RMSE (°C) | 0.0022 |
| 95th percentile Δ (°C) | 0.003 |
| Max Δ (°C) | 0.003 |
| Mean bias (primd − p3) | 0.0022°C |
| Pearson r | 1 |

### By primer length
| Length | n | MAE | Max Δ | Bias |
|--------|---|-----|-------|------|
| 18 bp | 5 | 0.003 | 0.003 | 0.003 |
| 20 bp | 4 | 0.002 | 0.003 | 0.002 |
| 22 bp | 5 | 0.002 | 0.002 | 0.002 |
| 24 bp | 4 | 0.002 | 0.002 | 0.002 |
| 27 bp | 2 | 0.002 | 0.002 | 0.002 |

### By GC content
| GC target | n | MAE | Max Δ | Bias |
|-----------|---|-----|-------|------|
| 45% | 2 | 0.002 | 0.002 | 0.002 |
| 55% | 7 | 0.002 | 0.002 | 0.002 |
| 65% | 10 | 0.002 | 0.003 | 0.002 |
| 74% | 1 | 0.003 | 0.003 | 0.003 |

## Condition B: 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs

**n = 14 oligos**

| Metric | Value |
|--------|-------|
| MAE (°C) | **0.1502** |
| RMSE (°C) | 0.1567 |
| 95th percentile Δ (°C) | 0.219 |
| Max Δ (°C) | 0.219 |
| Mean bias (primd − p3) | -0.1502°C |
| Pearson r | 0.99962 |

### By primer length
| Length | n | MAE | Max Δ | Bias |
|--------|---|-----|-------|------|
| 18 bp | 1 | 0.219 | 0.219 | -0.219 |
| 20 bp | 6 | 0.18 | 0.199 | -0.18 |
| 22 bp | 3 | 0.146 | 0.148 | -0.146 |
| 24 bp | 2 | 0.112 | 0.134 | -0.112 |
| 27 bp | 2 | 0.072 | 0.073 | -0.072 |

### By GC content
| GC target | n | MAE | Max Δ | Bias |
|-----------|---|-----|-------|------|
| 35% | 3 | 0.078 | 0.089 | -0.078 |
| 45% | 7 | 0.151 | 0.161 | -0.151 |
| 55% | 4 | 0.203 | 0.219 | -0.203 |

## Notes
- Condition A uses Owczarzy 2004 in both tools. Agreement here reflects NN parameter parity.
- Condition B uses Owczarzy 2004 divalent (primer3) vs Owczarzy 2008 (primd). Larger Δ here is expected and reflects the different Mg²⁺ correction formulas, not an error.
- Pearson r > 0.999 on identical oligos confirms the NN implementation is correct.
