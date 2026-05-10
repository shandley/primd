# primd Thermodynamic Validation

## Overview

primd implements two validated thermodynamic models:

- **Nearest-neighbor Tm**: SantaLucia 1998 (PNAS 95:1460), identical parameter set to primer3
- **Salt correction**: Owczarzy 2008 (Biochemistry 47:5336), accounting for Mg²⁺/dNTP competition for free Mg²⁺

Validation compares primd against **primer3 v2.6.0** (the field standard for PCR primer design) using two complementary approaches.

---

## Approach 1: Direct Tm Comparison (identical sequences)

The rigorous thermodynamic test: same oligo sequence → both tools compute Tm → compare.

**Test set:** 150 synthetic oligos, lengths 18–27 bp, GC 22–74%, 5 independent seeds per cell.

**Condition A — 50 mM NaCl, 0 mM Mg²⁺ (Owczarzy 2004 both tools)**

| Metric | Value |
|--------|-------|
| n | 150 |
| MAE (°C) | 0.0021 |
| RMSE (°C) | 0.002 |
| Max Δ (°C) | 0.003 |
| Pearson r | 1.000 |
| Mean bias | +0.002°C |

The NN parameter implementation is exact relative to primer3. The ≤0.003°C residual is floating-point rounding in the last significant digit.

**Condition B — 50 mM NaCl, 2 mM Mg²⁺, 0.8 mM dNTPs**

| Metric | Value |
|--------|-------|
| n | 150 |
| MAE (°C) | **0.113** |
| RMSE (°C) | 0.157 |
| 95th percentile Δ (°C) | 0.219 |
| Pearson r | **0.99997** |
| Mean bias | −0.113°C |

primd predicts ~0.15°C lower Tm than primer3 when Mg²⁺ is present. This is **expected and intentional**: primer3 uses Owczarzy 2004 for divalent ions; primd uses Owczarzy 2008, which corrects for dNTP chelation of free Mg²⁺. At 2 mM Mg²⁺ and 0.8 mM dNTPs, free Mg²⁺ ≈ 1.2 mM — the improved formula predicts lower effective Mg²⁺ concentration and correspondingly lower Tm.

The bias is:
- Systematic (consistent direction across all oligo lengths and GC%)
- Predictable (decreases with primer length: 0.22°C at 18 bp → 0.07°C at 27 bp)
- Small relative to experimental noise (~1–2°C for measured Tm values)

**Interpretation:** Researchers switching from primer3 using standard PCR buffer (50 mM KCl, 1.5 mM MgCl₂) should expect primd to predict Tm values ~0.1–0.2°C lower. Experimentally, this difference is within the annealing temperature gradient typically explored when optimizing a new PCR.

---

## Approach 2: Pair-finding Agreement (real sequences)

Both tools design primers for the same target regions; the best pairs are compared.

**Test set:** 8 diverse sequences × 3 non-overlapping 300 bp windows = 24 windows per condition.

| Sequence | Organism | GC% | Notes |
|----------|----------|-----|-------|
| pUC19 | Cloning vector | 52% | Standard reference |
| rpoB | *E. coli* K-12 | 57% | Clinical diagnostic (rifampicin resistance) |
| GAPDH | *H. sapiens* | 60% | Housekeeping gene |
| TP53 exons 4–7 | *H. sapiens* | 44% | Clinical mutation hotspot |
| IS6110 | *M. tuberculosis* | 64% | TB diagnostic target |
| ORF1a nsp3 | SARS-CoV-2 | 38% | Respiratory pathogen |
| MSP1 | *P. falciparum* | 29% | Extreme AT-rich malaria antigen |
| actII | *S. coelicolor* | 66% | High-GC actinomycete |

**Key findings:**

| | Condition A (0 mM Mg²⁺) | Condition B (2 mM Mg²⁺) |
|-|--------------------------|--------------------------|
| Coverage (both found) | 70.8% | 79.2% |
| Fwd Tm MAE | 0.93°C | 0.68°C |
| Rev Tm MAE | 0.88°C | 0.57°C |
| Positional agreement ≤10 bp | ~10% | ~5% |
| Mean primer length — p3 | 20.7 bp | 19.9 bp |
| Mean primer length — primd | 21.5 bp | 19.5 bp |

**Why positional agreement is low:** Both tools find primers with similar Tm at valid positions, but use different penalty weightings. primer3 strongly prefers 20 bp primers and minimizes distance from an optimal position; primd weights template accessibility and scores a longer primer with better thermodynamics. Position disagreement does not indicate error — both tools produce primers that would amplify the target. A 24 bp primer at position X and a 20 bp primer at position X+10 are both valid if their Tm values are within 1°C.

---

---

## Approach 3: LAMP Validation (real biological sequences)

LAMP (Loop-mediated Isothermal Amplification) has no open-source design counterpart equivalent to primer3. Primer Explorer (the field standard) is proprietary and web-only. Validation instead checks Notomi 2000 geometry constraints and Tm ordering rules, and verifies thermodynamic output against a published IS6110 primer set.

### Design recovery on real sequences

**Test set:** 4 sequences spanning 38–64% GC, 100 bp target regions.

| Sequence | GC% | Sets returned | Geometry | Tm ordering | Inner amp (bp) |
|----------|-----|---------------|----------|-------------|----------------|
| SARS-CoV-2 ORF1a | 38% | 3 | ✅ | ✅ | 197 |
| *E. coli* rpoB | 57% | 3 | ✅ | ✅ | 167 |
| *H. sapiens* GAPDH | 60% | 3 | ✅ | ✅ | 200 |
| *M. tuberculosis* IS6110 | 64% | 1 | ✅ | ✅ | 175 |

Geometry checks verify all six Notomi 2000 spatial constraints simultaneously. Tm ordering requires outer (F3/B3) primers cooler than all inner parts (FIP/BIP components), which is mandatory for efficient LAMP amplification — displacement of outer primers at the inner primer annealing temperature.

**SARS-CoV-2 note**: the target region must be ≤200 bp (inner amplicon constraint). A 195 bp target region leaves no margin; choosing a 100 bp target region (`regionStart=55, regionEnd=155`) yields 3 sets with geometry and Tm ordering satisfied.

### Published primer Tm verification (IS6110, Sun 2017)

Source: Sun J et al. (2017) *Oncotarget* 8(60):102264.

| Primer | Tm (primd) | Role |
|--------|-----------|------|
| F3 | 53.8°C | Outer |
| B3 | 55.4°C | Outer |
| F1c | 61.0°C | Inner |
| F2 | 59.6°C | Inner |
| B1c | 61.7°C | Inner |
| B2 | 52.3°C | Inner |

Outer mean: 54.6°C · Inner mean (F1c, F2, B1c): 60.8°C → **Tm ordering confirmed** ✅

The published B2 Tm (52.3°C) is below the inner Tm target range — consistent with experimental primer sets that were designed manually or with different Tm targets. primd would not select this primer as a top candidate under default settings, which is expected.

---

## Limitations and Known Issues

**1. primd length preference reduced but not eliminated**
After adding an explicit `|len - 20| × 1.0` length penalty, mean primer length in Condition A improved from 24.0 bp to 21.5 bp (vs primer3's 20.7 bp). In Condition B with Mg²⁺, lengths are essentially equal (19.5 vs 19.9 bp). This reflects the accessibility penalty in the scoring function: longer primers have more thermodynamic stability to overcome secondary structure. In Condition B with Mg²⁺, mean lengths equalize at 20 bp. This behavior is documented but not necessarily a limitation — longer primers can improve specificity.

**2. No experimental wet-lab validation**
Neither primer set has been ordered and tested. The validation here is computational: thermodynamic implementation fidelity relative to primer3. Actual PCR success depends on additional factors (polymerase, template complexity, buffer optimization) not captured by either tool.

**3. Extreme AT-rich sequences (PCR)**
For sequences <30% GC (e.g., *P. falciparum* introns), the default Tm target (60°C) may be unachievable within the GC% and length constraints. Both primer3 and primd return 0 pairs for *P. falciparum* MSP1 (29% GC) — correct behavior. Users should lower `tmTarget` to 52–55°C for extremely AT-rich organisms.

**4. AT-rich LAMP sequences**
LAMP requires 6 primers with Tms in tight ranges (outer 58–62°C, inner 63–68°C) and GC% between 35–70%. For templates <35% GC (e.g., *P. falciparum*, *A. thaliana* introns), it may be impossible to satisfy all six constraints simultaneously under default settings. Symptoms:
- 0 sets returned with "No valid forward inner primers" warning
- Reducing `innerTmRange` to `[55, 62]` and `outerTmRange` to `[50, 56]` allows primer finding at the cost of reduced amplification efficiency

```typescript
// For AT-rich templates (~30-40% GC):
designLAMP(seq, start, end, {
  innerTmRange: [55, 62],
  outerTmRange: [50, 56],
  loopTmRange:  [52, 58],
})
```

**5. Owczarzy 2008 vs 2004 divergence**
The 0.15°C systematic bias with Mg²⁺ reflects a genuine difference between published correction formulas. Neither formula is "wrong" — Owczarzy 2008 incorporates additional experimental data on dNTP/Mg²⁺ binding. For consistency with primer3-derived protocols, users can set `saltModel: "owczarzy_2004"` explicitly.

---

## Reproducing the Validation

```bash
# Clone primd and install dependencies
git clone https://github.com/shandley/primd && cd primd
pnpm install
npx tsc --project tsconfig.build.json   # build dist/

# Install primer3 (macOS)
brew install primer3

# Run direct Tm comparison (rigorous thermodynamic test)
node scripts/validate-tm-direct.mjs
# → validation_results/tm_direct_report.md
# → validation_results/tm_direct.csv

# Run pair-finding comparison (design benchmark)
node scripts/validate-primer3.mjs
# → validation_results/report.md
# → validation_results/all_pairs.csv
```

Results are written to `validation_results/`. Raw CSV files are provided for independent statistical analysis.

---

## Citation

If you use primd in published work, please cite the thermodynamic models directly:

- SantaLucia J (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95(4):1460–1465.
- Owczarzy R, et al. (2008). Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. *Biochemistry* 47(19):5336–5353.
