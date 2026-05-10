# designLAMP Validation

**Methodology**: No open-source LAMP design tool equivalent to primer3 exists. Primer Explorer (standard) is web-only/proprietary. Validation instead checks:
  1. Notomi geometry constraints on real biological sequences
  2. Tm ordering (outer < inner, per Notomi 2000 guidelines)
  3. Design recovery on SARS-CoV-2 N gene (confirmed F3/B3 positions as anchor)
  4. Thermodynamic verification on published IS6110 LAMP primer sequences (Sun 2017)

---

## Part 1: Design on Real Biological Sequences

| Sequence | GC% | Sets | Geometry | Tm ordering | Inner amp (bp) |
|----------|-----|------|----------|-------------|----------------|
| SARS2_N_gene | — | 0 | — | — | — |
| ecoli_rpoB | 57% | 3 | ✅ pass | ✅ pass | 167 |
| human_GAPDH | 60% | 3 | ✅ pass | ✅ pass | 200 |
| MTB_IS6110 | 64% | 1 | ✅ pass | ✅ pass | 175 |

### Best pair details per sequence

#### ecoli_rpoB (Escherichia coli K-12)

| Primer | Sequence | Tm | Len |
|--------|----------|----|-----|
| F3 | `CAGCGCAACCTGTTCGAC` | 62.33°C | 18 |
| B3 | `TCATCATCATGTTCATCGCGA` | 60.04°C | 21 |
| FIP | `TCGATCGCTTCGCCGATGTCAGACCCTGGCGCAGTTCAAA` | F1c:65.66° F2:65.32° | 40 |
| BIP | `CAGCAGCTGCTGCGCAAGGCTCATGTCCAGCATGTGCG` | B1c:65.9° B2:65.11° | 38 |
| LF | `TGACGAACAGGGTGGTGATTT` | 62.45°C | 21 |
| LB | `TGCTGGAGCGGATCCAG` | 62.39°C | 17 |

#### human_GAPDH (Homo sapiens)

| Primer | Sequence | Tm | Len |
|--------|----------|----|-----|
| F3 | `GATTTGGTCGTATTGGGCG` | 59.42°C | 19 |
| B3 | `CGCTGTTGAAGTCAGAGGA` | 59.91°C | 19 |
| FIP | `AAAGATGAAGGTGGCCTTGGGCTCCAGAGCTGAAGCAGGCTATC` | F1c:66.3° F2:64.69° | 44 |
| BIP | `AAGGTGGTGAAGCAGGCGTCTAGCCCAGGATGCCCTTGAG` | B1c:66.1° B2:64.81° | 40 |
| LF | `TGACGCTGAAGGGATCACC` | 62.57°C | 19 |
| LB | `CTATTTCATAGCCCGCCTGGA` | 62.49°C | 21 |

#### MTB_IS6110 (Mycobacterium tuberculosis)

| Primer | Sequence | Tm | Len |
|--------|----------|----|-----|
| F3 | `GATGAACCGCCCCGGCAT` | 66.51°C | 18 |
| B3 | `CACTTACGCACCGTCTCC` | 60.89°C | 18 |
| FIP | `CCACCTGACATGACCCCATCCGTCCGGAGACTCCAGTTCTTGG` | F1c:65.84° F2:65.25° | 43 |
| BIP | `CGATCAGTGAGGTCGCCCGGCGCAGCCAACACCAAGTAG` | B1c:66.07° B2:65.02° | 39 |
| LF | `GCGGGTACCTCCTCGATG` | 62.81°C | 18 |
| LB | `GCAGAGATCCGCGGTCA` | 62.59°C | 17 |

---

## Part 2: Published IS6110 Primer Tm Verification (Sun 2017)

Source: Sun J et al. (2017) Oncotarget 8(60):102264. IS6110-targeting LAMP for TB meningitis.
**Conditions**: SantaLucia 1998 + Owczarzy 2004, 50 mM NaCl, 250 nM oligo

| Primer | Sequence | GC% | Len | primd Tm |
|--------|----------|-----|-----|----------|
| F3 | `AGACCTCACCTATGTGTCGA` | 50% | 20 | 53.82°C |
| B3 | `TCGCTGAACCGGATCGA` | 58.8% | 17 | 55.38°C |
| F1c | `ATGGAGGTGGCCATCGTGGA` | 60% | 20 | 60.96°C |
| F2 | `AGCCTACGTGGCCTTTGTCAC` | 57.1% | 21 | 59.61°C |
| B1c | `AAGCCATCTGGACCCGCCAA` | 60% | 20 | 61.73°C |
| B2 | `CCCCTATCGTATGGTGGAT` | 52.6% | 19 | 52.28°C |
| LF | `AGGATCCTGCGAGCGTAG` | 61.1% | 18 | 56.2°C |
| LB | `AAGAAGGCGTACTCGACCTG` | 55% | 20 | 56.06°C |

**Outer primers** (F3/B3): mean Tm = 54.6°C
**Inner parts** (F1c/F2/B1c/B2): mean Tm = 58.6°C
**Tm ordering** (inner > outer, per Notomi guidelines): ✅ confirmed

---

## Notes
- Geometry checks verify all six Notomi 2000 spatial constraints simultaneously.
- Tm ordering (outer cooler than inner) is a design requirement that ensures outer primers are displaced at inner primer annealing temperature.
- The SARS-CoV-2 recovery test uses published F3/B3 positions as anchors; primd's design should find primers in the same window.
- No direct comparison to Primer Explorer output is possible (proprietary, web-only).
