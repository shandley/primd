# primd vs primer3 benchmark — pUC19 (L09137.2, 2686 bp)

**Conditions (both tools):** SantaLucia 1998 NN · Owczarzy 2004 salt correction · 50 mM NaCl · 0 mM Mg²⁺ · 250 nM oligo · Tm target 60°C · primer 18–27 bp · product 150–500 bp · GC 40–65%

> primer3 reports hairpin and dimer thresholds in °C (thermodynamic ΔG converted internally).  
> primd reports hairpin and heterodimer ΔG directly in kcal/mol.  
> Negative primd hairpin ΔG < −2 kcal/mol corresponds roughly to primer3 hairpin_th > 47°C.

## lacZ-alpha (200–400), 50% GC

**primer3** (5 pairs returned)  
**primd** (5 pairs returned)

### Best pair comparison

| | primer3 | primd |
|-|---------|-------|
| Fwd seq    | `GTGTCGGGGCTGGCTTAACT` | `GGCTTAACTATGCGGCATCAGAGC` |
| Rev seq    | `ACTGGAAAGCGGGCAGTGAG` | `CCGGCTCGTATGTTGTGTGGAATTG` |
| Fwd Tm     | 59.7°C | 60.0°C |
| Rev Tm     | 59.9°C | 60.0°C |
| Fwd Tm Δ60 | -0.26°C | 0.03°C |
| Rev Tm Δ60 | -0.06°C | 0.02°C |
| Fwd GC     | 60% | 54% |
| Rev GC     | 60% | 52% |
| Fwd len    | 20 bp | 24 bp |
| Rev len    | 20 bp | 25 bp |
| Product    | 483 bp | 389 bp |
| Fwd start  | 127 | 138 |
| Rev start  | 609 | 502 |
| Fwd hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Rev hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Pair 3′ dimer | 0.0°C (th) | -0.3 kcal/mol |

### Positional agreement
Fwd pos Δ11 bp, Rev pos Δ107 bp | Fwd Tm Δ0.29°C, Rev Tm Δ0.08°C

### All returned pairs — Tm summary
| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |
|------|-----------|-----------|-----------|-----------|------------|------------|
| 1 | 59.7 | 59.9 | 60.0 | 60.0 | 483 | 389 |
| 2 | 59.7 | 59.9 | 59.9 | 59.9 | 484 | 345 |
| 3 | 59.7 | 59.9 | 59.9 | 59.9 | 428 | 350 |
| 4 | 59.7 | 59.9 | 59.9 | 60.0 | 429 | 431 |
| 5 | 59.7 | 59.9 | 59.7 | 59.7 | 490 | 435 |

## bla/ampR (1800–2000), 55% GC

**primer3** (5 pairs returned)  
**primd** (5 pairs returned)

### Best pair comparison

| | primer3 | primd |
|-|---------|-------|
| Fwd seq    | `CGCTCACCGGCTCCAGATTT` | `CGCTCACCGGCTCCAGATTTATC` |
| Rev seq    | `CAACTCGGTCGCCGCATACA` | `GTGGCGCGGTATTATCCCGTATTG` |
| Fwd Tm     | 60.0°C | 59.8°C |
| Rev Tm     | 60.0°C | 59.8°C |
| Fwd Tm Δ60 | -0.03°C | -0.20°C |
| Rev Tm Δ60 | 0.02°C | -0.21°C |
| Fwd GC     | 60% | 57% |
| Rev GC     | 60% | 54% |
| Fwd len    | 20 bp | 23 bp |
| Rev len    | 20 bp | 24 bp |
| Product    | 452 bp | 490 bp |
| Fwd start  | 1773 | 1773 |
| Rev start  | 2224 | 2239 |
| Fwd hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Rev hairpin | 44.0°C (th) | 0.0 kcal/mol |
| Pair 3′ dimer | 2.8°C (th) | -0.4 kcal/mol |

### Positional agreement
Fwd pos Δ0 bp, Rev pos Δ15 bp | Fwd Tm Δ0.17°C, Rev Tm Δ0.24°C

### All returned pairs — Tm summary
| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |
|------|-----------|-----------|-----------|-----------|------------|------------|
| 1 | 60.0 | 60.0 | 59.8 | 59.8 | 452 | 490 |
| 2 | 60.0 | 60.0 | 59.8 | 59.8 | 451 | 491 |
| 3 | 60.0 | 60.1 | 59.9 | 59.9 | 473 | 422 |
| 4 | 60.0 | 60.2 | 59.6 | 59.6 | 471 | 330 |
| 5 | 60.0 | 59.8 | 59.8 | 59.9 | 490 | 423 |

## ori/AT-rich (700–900), 45% GC

**primer3** (5 pairs returned)  
**primd** (5 pairs returned)

### Best pair comparison

| | primer3 | primd |
|-|---------|-------|
| Fwd seq    | `CTCACTGCCCGCTTTCCAGT` | `GCCGGAAGCATAAAGTGTAAAGCCTG` |
| Rev seq    | `CAGGGTCGGAACAGGAGAGC` | `GGTTTCGCCACCTCTGACTTGAG` |
| Fwd Tm     | 59.9°C | 59.9°C |
| Rev Tm     | 60.1°C | 59.9°C |
| Fwd Tm Δ60 | -0.06°C | -0.12°C |
| Rev Tm Δ60 | 0.13°C | -0.13°C |
| Fwd GC     | 60% | 50% |
| Rev GC     | 65% | 57% |
| Fwd len    | 20 bp | 26 bp |
| Rev len    | 20 bp | 23 bp |
| Product    | 415 bp | 412 bp |
| Fwd start  | 590 | 522 |
| Rev start  | 1004 | 911 |
| Fwd hairpin | 36.9°C (th) | 0.0 kcal/mol |
| Rev hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Pair 3′ dimer | 12.9°C (th) | -0.3 kcal/mol |

### Positional agreement
Fwd pos Δ68 bp, Rev pos Δ93 bp | Fwd Tm Δ0.06°C, Rev Tm Δ0.25°C

### All returned pairs — Tm summary
| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |
|------|-----------|-----------|-----------|-----------|------------|------------|
| 1 | 59.9 | 60.1 | 59.9 | 59.9 | 415 | 412 |
| 2 | 59.9 | 60.1 | 59.8 | 59.8 | 416 | 333 |
| 3 | 59.9 | 60.1 | 60.2 | 60.2 | 414 | 329 |
| 4 | 59.9 | 60.1 | 60.2 | 60.2 | 415 | 484 |
| 5 | 59.9 | 60.1 | 59.8 | 59.9 | 469 | 492 |

## MCS region (390–570), 52% GC

**primer3** (5 pairs returned)  
**primd** (5 pairs returned)

### Best pair comparison

| | primer3 | primd |
|-|---------|-------|
| Fwd seq    | `AAGGGGGATGTGCTGCAAGG` | `TGCTGCAAGGCGATTAAGTTGGG` |
| Rev seq    | `AGCGGAAGAGCGCCCAATAC` | `CGCCTTTGAGTGAGCTGATACCG` |
| Fwd Tm     | 60.2°C | 59.9°C |
| Rev Tm     | 60.0°C | 59.9°C |
| Fwd Tm Δ60 | 0.23°C | -0.12°C |
| Rev Tm Δ60 | 0.01°C | -0.12°C |
| Fwd GC     | 60% | 52% |
| Rev GC     | 60% | 57% |
| Fwd len    | 20 bp | 23 bp |
| Rev len    | 20 bp | 23 bp |
| Product    | 378 bp | 437 bp |
| Fwd start  | 315 | 325 |
| Rev start  | 692 | 739 |
| Fwd hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Rev hairpin | 41.0°C (th) | 0.0 kcal/mol |
| Pair 3′ dimer | 0.0°C (th) | -1.7 kcal/mol |

### Positional agreement
Fwd pos Δ10 bp, Rev pos Δ47 bp | Fwd Tm Δ0.35°C, Rev Tm Δ0.13°C

### All returned pairs — Tm summary
| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |
|------|-----------|-----------|-----------|-----------|------------|------------|
| 1 | 60.2 | 60.0 | 59.9 | 59.9 | 378 | 437 |
| 2 | 59.7 | 60.0 | 59.9 | 59.9 | 432 | 442 |
| 3 | 60.2 | 59.9 | 59.9 | 60.0 | 295 | 431 |
| 4 | 59.7 | 59.9 | 60.2 | 60.2 | 483 | 463 |
| 5 | 60.2 | 59.9 | 60.2 | 60.2 | 296 | 334 |

## rep GC-rich (1100–1300), 62% GC

**primer3** (5 pairs returned)  
**primd** (5 pairs returned)

### Best pair comparison

| | primer3 | primd |
|-|---------|-------|
| Fwd seq    | `CTCTCCTGTTCCGACCCTGC` | `GGTATCTCAGTTCGGTGTAGGTCGTTC` |
| Rev seq    | `CCAGCGGTGGTTTGTTTGCC` | `TCAGCAGAGCGCAGATACCAAATACTG` |
| Fwd Tm     | 60.1°C | 60.0°C |
| Rev Tm     | 60.0°C | 60.0°C |
| Fwd Tm Δ60 | 0.13°C | -0.03°C |
| Rev Tm Δ60 | -0.05°C | -0.01°C |
| Fwd GC     | 65% | 52% |
| Rev GC     | 60% | 48% |
| Fwd len    | 20 bp | 27 bp |
| Rev len    | 20 bp | 27 bp |
| Product    | 411 bp | 261 bp |
| Fwd start  | 986 | 1075 |
| Rev start  | 1396 | 1309 |
| Fwd hairpin | 0.0°C (th) | 0.0 kcal/mol |
| Rev hairpin | 35.5°C (th) | 0.0 kcal/mol |
| Pair 3′ dimer | 0.0°C (th) | -0.8 kcal/mol |

### Positional agreement
Fwd pos Δ89 bp, Rev pos Δ87 bp | Fwd Tm Δ0.15°C, Rev Tm Δ0.04°C

### All returned pairs — Tm summary
| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |
|------|-----------|-----------|-----------|-----------|------------|------------|
| 1 | 60.1 | 60.0 | 60.0 | 60.0 | 411 | 261 |
| 2 | 60.1 | 60.0 | 60.0 | 60.0 | 412 | 257 |
| 3 | 60.1 | 60.1 | 59.9 | 59.9 | 410 | 420 |
| 4 | 60.1 | 60.1 | 59.8 | 59.8 | 411 | 322 |
| 5 | 60.1 | 59.9 | 60.0 | 59.9 | 451 | 419 |

---

## Summary: primd vs primer3 top-1 pair agreement

| Region | p3 pairs | primd pairs | Fwd Tm Δ (°C) | Rev Tm Δ (°C) | Fwd pos Δ (bp) | Same position? |
|--------|----------|-------------|---------------|---------------|----------------|----------------|
| lacZ-alpha | 5 | 5 | 0.29 | 0.08 | 11 | ❌ |
| bla-ampR | 5 | 5 | 0.17 | 0.24 | 0 | ❌ |
| ori-AT | 5 | 5 | 0.06 | 0.25 | 68 | ❌ |
| MCS | 5 | 5 | 0.35 | 0.13 | 10 | ❌ |
| rep-GC | 5 | 5 | 0.15 | 0.04 | 89 | ❌ |

**Position match criterion:** fwd Δ ≤ 2 bp AND rev Δ ≤ 2 bp.
