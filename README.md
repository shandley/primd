# primd

**Thermodynamically-accurate primer design for the browser and Node.js.**

primd is a pure-TypeScript library for designing PCR and LAMP primers with no native dependencies and no WebAssembly. It implements the SantaLucia 1998 nearest-neighbor model with Owczarzy 2008 Mg²⁺ salt correction — the same thermodynamic foundation used by Primer3 — and adds a Boltzmann-based model for template accessibility scoring. It is designed for web applications that need publication-grade Tm calculations without shipping native binaries or WASM bundles.

## Features

- **SantaLucia 1998 nearest-neighbor Tm** — full ΔH/ΔS/ΔG parameter set, not the Wallace rule
- **Owczarzy 2008 Mg²⁺ salt correction** — accounts for dNTP chelation of free Mg²⁺
- **Secondary structure screening** — hairpin ΔG, self-dimer ΔG, and heterodimer ΔG at the 3′ end
- **Template accessibility scoring** — Boltzmann two-state model; avoids primers that land in template hairpin arms
- **PCR primer design** — forward/reverse pairs ranked by a composite penalty score
- **LAMP primer design** — full six-primer sets (F3/FIP/B3/BIP + optional LoopF/LoopB) following the Notomi 2000 geometry
- **Zero runtime dependencies** — pure TypeScript, ESM, works in modern browsers and Node.js ≥ 18

## Installation

```bash
npm install primd
# pnpm
pnpm add primd
# yarn
yarn add primd
```

## Quick start

```typescript
import { designPCR } from "primd";

const template =
  "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA" +
  "CGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTCGACCAAAGGT";

// Design primers flanking bases 20–60
const result = designPCR(template, 20, 60);

if (result.pairs.length === 0) {
  console.warn(result.warning);
} else {
  const best = result.pairs[0];
  console.log("Forward:", best.fwd.seq, `Tm=${best.fwd.tm.toFixed(1)}°C`);
  console.log("Reverse:", best.rev.seq, `Tm=${best.rev.tm.toFixed(1)}°C`);
  console.log("Product:", best.productSize, "bp");
  console.log("Pair penalty:", best.penalty.toFixed(2));
}
```

With custom conditions:

```typescript
const result = designPCR(template, 20, 60, {
  productSizeRange: [150, 500],
  tmTarget: 62,
  maxTmDiff: 2,
  mgConc: 0.003,    // 3 mM Mg²⁺
  dntpConc: 0.0002, // 0.2 mM dNTPs
  numReturn: 3,
});
```

## API reference

### `designPCR(template, regionStart, regionEnd, opts?)`

Designs primer pairs that flank the target region. Primers are placed outside `[regionStart, regionEnd)` and ranked by a composite penalty (Tm deviation, GC content, secondary structure, template accessibility, pair Tm balance, and 3′ heterodimer ΔG).

```typescript
function designPCR(
  template: string,       // full template, 5'→3', top strand
  regionStart: number,    // 0-indexed start of region to amplify
  regionEnd: number,      // 0-indexed end (exclusive)
  opts?: PCROptions,
): PCRResult
```

**`PCROptions`** (all optional):

| Option | Default | Description |
|--------|---------|-------------|
| `productSizeRange` | `[100, 1000]` | Amplicon size range in bp |
| `primerLenRange` | `[18, 27]` | Primer length range in bp |
| `tmTarget` | `60` | Target Tm in °C |
| `maxTmDiff` | `3` | Max allowed \|Tm_fwd − Tm_rev\| in °C |
| `gcRange` | `[0.40, 0.65]` | GC fraction range |
| `maxPolyRun` | `4` | Max homopolymer run length |
| `numReturn` | `5` | Number of top pairs to return |
| `monoConc` | `0.05` | Monovalent salt concentration (mol/L) |
| `mgConc` | `0.002` | Mg²⁺ concentration (mol/L) |
| `dntpConc` | `0.0008` | dNTP concentration (mol/L); reduces free Mg²⁺ |
| `annealTempC` | `tmTarget − 5` | Annealing temp for accessibility scoring |

**`PCRResult`**:

```typescript
interface PCRResult {
  pairs: PrimerPair[];
  warning?: string; // present when no pairs found — explains why
}
```

When `pairs` is empty, `warning` describes which filter eliminated candidates (e.g., no compatible Tm pairs, region too short) so you can adjust parameters.

---

### `designLAMP(template, regionStart, regionEnd, opts?)`

Designs a LAMP primer set following the Notomi et al. 2000 six-region geometry: outer primers F3/B3, inner chimeric primers FIP (F1c+F2) and BIP (B1c+B2), and optional loop primers LoopF/LoopB. Inner amplicon (F2_start → B2c_end) is constrained to 120–200 bp, so `regionEnd − regionStart` must be ≤ 200 bp.

```typescript
function designLAMP(
  template: string,
  regionStart: number,
  regionEnd: number,
  opts?: LAMPOptions,
): LAMPResult
```

**`LAMPOptions`** (all optional):

| Option | Default | Description |
|--------|---------|-------------|
| `outerTmRange` | `[59, 67]` | Tm range for F3/B3 in °C |
| `innerTmRange` | `[60, 65]` | Tm range for FIP/BIP parts in °C |
| `loopTmRange` | `[62, 68]` | Tm range for loop primers in °C |
| `designLoops` | `true` | Whether to attempt loop primer design |
| `numReturn` | `3` | Number of top sets to return |
| `monoConc` / `mgConc` / `dntpConc` | same as PCR | Salt conditions |

**`LAMPResult`**:

```typescript
interface LAMPResult {
  sets: LAMPPrimerSet[];
  warning?: string;
}

interface LAMPPrimerSet {
  F3: PrimerCandidate;
  B3: PrimerCandidate;
  FIP: LAMPInnerPrimer; // { seq, part1 (F1c), part2 (F2), tm1, tm2 }
  BIP: LAMPInnerPrimer; // { seq, part1 (B1c), part2 (B2), tm1, tm2 }
  LoopF: PrimerCandidate | null;
  LoopB: PrimerCandidate | null;
  penalty: number;
}
```

---

### `calcTm(seq, opts?)`

Returns the salt-corrected melting temperature and full thermodynamic parameters.

```typescript
function calcTm(seq: string, opts?: ThermoOptions): TmResult
// TmResult: { tm, tm1M, dH, dS, dG37 }
```

`tm` is the Owczarzy 2008-corrected Tm at the specified salt concentrations. `tm1M` is the Tm at 1 M NaCl before correction.

---

### `calcAccessibilityProfile(template, primerLen?, opts?)`

Pre-computes template accessibility for every start position. Returns a `Float32Array` of length `template.length` where `profile[i]` is the probability [0, 1] that a primer of `primerLen` bases starting at position `i` binds to a single-stranded (unstructured) template at `annealTempC`.

```typescript
function calcAccessibilityProfile(
  template: string,
  primerLen?: number,          // default 20
  opts?: AccessibilityOpts,    // { annealTempC?, windowExtra? }
): Float32Array
```

Values close to 1.0 indicate a freely accessible site; values below ~0.2 suggest the primer would compete with a stable template hairpin.

---

### Other exported utilities

| Function | Description |
|----------|-------------|
| `calcGC(seq)` | GC fraction [0, 1] |
| `reverseComplement(seq)` | Reverse complement |
| `calcHairpinDG(seq)` | Most stable hairpin ΔG (kcal/mol) |
| `calcSelfDimerDG(seq)` | Most stable 3′ self-dimer ΔG |
| `calcHeteroDimerDG(a, b)` | Most stable 3′ heterodimer ΔG |
| `calcAccessibility(seq, pos, len, opts?)` | Single-site accessibility score |
| `calcNNThermo(seq)` | Raw ΔH/ΔS from SantaLucia 1998 tables |
| `applySaltCorrection(tm1M, gc, opts)` | Apply salt correction to a 1M Tm |

> **Not yet available:** `designQPCR` and `designAssembly` are defined in the type system but not implemented.

## Thermodynamics

Melting temperatures are calculated using the unified nearest-neighbor parameter set from SantaLucia (1998) *PNAS* 95:1460–1465. Salt correction uses the Owczarzy et al. (2008) *Biochemistry* 47:5336–5353 model, which explicitly accounts for the chelation of free Mg²⁺ by dNTPs — this matters when dNTP concentrations are comparable to Mg²⁺, as in standard PCR mixes.

Template accessibility is estimated with a two-state Boltzmann model: for each candidate binding site, the most stable overlapping hairpin (identified by an exhaustive O(n³) stem scan with SantaLucia 1998 stacking parameters and Turner 2004 loop penalties) defines a folding free energy ΔG at the annealing temperature. The accessible fraction is `1 / (1 + exp(−ΔG / RT))`.

## Validation

primd has been validated against primer3 (the field standard) using two complementary approaches. Full details are in [VALIDATION.md](VALIDATION.md).

**Direct Tm comparison** — 150 synthetic oligos (18–27 bp, GC 22–74%) computed by both tools on identical sequences:

| Condition | MAE (°C) | Pearson r | Bias |
|-----------|----------|-----------|------|
| 50 mM NaCl, 0 mM Mg²⁺ (Owczarzy 2004) | **0.002** | **1.000** | +0.002°C |
| 50 mM NaCl, 2 mM Mg²⁺ (Owczarzy 2008 vs 2004) | 0.150 | 0.9996 | −0.150°C |

The NN implementation is exact relative to primer3. The 0.15°C systematic bias under Mg²⁺ reflects the intentional use of Owczarzy 2008 vs primer3's Owczarzy 2004 divalent correction — Owczarzy 2008 accounts for dNTP chelation of free Mg²⁺. Users switching from primer3 with standard PCR buffer should expect primd Tm values to be ~0.15°C lower, which is within the 1–2°C experimental noise of oligonucleotide Tm measurements.

## Browser and Node.js compatibility

primd is ESM-only and has no runtime dependencies. It uses no browser APIs in its computation path, so it runs identically in:

- Modern browsers (Chrome, Firefox, Safari, Edge)
- Web Workers (recommended for long templates to keep the UI thread free)
- Node.js ≥ 18
- Deno and Bun

## License

MIT
