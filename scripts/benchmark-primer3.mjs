#!/usr/bin/env node
/**
 * Benchmark primd designPCR against primer3 on five pUC19 regions.
 *
 * Both tools use identical conditions:
 *   - SantaLucia 1998 nearest-neighbor parameters
 *   - Owczarzy 2004 monovalent salt correction
 *   - 50 mM NaCl, 0 mM Mg²⁺ (removes divalent handling as a variable)
 *   - 250 nM oligo concentration
 *   - Tm target 60°C, primer length 18–27 bp, product 150–500 bp
 *
 * Usage: node scripts/benchmark-primer3.mjs
 * Writes results to benchmark_results/primer3_comparison.md
 */

import { execSync, spawnSync } from "node:child_process";
import { mkdirSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const __dir = dirname(fileURLToPath(import.meta.url));
const ROOT = join(__dir, "..");

// ── pUC19 sequence (L09137.2, 2686 bp) ───────────────────────────────────────
const PUC19 =
  "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCA" +
  "CAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTG" +
  "TTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGC" +
  "ACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCC" +
  "ATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTAT" +
  "TACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGT" +
  "TTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGAT" +
  "CCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCT" +
  "GTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGT" +
  "AAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCC" +
  "GCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGG" +
  "AGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCG" +
  "GTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACA" +
  "GAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAAC" +
  "CGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCAC" +
  "AAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCG" +
  "TTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATAC" +
  "CTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTAT" +
  "CTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAG" +
  "CCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGAC" +
  "TTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGT" +
  "GCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGT" +
  "ATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGC" +
  "AAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGA" +
  "AAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAAC" +
  "GAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATC" +
  "CTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCT" +
  "GACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCA" +
  "TCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCT" +
  "GGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCA" +
  "ATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCC" +
  "ATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTG" +
  "CGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCT" +
  "TCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAA" +
  "AAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTA" +
  "TCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGC" +
  "TTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCG" +
  "AGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAA" +
  "GTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTG" +
  "AGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTC" +
  "ACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGG" +
  "GCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTAT" +
  "CAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATA" +
  "GGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATC" +
  "ATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC";

// ── Benchmark regions ─────────────────────────────────────────────────────────
const REGIONS = [
  { id: "lacZ-alpha",    start: 200,  end: 400,  label: "lacZ-alpha (200–400), 50% GC" },
  { id: "bla-ampR",      start: 1800, end: 2000, label: "bla/ampR (1800–2000), 55% GC" },
  { id: "ori-AT",        start: 700,  end: 900,  label: "ori/AT-rich (700–900), 45% GC" },
  { id: "MCS",           start: 390,  end: 570,  label: "MCS region (390–570), 52% GC" },
  { id: "rep-GC",        start: 1100, end: 1300, label: "rep GC-rich (1100–1300), 62% GC" },
];

const P3_CONFIG = "/opt/homebrew/share/primer3/primer3_config/";

// ── primer3 runner ────────────────────────────────────────────────────────────
function runPrimer3(seq, regionStart, regionEnd) {
  const targetLen = regionEnd - regionStart;
  const input = [
    `SEQUENCE_ID=benchmark`,
    `SEQUENCE_TEMPLATE=${seq}`,
    `SEQUENCE_TARGET=${regionStart},${targetLen}`,
    `PRIMER_TASK=generic`,
    `PRIMER_PICK_LEFT_PRIMER=1`,
    `PRIMER_PICK_RIGHT_PRIMER=1`,
    `PRIMER_OPT_SIZE=20`,
    `PRIMER_MIN_SIZE=18`,
    `PRIMER_MAX_SIZE=27`,
    `PRIMER_OPT_TM=60.0`,
    `PRIMER_MIN_TM=57.0`,
    `PRIMER_MAX_TM=63.0`,
    `PRIMER_MIN_GC=40.0`,
    `PRIMER_MAX_GC=65.0`,
    `PRIMER_SALT_MONOVALENT=50.0`,
    `PRIMER_SALT_DIVALENT=0.0`,
    `PRIMER_DNTP_CONC=0.0`,
    `PRIMER_DNA_CONC=250.0`,
    `PRIMER_SALT_CORRECTIONS=2`,
    `PRIMER_TM_FORMULA=1`,
    `PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1`,
    `PRIMER_THERMODYNAMIC_PARAMETERS_PATH=${P3_CONFIG}`,
    `PRIMER_NUM_RETURN=5`,
    `PRIMER_PRODUCT_SIZE_RANGE=150-500`,
    `=`,
  ].join("\n");

  const result = spawnSync("primer3_core", [], {
    input,
    encoding: "utf8",
    timeout: 15000,
  });

  if (result.error) throw new Error(`primer3_core failed: ${result.error.message}`);

  const pairs = [];
  const lines = result.stdout.split("\n");
  const kv = {};
  for (const line of lines) {
    const eq = line.indexOf("=");
    if (eq > 0) kv[line.slice(0, eq)] = line.slice(eq + 1);
  }

  for (let i = 0; i < 5; i++) {
    const leftSeq  = kv[`PRIMER_LEFT_${i}_SEQUENCE`];
    const rightSeq = kv[`PRIMER_RIGHT_${i}_SEQUENCE`];
    if (!leftSeq || !rightSeq) break;

    const leftPos  = (kv[`PRIMER_LEFT_${i}`]  || "").split(",");
    const rightPos = (kv[`PRIMER_RIGHT_${i}`] || "").split(",");

    pairs.push({
      rank: i,
      fwdSeq:  leftSeq,
      revSeq:  rightSeq,
      fwdTm:   parseFloat(kv[`PRIMER_LEFT_${i}_TM`]  || "0"),
      revTm:   parseFloat(kv[`PRIMER_RIGHT_${i}_TM`] || "0"),
      fwdGC:   parseFloat(kv[`PRIMER_LEFT_${i}_GC_PERCENT`]  || "0"),  // already %
      revGC:   parseFloat(kv[`PRIMER_RIGHT_${i}_GC_PERCENT`] || "0"),
      fwdLen:  leftSeq.length,
      revLen:  rightSeq.length,
      fwdStart: parseInt(leftPos[0] || "0"),
      revStart: parseInt(rightPos[0] || "0"),
      productSize: parseInt(kv[`PRIMER_PAIR_${i}_PRODUCT_SIZE`] || "0"),
      fwdHairpin:  parseFloat(kv[`PRIMER_LEFT_${i}_HAIRPIN_TH`]  || "0"),
      revHairpin:  parseFloat(kv[`PRIMER_RIGHT_${i}_HAIRPIN_TH`] || "0"),
      pairAny:     parseFloat(kv[`PRIMER_PAIR_${i}_COMPL_ANY_TH`] || "0"),
      pairEnd:     parseFloat(kv[`PRIMER_PAIR_${i}_COMPL_END_TH`] || "0"),
    });
  }

  const warning = kv["PRIMER_ERROR"] || kv["PRIMER_WARNING"] || null;
  return { pairs, warning };
}

// ── primd runner ──────────────────────────────────────────────────────────────
async function runPrimd(seq, regionStart, regionEnd) {
  const { designPCR } = await import(join(ROOT, "dist/modes/pcr.js"));
  const result = designPCR(seq, regionStart, regionEnd, {
    saltModel:    "owczarzy_2004",
    monoConc:     0.05,
    mgConc:       0,
    dntpConc:     0,
    oligoConc:    250e-9,
    tmTarget:     60,
    primerLenRange:    [18, 27],
    productSizeRange:  [150, 500],
    gcRange:           [0.40, 0.65],
    maxTmDiff:         3,
    numReturn:         5,
  });

  const pairs = result.pairs.map((p, i) => ({
    rank: i,
    fwdSeq:  p.fwd.seq,
    revSeq:  p.rev.seq,
    fwdTm:   p.fwd.tm,
    revTm:   p.rev.tm,
    fwdGC:   p.fwd.gc,
    revGC:   p.rev.gc,
    fwdLen:  p.fwd.len,
    revLen:  p.rev.len,
    fwdStart: p.fwd.start,
    revStart: p.rev.start,
    productSize: p.productSize,
    fwdHairpin: p.fwd.hairpinDG,
    revHairpin: p.rev.hairpinDG,
    pairEnd:    p.heteroDimerDG,
    tmDiff:     p.tmDiff,
    penalty:    p.penalty,
  }));

  return { pairs, warning: result.warning };
}

// ── comparison helpers ────────────────────────────────────────────────────────
function pct(v) { return `${(v * 100).toFixed(0)}%`; }
function tm(v)  { return v.toFixed(1); }
function dg(v)  { return v.toFixed(1); }

function fmtPair(p, tool) {
  if (!p) return "—";
  const tmDiff = Math.abs(p.fwdTm - p.revTm).toFixed(1);
  return [
    `Fwd: \`${p.fwdSeq}\``,
    `Rev: \`${p.revSeq}\``,
    `Tm: ${tm(p.fwdTm)} / ${tm(p.revTm)}°C (Δ${tmDiff}°C)`,
    `GC: ${pct(p.fwdGC)} / ${pct(p.revGC)}`,
    `Len: ${p.fwdLen} / ${p.revLen} bp`,
    `Product: ${p.productSize} bp`,
    tool === "primer3"
      ? `Hairpin (°C): ${p.fwdHairpin.toFixed(1)} / ${p.revHairpin.toFixed(1)}  Pair end: ${p.pairEnd.toFixed(1)}°C`
      : `Hairpin ΔG: ${dg(p.fwdHairpin)} / ${dg(p.revHairpin)} kcal/mol  Het-dimer: ${dg(p.pairEnd)} kcal/mol`,
  ].join("  \n  ");
}

function overlap(p3, pm) {
  if (!p3 || !pm) return "—";
  const fwdDelta = Math.abs(p3.fwdStart - pm.fwdStart);
  const revDelta = Math.abs(p3.revStart - pm.revStart);
  const tmFwdDelta = Math.abs(p3.fwdTm - pm.fwdTm).toFixed(2);
  const tmRevDelta = Math.abs(p3.revTm - pm.revTm).toFixed(2);
  return `Fwd pos Δ${fwdDelta} bp, Rev pos Δ${revDelta} bp | Fwd Tm Δ${tmFwdDelta}°C, Rev Tm Δ${tmRevDelta}°C`;
}

// ── main ──────────────────────────────────────────────────────────────────────
const lines = [
  "# primd vs primer3 benchmark — pUC19 (L09137.2, 2686 bp)",
  "",
  "**Conditions (both tools):** SantaLucia 1998 NN · Owczarzy 2004 salt correction · 50 mM NaCl · 0 mM Mg²⁺ · 250 nM oligo · Tm target 60°C · primer 18–27 bp · product 150–500 bp · GC 40–65%",
  "",
  "> primer3 reports hairpin and dimer thresholds in °C (thermodynamic ΔG converted internally).  ",
  "> primd reports hairpin and heterodimer ΔG directly in kcal/mol.  ",
  "> Negative primd hairpin ΔG < −2 kcal/mol corresponds roughly to primer3 hairpin_th > 47°C.",
  "",
];

const summary = [];

for (const region of REGIONS) {
  process.stderr.write(`Running ${region.id}...\n`);

  const [p3result, pmResult] = await Promise.all([
    Promise.resolve(runPrimer3(PUC19, region.start, region.end)),
    runPrimd(PUC19, region.start, region.end),
  ]);

  const p3 = p3result.pairs[0];
  const pm = pmResult.pairs[0];
  const n  = pmResult.pairs.length;

  lines.push(`## ${region.label}`, "");
  lines.push(`**primer3** (${p3result.pairs.length} pairs returned)  `);
  lines.push(`**primd** (${n} pairs returned)`, "");

  lines.push("### Best pair comparison", "");
  lines.push("| | primer3 | primd |");
  lines.push("|-|---------|-------|");
  lines.push(`| Fwd seq    | \`${p3?.fwdSeq ?? "—"}\` | \`${pm?.fwdSeq ?? "—"}\` |`);
  lines.push(`| Rev seq    | \`${p3?.revSeq ?? "—"}\` | \`${pm?.revSeq ?? "—"}\` |`);
  lines.push(`| Fwd Tm     | ${p3 ? tm(p3.fwdTm)+"°C" : "—"} | ${pm ? tm(pm.fwdTm)+"°C" : "—"} |`);
  lines.push(`| Rev Tm     | ${p3 ? tm(p3.revTm)+"°C" : "—"} | ${pm ? tm(pm.revTm)+"°C" : "—"} |`);
  lines.push(`| Fwd Tm Δ60 | ${p3 ? (p3.fwdTm-60).toFixed(2)+"°C" : "—"} | ${pm ? (pm.fwdTm-60).toFixed(2)+"°C" : "—"} |`);
  lines.push(`| Rev Tm Δ60 | ${p3 ? (p3.revTm-60).toFixed(2)+"°C" : "—"} | ${pm ? (pm.revTm-60).toFixed(2)+"°C" : "—"} |`);
  lines.push(`| Fwd GC     | ${p3 ? p3.fwdGC.toFixed(0)+"%" : "—"} | ${pm ? pct(pm.fwdGC) : "—"} |`);
  lines.push(`| Rev GC     | ${p3 ? p3.revGC.toFixed(0)+"%" : "—"} | ${pm ? pct(pm.revGC) : "—"} |`);
  lines.push(`| Fwd len    | ${p3?.fwdLen ?? "—"} bp | ${pm?.fwdLen ?? "—"} bp |`);
  lines.push(`| Rev len    | ${p3?.revLen ?? "—"} bp | ${pm?.revLen ?? "—"} bp |`);
  lines.push(`| Product    | ${p3?.productSize ?? "—"} bp | ${pm?.productSize ?? "—"} bp |`);
  lines.push(`| Fwd start  | ${p3?.fwdStart ?? "—"} | ${pm?.fwdStart ?? "—"} |`);
  lines.push(`| Rev start  | ${p3?.revStart ?? "—"} | ${pm?.revStart ?? "—"} |`);
  lines.push(`| Fwd hairpin | ${p3 ? p3.fwdHairpin.toFixed(1)+"°C (th)" : "—"} | ${pm ? dg(pm.fwdHairpin)+" kcal/mol" : "—"} |`);
  lines.push(`| Rev hairpin | ${p3 ? p3.revHairpin.toFixed(1)+"°C (th)" : "—"} | ${pm ? dg(pm.revHairpin)+" kcal/mol" : "—"} |`);
  lines.push(`| Pair 3′ dimer | ${p3 ? p3.pairEnd.toFixed(1)+"°C (th)" : "—"} | ${pm ? dg(pm.pairEnd)+" kcal/mol" : "—"} |`);
  lines.push("");

  // Positional overlap
  lines.push("### Positional agreement");
  lines.push(overlap(p3, pm), "");

  // All top-5 Tm values
  lines.push("### All returned pairs — Tm summary");
  lines.push("| Rank | p3 Fwd Tm | p3 Rev Tm | pm Fwd Tm | pm Rev Tm | p3 product | pm product |");
  lines.push("|------|-----------|-----------|-----------|-----------|------------|------------|");
  for (let i = 0; i < Math.max(p3result.pairs.length, pmResult.pairs.length); i++) {
    const a = p3result.pairs[i];
    const b = pmResult.pairs[i];
    lines.push(`| ${i+1} | ${a ? tm(a.fwdTm) : "—"} | ${a ? tm(a.revTm) : "—"} | ${b ? tm(b.fwdTm) : "—"} | ${b ? tm(b.revTm) : "—"} | ${a?.productSize ?? "—"} | ${b?.productSize ?? "—"} |`);
  }
  lines.push("");

  if (p3result.warning) lines.push(`> primer3 warning: ${p3result.warning}`, "");
  if (pmResult.warning) lines.push(`> primd warning: ${pmResult.warning}`, "");

  // Summary row for aggregate table
  const fwdTmDelta = (p3 && pm) ? Math.abs(p3.fwdTm - pm.fwdTm).toFixed(2) : "—";
  const revTmDelta = (p3 && pm) ? Math.abs(p3.revTm - pm.revTm).toFixed(2) : "—";
  const fwdPosDelta = (p3 && pm) ? Math.abs(p3.fwdStart - pm.fwdStart) : "—";
  const samePos = (p3 && pm) ? (Math.abs(p3.fwdStart - pm.fwdStart) <= 2 && Math.abs(p3.revStart - pm.revStart) <= 2 ? "✅" : "❌") : "—";
  summary.push({ region: region.id, p3n: p3result.pairs.length, pmn: n, fwdTmDelta, revTmDelta, fwdPosDelta, samePos });
}

// ── Summary table ─────────────────────────────────────────────────────────────
lines.push("---", "");
lines.push("## Summary: primd vs primer3 top-1 pair agreement", "");
lines.push("| Region | p3 pairs | primd pairs | Fwd Tm Δ (°C) | Rev Tm Δ (°C) | Fwd pos Δ (bp) | Same position? |");
lines.push("|--------|----------|-------------|---------------|---------------|----------------|----------------|");
for (const s of summary) {
  lines.push(`| ${s.region} | ${s.p3n} | ${s.pmn} | ${s.fwdTmDelta} | ${s.revTmDelta} | ${s.fwdPosDelta} | ${s.samePos} |`);
}
lines.push("");
lines.push("**Position match criterion:** fwd Δ ≤ 2 bp AND rev Δ ≤ 2 bp.");

const outDir = join(ROOT, "benchmark_results");
mkdirSync(outDir, { recursive: true });
const outFile = join(outDir, "primer3_comparison.md");
writeFileSync(outFile, lines.join("\n") + "\n");
process.stderr.write(`\nResults written to ${outFile}\n`);
console.log(lines.join("\n"));
