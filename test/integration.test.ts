/**
 * Integration tests using the real pUC19 plasmid sequence (L09137, 2686 bp).
 *
 * Fetched from NCBI:
 *   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=fasta&retmode=text
 *
 * Two biologically meaningful target regions:
 *   Region A — lacZ-alpha fragment (positions 200–400, 0-indexed)
 *   Region B — bla ampicillin-resistance gene (positions 1800–2000, 0-indexed)
 *
 * LAMP sub-region — a 120 bp window inside lacZ-alpha (280–400, 0-indexed),
 * chosen so the inner amplicon constraint (120–200 bp) is satisfiable and
 * there is ample flanking sequence for F3 and B3.
 */

import { describe, it, expect } from "vitest";
import { designPCR } from "../src/modes/pcr.js";
import { designLAMP } from "../src/modes/lamp.js";

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

// ── Sanity check ─────────────────────────────────────────────────────────────
describe("pUC19 sequence sanity", () => {
	it("is 2686 bp", () => {
		expect(PUC19.length).toBe(2686);
	});

	it("starts with the expected prefix", () => {
		expect(PUC19.startsWith("TCGCGCGTTTCGGTGATGAC")).toBe(true);
	});
});

// ── Region A: lacZ-alpha (positions 200–400, 0-indexed) ──────────────────────
// This window lies within the lacZ-alpha coding sequence.
// Product size with primers just outside: expect 150–500 bp.
const REGION_A_START = 200;
const REGION_A_END = 400;

describe("designPCR — Region A (lacZ-alpha, 200–400)", () => {
	it("returns at least one primer pair", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		expect(result.warning).toBeUndefined();
		expect(result.pairs.length).toBeGreaterThan(0);
	});

	it("best pair has product size within [150, 500] bp", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		const best = result.pairs[0];
		expect(best).toBeDefined();
		expect(best!.productSize).toBeGreaterThanOrEqual(150);
		expect(best!.productSize).toBeLessThanOrEqual(500);
	});

	it("both primers have Tm 55–70°C", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.tm).toBeGreaterThanOrEqual(55);
		expect(fwd.tm).toBeLessThanOrEqual(70);
		expect(rev.tm).toBeGreaterThanOrEqual(55);
		expect(rev.tm).toBeLessThanOrEqual(70);
	});

	it("both primers have GC content 35–65%", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.gc).toBeGreaterThanOrEqual(0.35);
		expect(fwd.gc).toBeLessThanOrEqual(0.65);
		expect(rev.gc).toBeGreaterThanOrEqual(0.35);
		expect(rev.gc).toBeLessThanOrEqual(0.65);
	});

	it("no extreme hairpins (hairpinDG > -5 kcal/mol)", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.hairpinDG).toBeGreaterThan(-5);
		expect(rev.hairpinDG).toBeGreaterThan(-5);
	});

	it("pairs are sorted by ascending penalty", () => {
		const result = designPCR(PUC19, REGION_A_START, REGION_A_END, {
			productSizeRange: [150, 500],
		});
		const penalties = result.pairs.map((p) => p.penalty);
		for (let i = 1; i < penalties.length; i++) {
			expect(penalties[i]!).toBeGreaterThanOrEqual(penalties[i - 1]!);
		}
	});
});

// ── Region B: bla ampicillin-resistance gene (positions 1800–2000, 0-indexed) ─
// This window targets the beta-lactamase coding sequence.
const REGION_B_START = 1800;
const REGION_B_END = 2000;

describe("designPCR — Region B (bla/ampR, 1800–2000)", () => {
	it("returns at least one primer pair", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		expect(result.warning).toBeUndefined();
		expect(result.pairs.length).toBeGreaterThan(0);
	});

	it("best pair has product size within [150, 500] bp", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		const best = result.pairs[0]!;
		expect(best.productSize).toBeGreaterThanOrEqual(150);
		expect(best.productSize).toBeLessThanOrEqual(500);
	});

	it("both primers have Tm 55–70°C", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.tm).toBeGreaterThanOrEqual(55);
		expect(fwd.tm).toBeLessThanOrEqual(70);
		expect(rev.tm).toBeGreaterThanOrEqual(55);
		expect(rev.tm).toBeLessThanOrEqual(70);
	});

	it("both primers have GC content 35–65%", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.gc).toBeGreaterThanOrEqual(0.35);
		expect(fwd.gc).toBeLessThanOrEqual(0.65);
		expect(rev.gc).toBeGreaterThanOrEqual(0.35);
		expect(rev.gc).toBeLessThanOrEqual(0.65);
	});

	it("no extreme hairpins (hairpinDG > -5 kcal/mol)", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		const { fwd, rev } = result.pairs[0]!;
		expect(fwd.hairpinDG).toBeGreaterThan(-5);
		expect(rev.hairpinDG).toBeGreaterThan(-5);
	});

	it("pairs are sorted by ascending penalty", () => {
		const result = designPCR(PUC19, REGION_B_START, REGION_B_END, {
			productSizeRange: [150, 500],
		});
		const penalties = result.pairs.map((p) => p.penalty);
		for (let i = 1; i < penalties.length; i++) {
			expect(penalties[i]!).toBeGreaterThanOrEqual(penalties[i - 1]!);
		}
	});
});

// ── LAMP: sub-region of lacZ-alpha (positions 280–400, 0-indexed, 120 bp) ────
// 120 bp fits within the inner amplicon constraint (120–200 bp).
// Ample upstream sequence (280 bp) for F3; ample downstream (2286 bp) for B3.
const LAMP_START = 280;
const LAMP_END = 400;

describe("designLAMP — lacZ-alpha sub-region (280–400, 120 bp)", () => {
	it("returns at least one primer set (no warning)", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		expect(result.sets.length).toBeGreaterThan(0);
		expect(result.warning).toBeUndefined();
	});

	it("inner amplicon (F2_start → B2c_end) is 120–200 bp", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		const set = result.sets[0]!;
		const innerAmp = set.BIP.part2.end - set.FIP.part2.start;
		expect(innerAmp).toBeGreaterThanOrEqual(120);
		expect(innerAmp).toBeLessThanOrEqual(200);
	});

	it("F3 starts before FIP F2 on the template", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		const set = result.sets[0]!;
		expect(set.F3.start).toBeLessThan(set.FIP.part2.start);
	});

	it("FIP F2 starts before FIP F1c on the template", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		const set = result.sets[0]!;
		// part2 = F2 (upstream on template); part1 = F1c (downstream on template)
		expect(set.FIP.part2.start).toBeLessThan(set.FIP.part1.start);
	});

	it("BIP B1c region precedes BIP B2c on the template", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		const set = result.sets[0]!;
		// part1 = B1c (upstream / closer to target), part2 = B2 from B2c (further)
		expect(set.BIP.part1.start).toBeLessThan(set.BIP.part2.start);
	});

	it("B3 starts at or after BIP B2c end on the template", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END);
		const set = result.sets[0]!;
		expect(set.B3.start).toBeGreaterThanOrEqual(set.BIP.part2.end);
	});

	it("sets are sorted by ascending penalty", () => {
		const result = designLAMP(PUC19, LAMP_START, LAMP_END, { numReturn: 3 });
		const penalties = result.sets.map((s) => s.penalty);
		for (let i = 1; i < penalties.length; i++) {
			expect(penalties[i]!).toBeGreaterThanOrEqual(penalties[i - 1]!);
		}
	});
});
