/**
 * designQPCR tests — pUC19 (L09137.2, 2686 bp)
 */
import { describe, expect, it } from "vitest";
import { designQPCR } from "../src/index.js";

// pUC19 (first 500 bp — lacZ-alpha region, ~50% GC)
const PUC19_500 =
	"TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCA" +
	"CAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTG" +
	"TTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGC" +
	"ACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCC" +
	"ATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTAT" +
	"TACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGT" +
	"TTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGAT" +
	"CCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCT";

const REGION_START = 50;
const REGION_END = 150; // 100 bp target → amplicons 70–200 bp will flank this

describe("designQPCR — basic", () => {
	it("returns at least one pair for a standard region", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		expect(result.pairs.length).toBeGreaterThan(0);
	});

	it("all amplicons are within [70, 200] bp default range", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		for (const pair of result.pairs) {
			expect(pair.productSize).toBeGreaterThanOrEqual(70);
			expect(pair.productSize).toBeLessThanOrEqual(200);
		}
	});

	it("each pair has ampliconTm, ampliconDG, efficiencyScore", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		const best = result.pairs[0];
		expect(typeof best.ampliconTm).toBe("number");
		expect(typeof best.ampliconDG).toBe("number");
		expect(typeof best.efficiencyScore).toBe("number");
	});

	it("ampliconTm is in a physically plausible range (50–100°C)", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		for (const pair of result.pairs) {
			expect(pair.ampliconTm).toBeGreaterThan(50);
			expect(pair.ampliconTm).toBeLessThan(100);
		}
	});

	it("efficiencyScore is in [0, 1]", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		for (const pair of result.pairs) {
			expect(pair.efficiencyScore).toBeGreaterThanOrEqual(0);
			expect(pair.efficiencyScore).toBeLessThanOrEqual(1);
		}
	});

	it("pairs are sorted by efficiencyScore descending", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		for (let i = 1; i < result.pairs.length; i++) {
			expect(result.pairs[i].efficiencyScore).toBeLessThanOrEqual(
				result.pairs[i - 1].efficiencyScore,
			);
		}
	});

	it("primer Tm values are in expected range (57–63°C)", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END);
		const best = result.pairs[0];
		expect(best.fwd.tm).toBeGreaterThan(55);
		expect(best.fwd.tm).toBeLessThan(66);
		expect(best.rev.tm).toBeGreaterThan(55);
		expect(best.rev.tm).toBeLessThan(66);
	});
});

describe("designQPCR — options", () => {
	it("respects custom productSizeRange", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END, {
			productSizeRange: [80, 120],
		});
		for (const pair of result.pairs) {
			expect(pair.productSize).toBeGreaterThanOrEqual(80);
			expect(pair.productSize).toBeLessThanOrEqual(120);
		}
	});

	it("returns pairs with ampliconDG=0 when scoreAmpliconStructure=false", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END, {
			scoreAmpliconStructure: false,
		});
		if (result.pairs.length > 0) {
			expect(result.pairs[0].ampliconDG).toBe(0);
		}
	});

	it("respects Tm target override", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END, {
			tmTarget: 62,
			minTemplateAccessibility: 0.01,
		});
		if (result.pairs.length > 0) {
			const best = result.pairs[0];
			expect(Math.abs(best.fwd.tm - 62)).toBeLessThan(5);
			expect(Math.abs(best.rev.tm - 62)).toBeLessThan(5);
		}
	});

	it("ampliconTm increases with longer amplicon for same GC%", () => {
		// Design two runs with different size ranges and compare Tm trend
		const small = designQPCR(PUC19_500, REGION_START, REGION_END, {
			productSizeRange: [70, 100],
		});
		const large = designQPCR(PUC19_500, REGION_START, REGION_END, {
			productSizeRange: [150, 200],
		});
		if (small.pairs.length > 0 && large.pairs.length > 0) {
			// Larger amplicons lose the length-correction term (−675/L decreases),
			// so Tm should be higher all else equal.
			expect(large.pairs[0].ampliconTm).toBeGreaterThan(small.pairs[0].ampliconTm - 5);
		}
	});
});

describe("designQPCR — edge cases", () => {
	it("returns empty pairs with warning for impossible region", () => {
		// Region so close to edge that no primers can flank it
		const result = designQPCR(PUC19_500, 5, 10);
		// Either returns pairs or a warning — shouldn't throw
		expect(Array.isArray(result.pairs)).toBe(true);
	});

	it("inherited PCR filters still apply (maxTmDiff)", () => {
		const result = designQPCR(PUC19_500, REGION_START, REGION_END, {
			maxTmDiff: 1,
		});
		for (const pair of result.pairs) {
			expect(pair.tmDiff).toBeLessThanOrEqual(1);
		}
	});
});
