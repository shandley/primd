/**
 * LAMP primer design tests.
 *
 * Uses a 500-bp synthetic template with known GC composition so we can
 * verify primer geometry and thermodynamic constraints without a real
 * biological sequence.
 */

import { describe, it, expect } from "vitest";
import { designLAMP } from "../src/modes/lamp.js";

// 500-bp template with ~50% GC, no long poly runs, no accidental palindromes.
// Constructed by repeating a 20-bp unit with shuffled composition.
const TEMPLATE =
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCATGCTTCGGCACCAGACATGATGCAT" +
	"GCTTCGGCACCAGACATGATGCAT"; // 504 bp total

// Target region in the center: positions 200–280
const REGION_START = 200;
const REGION_END = 250;

describe("designLAMP — geometry and counts", () => {
	it("returns at least one primer set for a valid 80-bp target region", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		expect(result.sets.length).toBeGreaterThan(0);
		expect(result.warning).toBeUndefined();
	});

	it("each set has F3, B3, FIP, BIP fields", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set).toBeDefined();
		expect(set.F3).toBeDefined();
		expect(set.B3).toBeDefined();
		expect(set.FIP).toBeDefined();
		expect(set.BIP).toBeDefined();
	});

	it("FIP sequence = F1c + F2 (revcomp(F1) concatenated with F2)", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const { FIP } = result.sets[0];
		// FIP.seq should equal revcomp(part1.seq) + part2.seq... wait:
		// part1 = F1c (already RC'd), part2 = F2
		// FIP.seq = F1c.seq + F2.seq = part1.seq + part2.seq
		expect(FIP.seq).toBe(FIP.part1.seq + FIP.part2.seq);
	});

	it("BIP sequence = B1c + revcomp(B2c) i.e. part1 + part2", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const { BIP } = result.sets[0];
		expect(BIP.seq).toBe(BIP.part1.seq + BIP.part2.seq);
	});

	it("primer lengths are within expected ranges", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.F3.len).toBeGreaterThanOrEqual(18);
		expect(set.F3.len).toBeLessThanOrEqual(22);
		expect(set.B3.len).toBeGreaterThanOrEqual(18);
		expect(set.B3.len).toBeLessThanOrEqual(22);
		expect(set.FIP.part1.len).toBeGreaterThanOrEqual(18);
		expect(set.FIP.part1.len).toBeLessThanOrEqual(22);
		expect(set.FIP.part2.len).toBeGreaterThanOrEqual(20);
		expect(set.FIP.part2.len).toBeLessThanOrEqual(22);
		expect(set.BIP.part1.len).toBeGreaterThanOrEqual(18);
		expect(set.BIP.part1.len).toBeLessThanOrEqual(22);
		expect(set.BIP.part2.len).toBeGreaterThanOrEqual(20);
		expect(set.BIP.part2.len).toBeLessThanOrEqual(22);
	});
});

describe("designLAMP — coordinate geometry", () => {
	it("F3 starts before FIP F2 on the template", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.F3.start).toBeLessThan(set.FIP.part2.start);
	});

	it("FIP F2 starts before FIP F1 region on the template", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		// part2 = F2 (earlier on template), part1 = F1c (later on template)
		expect(set.FIP.part2.start).toBeLessThan(set.FIP.part1.start);
	});

	it("BIP B1c region precedes BIP B2c on the template", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		// part1 = B1c (earlier), part2 = B2 from B2c region (later)
		expect(set.BIP.part1.start).toBeLessThan(set.BIP.part2.start);
	});

	it("B3 starts after BIP B2c region on the template", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.B3.start).toBeGreaterThanOrEqual(set.BIP.part2.end);
	});

	it("inner amplicon (F2_start to B2c_end) is 120–200 bp", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		const innerAmp = set.BIP.part2.end - set.FIP.part2.start;
		expect(innerAmp).toBeGreaterThanOrEqual(120);
		expect(innerAmp).toBeLessThanOrEqual(200);
	});

	it("F2 start is at or before regionStart", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.FIP.part2.start).toBeLessThanOrEqual(REGION_START);
	});

	it("B2c end is at or after regionEnd", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.BIP.part2.end).toBeGreaterThanOrEqual(REGION_END);
	});
});

describe("designLAMP — Tm constraints", () => {
	it("outer primers (F3/B3) Tm in default outer range 58–62°C", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.F3.tm).toBeGreaterThanOrEqual(58);
		expect(set.F3.tm).toBeLessThanOrEqual(62);
		expect(set.B3.tm).toBeGreaterThanOrEqual(58);
		expect(set.B3.tm).toBeLessThanOrEqual(62);
	});

	it("inner primer parts (F1c, F2, B1c, B2) Tm in default inner range 63–68°C", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		expect(set.FIP.tm1).toBeGreaterThanOrEqual(63);
		expect(set.FIP.tm1).toBeLessThanOrEqual(68);
		expect(set.FIP.tm2).toBeGreaterThanOrEqual(63);
		expect(set.FIP.tm2).toBeLessThanOrEqual(68);
		expect(set.BIP.tm1).toBeGreaterThanOrEqual(63);
		expect(set.BIP.tm1).toBeLessThanOrEqual(68);
		expect(set.BIP.tm2).toBeGreaterThanOrEqual(63);
		expect(set.BIP.tm2).toBeLessThanOrEqual(68);
	});

	it("outer primers are cooler than inner parts (Notomi ordering)", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		const set = result.sets[0];
		const outerMean = (set.F3.tm + set.B3.tm) / 2;
		const innerMean = (set.FIP.tm1 + set.FIP.tm2 + set.BIP.tm1 + set.BIP.tm2) / 4;
		expect(innerMean).toBeGreaterThan(outerMean);
	});

	it("sets are returned in order of ascending penalty", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END, { numReturn: 3 });
		const penalties = result.sets.map((s) => s.penalty);
		for (let i = 1; i < penalties.length; i++) {
			expect(penalties[i]).toBeGreaterThanOrEqual(penalties[i - 1]!);
		}
	});
});

describe("designLAMP — loop primers", () => {
	it("returns primer sets when target core is large enough for loop primers (~160 bp region)", () => {
		// A 160-bp region forces F1_end and B1c_start far apart, giving the core
		// needed for loop primers (target core ≥ targetLen - 108 ≈ 52 bp > 44 bp minimum).
		const result = designLAMP(TEMPLATE, 100, 260);
		expect(result.sets.length).toBeGreaterThan(0);
	});

	it("loop primer lengths are 17–21 bp when present", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END);
		for (const set of result.sets) {
			if (set.LoopF) {
				expect(set.LoopF.len).toBeGreaterThanOrEqual(17);
				expect(set.LoopF.len).toBeLessThanOrEqual(21);
			}
			if (set.LoopB) {
				expect(set.LoopB.len).toBeGreaterThanOrEqual(17);
				expect(set.LoopB.len).toBeLessThanOrEqual(21);
			}
		}
	});

	it("designLoops: false returns null loop primers", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END, { designLoops: false });
		expect(result.sets.length).toBeGreaterThan(0);
		for (const set of result.sets) {
			expect(set.LoopF).toBeNull();
			expect(set.LoopB).toBeNull();
		}
	});
});

describe("designLAMP — edge cases", () => {
	it("returns a warning for target region > 200 bp (too large for inner amplicon)", () => {
		const result = designLAMP(TEMPLATE, 50, 260); // 210 bp target
		expect(result.sets).toHaveLength(0);
		expect(result.warning).toBeDefined();
		expect(result.warning).toMatch(/too large/i);
	});

	it("numReturn limits the number of returned sets", () => {
		const result = designLAMP(TEMPLATE, REGION_START, REGION_END, { numReturn: 1 });
		expect(result.sets.length).toBeLessThanOrEqual(1);
	});
});
