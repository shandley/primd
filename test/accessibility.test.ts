/**
 * Template accessibility tests.
 *
 * Physical laws tested:
 * 1. A site with no possible WC hairpin returns 1.0.
 * 2. A site whose binding region is in a GC-rich hairpin stem at the annealing
 *    temperature has significantly lower accessibility than at the same site
 *    without the stem (i.e., high ΔG → high accessibility, low ΔG → low).
 * 3. Accessibility increases with temperature (hotter = more unfolded).
 * 4. A site far from any hairpin in the window is fully accessible.
 * 5. The accessibility profile covers all positions.
 */

import { describe, it, expect } from "vitest";
import { calcAccessibility, calcAccessibilityProfile } from "../src/thermodynamics/accessibility.js";

// ── Helpers ───────────────────────────────────────────────────────────────────

// Build a template where positions [primerStart, primerStart+primerLen) are in
// a GC-rich hairpin stem. The stem bridges [left, left+stemLen) and [right, right+stemLen).
// Flanking is all-A to avoid incidental structure.
function makeHairpinTemplate(
	leftStem: string,
	loop: string,
	primerStart: number,
	primerLen: number,
): { template: string; primerStart: number; primerLen: number } {
	// revcomp of leftStem = rightStem (antiparallel complement)
	const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };
	const rightStem = leftStem.split("").reverse().map((b) => comp[b] ?? "N").join("");

	// Layout: [flankL][leftStem][loop][rightStem][flankR]
	// We want the primer site to cover [leftStemStart, leftStemStart+primerLen)
	const flankLen = 10;
	const prefix = "A".repeat(flankLen);
	const template = prefix + leftStem + loop + rightStem + "A".repeat(flankLen);

	return {
		template,
		primerStart: flankLen, // primer starts at the beginning of the left stem
		primerLen,
	};
}

// ── Tests ─────────────────────────────────────────────────────────────────────

describe("calcAccessibility — physical correctness", () => {
	it("returns 1.0 for a poly-A site (no hairpin possible)", () => {
		const template = "A".repeat(80);
		const score = calcAccessibility(template, 30, 20);
		expect(score).toBe(1.0);
	});

	it("returns 1.0 for a primer site far from any stem in the window", () => {
		// GC hairpin far away (position 0–19), primer site at position 60+
		// The hairpin is outside the default ±25 bp window around position 60
		const hairpin = "GCGCGCAAAAAGCGCGC"; // length 17
		const filler = "A".repeat(50);
		const template = hairpin + filler + "ATGATGATGATGATGATGATG"; // primer site at 67+
		const primerStart = hairpin.length + filler.length;
		const score = calcAccessibility(template, primerStart, 20);
		expect(score).toBe(1.0);
	});

	it("GC-rich hairpin overlapping primer site reduces accessibility below 0.5", () => {
		// 6-bp GC stem: GCGCGC | AAAAA (loop) | GCGCGC
		// Primer covers the left stem (positions 10–19, length 6+)
		const leftStem = "GCGCGC"; // 6 bp
		const loop = "AAAAA"; // 5-nt loop (valid, ≥3)
		const { template, primerStart, primerLen } = makeHairpinTemplate(
			leftStem, loop, 10, leftStem.length,
		);
		const score = calcAccessibility(template, primerStart, primerLen, { annealTempC: 55 });
		// A 6-bp GC hairpin at 55°C should have ΔG well below -1 kcal/mol → accessibility < 0.5
		expect(score).toBeLessThan(0.5);
	});

	it("AT-rich stem has higher accessibility than GC-rich stem of same length", () => {
		const loop = "AAAAA";

		const gcTemplate = makeHairpinTemplate("GCGCGC", loop, 10, 6);
		const atTemplate = makeHairpinTemplate("ATATAT", loop, 10, 6);

		const gcScore = calcAccessibility(gcTemplate.template, gcTemplate.primerStart, gcTemplate.primerLen, { annealTempC: 55 });
		const atScore = calcAccessibility(atTemplate.template, atTemplate.primerStart, atTemplate.primerLen, { annealTempC: 55 });

		// GC stem is more stable → lower accessibility
		expect(gcScore).toBeLessThan(atScore);
	});

	it("accessibility increases with annealing temperature (higher T → more unfolded)", () => {
		const { template, primerStart, primerLen } = makeHairpinTemplate(
			"GCGCGC", "AAAAA", 10, 6,
		);

		const score55 = calcAccessibility(template, primerStart, primerLen, { annealTempC: 55 });
		const score72 = calcAccessibility(template, primerStart, primerLen, { annealTempC: 72 });

		expect(score72).toBeGreaterThan(score55);
	});

	it("primer site in the hairpin loop (not stem) remains accessible", () => {
		// The loop is always single-stranded — only stem arms block primers.
		// Primer covers ONLY the loop region (positions after the left stem).
		const leftStem = "GCGCGCGC"; // 8 bp
		const loop = "AAAAAAAAAAAAAAAAAAAAAAAA"; // 24-nt loop, primer fits entirely within it
		const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };
		const rightStem = leftStem.split("").reverse().map((b) => comp[b] ?? "N").join("");
		const flankLen = 10;
		const template = "A".repeat(flankLen) + leftStem + loop + rightStem + "A".repeat(flankLen);

		// Primer site starts 2 bases into the loop (not in any stem)
		const primerStart = flankLen + leftStem.length + 2; // well inside the loop
		const primerLen = 18; // entirely within the 24-nt loop

		const score = calcAccessibility(template, primerStart, primerLen, { annealTempC: 55 });
		// Loop bases are unpaired — no stem-overlap → should be 1.0
		expect(score).toBe(1.0);
	});
});

describe("calcAccessibilityProfile", () => {
	it("returns a Float32Array of the same length as the template", () => {
		const template = "ATGATGATGATGATGATGATGATGATGATG";
		const profile = calcAccessibilityProfile(template, 20);
		expect(profile).toBeInstanceOf(Float32Array);
		expect(profile.length).toBe(template.length);
	});

	it("all values are in [0, 1]", () => {
		const template = "GCGCGCAAAAAGCGCGCATGATGATGATGATGATGATGATG";
		const profile = calcAccessibilityProfile(template, 10);
		for (const v of profile) {
			expect(v).toBeGreaterThanOrEqual(0);
			expect(v).toBeLessThanOrEqual(1);
		}
	});

	it("hairpin region has lower mean accessibility than poly-A region", () => {
		// First half: GC hairpin; second half: poly-A
		const hairpinRegion = "GCGCGCGCAAAAAGCGCGCGC"; // 21 bp with GC stem
		const polyARegion = "A".repeat(30);
		const template = hairpinRegion + polyARegion;

		const profile = calcAccessibilityProfile(template, 10, { annealTempC: 55 });

		const meanHairpin = [...profile.slice(0, hairpinRegion.length)].reduce((a, b) => a + b, 0) / hairpinRegion.length;
		const meanPolyA = [...profile.slice(hairpinRegion.length)].reduce((a, b) => a + b, 0) / polyARegion.length;

		expect(meanPolyA).toBeGreaterThan(meanHairpin);
	});
});
