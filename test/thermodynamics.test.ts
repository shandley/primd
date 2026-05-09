/**
 * Thermodynamics validation tests.
 *
 * Reference values from:
 * - SantaLucia 1998 PNAS 95:1460 (Table 2, propagation parameters)
 * - Owczarzy 2008 Biochemistry 47:5336 (Mg²⁺ correction behavior)
 * - Manual calculation against SL98 equations for regression
 */

import { describe, it, expect } from "vitest";
import { calcTm, calcNNThermo } from "../src/thermodynamics/index.js";
import { calcHairpinDG, calcSelfDimerDG } from "../src/thermodynamics/secondary-structure.js";

// ── SantaLucia 1998 NN parameter sums ────────────────────────────────────────
// Validate that our NN lookup matches manual sums from Table 2.
describe("SantaLucia 1998 NN parameters", () => {
	it("correctly sums ΔH for GCATGC", () => {
		const { dH } = calcNNThermo("GCATGC");
		// Propagation: GC(-9.8)+CA(-8.5)+AT(-7.2)+TG(-8.5)+GC(-9.8) = -43.8
		// Initiation:  GC(+0.1)+GC(+0.1) = +0.2 → total = -43.6 ≈ -43.7
		expect(dH).toBeCloseTo(-43.6, 0);
	});

	it("correctly sums ΔS for GCATGC", () => {
		const { dS } = calcNNThermo("GCATGC");
		// Propagation: GC(-24.4)+CA(-22.7)+AT(-20.4)+TG(-22.7)+GC(-24.4) = -114.6
		// Initiation:  GC(-2.8)+GC(-2.8) = -5.6 → total = -120.2
		expect(dS).toBeCloseTo(-120.2, 0);
	});

	it("computes Tm below 10°C for GCATGC at 250 nM (too short to be a primer)", () => {
		const result = calcTm("GCATGC", {
			oligoConc: 250e-9,
			monoConc: 0.05,
			mgConc: 0,
			dntpConc: 0,
			saltModel: "owczarzy_2004",
		});
		// 6-mer at 250 nM has very low Tm (~2°C) — validates that concentration matters
		expect(result.tm).toBeLessThan(10);
	});
});

// ── Self-consistency checks ───────────────────────────────────────────────────
// These test that ΔG37 = ΔH - 310.15 × (ΔS/1000) and that longer/GC-richer
// primers have higher Tm — not exact values, but physical laws.
describe("Physical consistency", () => {
	it("dG37 is consistent with dH and dS", () => {
		const result = calcTm("GCTTCGGCACCAGACATGAT");
		const expectedDG = result.dH - 310.15 * (result.dS / 1000);
		expect(result.dG37).toBeCloseTo(expectedDG, 4);
	});

	it("GC-rich primer has higher Tm than AT-rich primer of same length", () => {
		const gcRich = calcTm("GCGCGCGCGCGCGCGCGCGC", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		const atRich = calcTm("AAATTTAAATTTAAATTTAA", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		expect(gcRich.tm).toBeGreaterThan(atRich.tm);
	});

	it("longer primer has higher Tm than shorter one of similar GC", () => {
		const short = calcTm("GCTTCGGCAC", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		const long = calcTm("GCTTCGGCACCAGACATGAT", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		expect(long.tm).toBeGreaterThan(short.tm);
	});

	it("AT-rich 20-mer has Tm below 45°C at standard conditions", () => {
		const result = calcTm("AAATTTAAATTTAAATTTAA", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		expect(result.tm).toBeLessThan(45);
		expect(result.tm).toBeGreaterThan(15);
	});

	it("GC-rich 20-mer has Tm above 70°C at standard conditions", () => {
		const result = calcTm("GCGCGCGCGCGCGCGCGCGC", { saltModel: "owczarzy_2004", monoConc: 0.05, mgConc: 0, dntpConc: 0 });
		expect(result.tm).toBeGreaterThan(70);
	});
});

// ── Owczarzy 2008 Mg²⁺ model behavior ────────────────────────────────────────
describe("Owczarzy 2008 salt correction", () => {
	const seq = "GCTTCGGCACCAGACATGAT";

	it("higher Mg²⁺ increases Tm", () => {
		const low = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0.0005, dntpConc: 0, monoConc: 0.05 });
		const high = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0.005, dntpConc: 0, monoConc: 0.05 });
		expect(high.tm).toBeGreaterThan(low.tm);
	});

	it("dNTPs chelate Mg²⁺ and lower Tm", () => {
		const noDNTP = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0.002, dntpConc: 0, monoConc: 0.05 });
		const withDNTP = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0.002, dntpConc: 0.0008, monoConc: 0.05 });
		expect(withDNTP.tm).toBeLessThan(noDNTP.tm);
	});

	it("Mg²⁺ increases Tm compared to monovalent-only conditions", () => {
		// Monovalent only (Owczarzy 2004 proxy)
		const monoOnly = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0, dntpConc: 0, monoConc: 0.05 });
		// With typical PCR Mg²⁺ (2mM total, 0.8mM dNTP → ~1.2mM free)
		const withMg = calcTm(seq, { saltModel: "owczarzy_2008", mgConc: 0.002, dntpConc: 0.0008, monoConc: 0.05 });
		// Free Mg²⁺ stabilizes duplexes — expect higher Tm
		expect(withMg.tm).toBeGreaterThan(monoOnly.tm);
	});

	it("Tm at typical PCR conditions is in a physiologically reasonable range (50–75°C)", () => {
		const result = calcTm(seq, { saltModel: "owczarzy_2008", oligoConc: 250e-9, monoConc: 0.05, mgConc: 0.002, dntpConc: 0.0008 });
		expect(result.tm).toBeGreaterThan(50);
		expect(result.tm).toBeLessThan(75);
	});

	it("high-Mg regime (R ≥ 6) returns finite Tm in a reasonable range", () => {
		// R = sqrt(mg) / mono = sqrt(0.01) / 0.002 = 0.1 / 0.002 = 50 ≥ 6 → Mg-dominant branch
		const result = calcTm(seq, { saltModel: "owczarzy_2008", oligoConc: 250e-9, monoConc: 0.002, mgConc: 0.01, dntpConc: 0 });
		expect(Number.isFinite(result.tm)).toBe(true);
		expect(result.tm).toBeGreaterThan(55);
		expect(result.tm).toBeLessThan(80);
	});

	it("Tm increases monotonically with Mg²⁺ in the high-Mg regime (R ≥ 6)", () => {
		// All three conditions have R ≥ 6 (monoConc=0.002, mg from 0.005 to 0.025)
		// sqrt(0.005)/0.002 = 35.4, sqrt(0.015)/0.002 = 61.2, sqrt(0.025)/0.002 = 79.1
		const low  = calcTm(seq, { saltModel: "owczarzy_2008", oligoConc: 250e-9, monoConc: 0.002, mgConc: 0.005, dntpConc: 0 });
		const mid  = calcTm(seq, { saltModel: "owczarzy_2008", oligoConc: 250e-9, monoConc: 0.002, mgConc: 0.015, dntpConc: 0 });
		const high = calcTm(seq, { saltModel: "owczarzy_2008", oligoConc: 250e-9, monoConc: 0.002, mgConc: 0.025, dntpConc: 0 });
		expect(mid.tm).toBeGreaterThan(low.tm);
		expect(high.tm).toBeGreaterThan(mid.tm);
	});
});

// ── Secondary structure ───────────────────────────────────────────────────────
describe("Secondary structure prediction", () => {
	it("stable 6-bp GC hairpin has negative ΔG", () => {
		// GCGCGC pairs with GCGCGC across a 4-nt AAAA loop
		// Stem ΔG ≈ -9 kcal/mol; loop penalty +4.1 → net ≈ -5 kcal/mol
		const dG = calcHairpinDG("GCGCGCAAAAGCGCGC");
		expect(dG).toBeLessThan(0);
	});

	it("short AT-rich sequence has hairpin ΔG near 0 (no stable structure)", () => {
		const dG = calcHairpinDG("AAACCTTTGGG");
		expect(dG).toBeGreaterThan(-3);
	});

	it("3′ self-complementary sequence has negative dimer ΔG", () => {
		const dG = calcSelfDimerDG("AAAAAAAAAAAAAAGCGC");
		expect(dG).toBeLessThan(0);
	});

	it("non-self-complementary 3′ end has dimer ΔG near 0", () => {
		const dG = calcSelfDimerDG("GCGCGCGCAATCGATCGA");
		expect(dG).toBeGreaterThan(-4);
	});
});
