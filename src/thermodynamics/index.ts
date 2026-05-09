/**
 * High-level thermodynamics API.
 * Combines NN calculation + salt correction into a single calcTm() call.
 */

import { calcNNThermo, type NNModel } from "./nearest-neighbor.js";
import { applySaltCorrection, type SaltModel } from "./salt-correction.js";

export type { NNModel, SaltModel };
export { calcNNThermo } from "./nearest-neighbor.js";
export { applySaltCorrection } from "./salt-correction.js";

export interface TmOptions {
	nnModel?: NNModel;
	saltModel?: SaltModel;
	oligoConc?: number;   // mol/L, default 250e-9
	monoConc?: number;    // mol/L, default 0.05
	mgConc?: number;      // mol/L, default 0.002
	dntpConc?: number;    // mol/L, default 0.0008
}

export interface TmResult {
	tm: number;   // °C, salt-corrected
	tm1M: number; // °C, at 1M NaCl (before salt correction)
	dH: number;   // kcal/mol
	dS: number;   // cal/mol/K
	dG37: number; // kcal/mol at 37°C
}

const R = 1.987; // cal/(mol·K) — gas constant

/** Compute the full melting temperature of a DNA oligo with salt correction. */
export function calcTm(seq: string, opts: TmOptions = {}): TmResult {
	const nnModel = opts.nnModel ?? "santa_lucia_1998";
	const CT = opts.oligoConc ?? 250e-9;

	const { dH, dS } = calcNNThermo(seq, nnModel);

	// Tm at 1M NaCl (non-self-complementary: use CT/4)
	const tm1M = (dH * 1000) / (dS + R * Math.log(CT / 4)) - 273.15;

	const gc = [...seq.toUpperCase()].filter((b) => b === "G" || b === "C").length / seq.length;
	const tm = applySaltCorrection(tm1M, gc, {
		model: opts.saltModel ?? "owczarzy_2008",
		monoConc: opts.monoConc ?? 0.05,
		mgConc: opts.mgConc ?? 0.002,
		dntpConc: opts.dntpConc ?? 0.0008,
		oligoLen: seq.length,
	});

	const dG37 = dH - (37 + 273.15) * (dS / 1000);

	return { tm, tm1M, dH, dS, dG37 };
}

/** GC fraction [0,1] */
export function calcGC(seq: string): number {
	const s = seq.toUpperCase();
	return [...s].filter((b) => b === "G" || b === "C").length / s.length;
}

/** Reverse complement of a DNA sequence */
export function reverseComplement(seq: string): string {
	const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N" };
	return seq
		.toUpperCase()
		.split("")
		.reverse()
		.map((b) => comp[b] ?? "N")
		.join("");
}
