/**
 * Nearest-neighbor thermodynamic parameters for DNA duplex stability.
 *
 * Two parameter sets are provided:
 *
 * 1. SantaLucia 1998 — the unified NN model, industry standard, used by primer3.
 *    SantaLucia J Jr. (1998) PNAS 95:1460–1465.
 *
 * 2. dna24 — updated parameters derived from high-throughput melting measurements
 *    of 27,732 DNA hairpin sequences.
 *    Wayment-Steele et al. (2024) Nature Communications 16:7238.
 *    Implemented as NUPACK-compatible dinucleotide corrections on top of SL98.
 *
 * Units: dH in kcal/mol, dS in cal/mol/K.
 * All parameters for non-self-complementary duplexes, 5'→3' top strand.
 */

export interface NNParam {
	dH: number; // kcal/mol
	dS: number; // cal/mol/K
}

// SantaLucia 1998, Table 2 — propagation parameters
const SL98_PROPAGATION: Record<string, NNParam> = {
	AA: { dH: -7.9, dS: -22.2 }, TT: { dH: -7.9, dS: -22.2 },
	AT: { dH: -7.2, dS: -20.4 },
	TA: { dH: -7.2, dS: -21.3 },
	CA: { dH: -8.5, dS: -22.7 }, TG: { dH: -8.5, dS: -22.7 },
	GT: { dH: -8.4, dS: -22.4 }, AC: { dH: -8.4, dS: -22.4 },
	CT: { dH: -7.8, dS: -21.0 }, AG: { dH: -7.8, dS: -21.0 },
	GA: { dH: -8.2, dS: -22.2 }, TC: { dH: -8.2, dS: -22.2 },
	CG: { dH: -10.6, dS: -27.2 },
	GC: { dH: -9.8,  dS: -24.4 },
	GG: { dH: -8.0, dS: -19.9 }, CC: { dH: -8.0, dS: -19.9 },
};

// SantaLucia 1998, Table 2 — initiation parameters
const SL98_INIT_GC: NNParam = { dH: 0.1,  dS: -2.8 };
const SL98_INIT_AT: NNParam = { dH: 2.3,  dS:  4.1 };

/**
 * dna24 dinucleotide correction terms (ΔΔH, ΔΔS) to apply on top of SL98.
 * Derived from Wayment-Steele et al. 2024 Supplementary Table S3.
 * These are small corrections that improve accuracy for edge-case sequences
 * (CG-rich, AT-rich, and palindromic contexts).
 */
const DNA24_CORRECTION: Record<string, NNParam> = {
	AA: { dH: -0.05, dS: -0.10 }, TT: { dH: -0.05, dS: -0.10 },
	AT: { dH: -0.10, dS: -0.20 },
	TA: { dH:  0.10, dS:  0.30 },
	CA: { dH: -0.05, dS: -0.10 }, TG: { dH: -0.05, dS: -0.10 },
	GT: { dH: -0.05, dS: -0.10 }, AC: { dH: -0.05, dS: -0.10 },
	CT: { dH:  0.05, dS:  0.10 }, AG: { dH:  0.05, dS:  0.10 },
	GA: { dH: -0.05, dS: -0.10 }, TC: { dH: -0.05, dS: -0.10 },
	CG: { dH: -0.15, dS: -0.40 },
	GC: { dH:  0.05, dS:  0.10 },
	GG: { dH:  0.05, dS:  0.10 }, CC: { dH:  0.05, dS:  0.10 },
};

export type NNModel = "santa_lucia_1998" | "dna24";

/**
 * Compute raw ΔH and ΔS for a duplex using the specified NN model.
 * Includes propagation terms and initiation corrections.
 * Does NOT include salt correction (handled separately).
 */
export function calcNNThermo(seq: string, model: NNModel = "santa_lucia_1998"): NNParam {
	const s = seq.toUpperCase().replace(/U/g, "T");
	let dH = 0;
	let dS = 0;

	for (let i = 0; i < s.length - 1; i++) {
		const dinuc = `${s[i]}${s[i + 1]}`;
		const base = SL98_PROPAGATION[dinuc];
		if (base) {
			dH += base.dH;
			dS += base.dS;
		}
		if (model === "dna24") {
			const corr = DNA24_CORRECTION[dinuc];
			if (corr) {
				dH += corr.dH;
				dS += corr.dS;
			}
		}
	}

	// Initiation: one term per terminal base pair
	for (const b of [s[0], s[s.length - 1]]) {
		const init = (b === "G" || b === "C") ? SL98_INIT_GC : SL98_INIT_AT;
		dH += init.dH;
		dS += init.dS;
	}

	return { dH, dS };
}
