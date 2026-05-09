/**
 * Nearest-neighbor thermodynamic parameters for DNA duplex stability.
 *
 * SantaLucia 1998 — the unified NN model, industry standard, used by primer3.
 * SantaLucia J Jr. (1998) PNAS 95:1460–1465.
 *
 * Units: dH in kcal/mol, dS in cal/mol/K.
 * All parameters for non-self-complementary duplexes, 5'→3' top strand.
 *
 * Note on alternative parameter sets: Ke et al. (2025) Nature Commun 16:5572
 * (DOI: 10.1038/s41467-025-60455-4) updated NUPACK hairpin loop and mismatch
 * parameters but explicitly held WC propagation stacks fixed at dna04 values
 * (which are themselves consistent with SL98). There is therefore no published
 * correction table for the WC stacking terms used in primer Tm calculation.
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

export type NNModel = "santa_lucia_1998";

/**
 * Compute raw ΔH and ΔS for a duplex using the SantaLucia 1998 NN model.
 * Includes propagation terms and initiation corrections.
 * Does NOT include salt correction (handled separately by applySaltCorrection).
 */
export function calcNNThermo(seq: string, model: NNModel = "santa_lucia_1998"): NNParam {
	void model; // single model for now; parameter kept for API forward-compatibility
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
	}

	// Initiation: one term per terminal base pair
	for (const b of [s[0], s[s.length - 1]]) {
		const init = (b === "G" || b === "C") ? SL98_INIT_GC : SL98_INIT_AT;
		dH += init.dH;
		dS += init.dS;
	}

	return { dH, dS };
}
