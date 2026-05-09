/**
 * Salt correction models for DNA duplex melting temperature.
 *
 * Two models are provided:
 *
 * 1. Owczarzy 2004 — monovalent ions + simplified Mg²⁺ correction.
 *    Owczarzy et al. (2004) Biochemistry 43:3537–3554.
 *    This is what primer3 uses.
 *
 * 2. Owczarzy 2008 — improved Mg²⁺ model that accounts for free Mg²⁺
 *    after chelation by dNTPs, and uses a ratio-based correction for
 *    mixed mono/divalent conditions. Recommended for typical PCR conditions.
 *    Owczarzy et al. (2008) Biochemistry 47:5336–5353.
 *
 * Input Tm must be the value at 1M NaCl (output of the NN calculation).
 * All concentrations in mol/L.
 */

export type SaltModel = "owczarzy_2004" | "owczarzy_2008";

export interface SaltOptions {
	model?: SaltModel;
	/** Monovalent cation concentration (mol/L). Default: 0.05 */
	monoConc?: number;
	/** Total Mg²⁺ concentration (mol/L). Default: 0.002 */
	mgConc?: number;
	/** Total dNTP concentration (mol/L). Default: 0.0008 */
	dntpConc?: number;
	/** Oligo/primer length in bases (needed for 2008 model). */
	oligoLen?: number;
}

/**
 * Apply the Owczarzy 2004 monovalent correction.
 * Eq. 22 from the paper: 1/Tm = 1/Tm(1M) + (4.29*fGC - 3.95)*1e-5*ln[Na+] + 9.40e-6*(ln[Na+])²
 */
function applyMono2004(tm1M: number, monoConc: number, fGC: number): number {
	const lnNa = Math.log(monoConc);
	const invTm = 1 / (tm1M + 273.15) + (4.29 * fGC - 3.95) * 1e-5 * lnNa + 9.40e-6 * lnNa * lnNa;
	return 1 / invTm - 273.15;
}

/**
 * Compute free Mg²⁺ concentration after dNTP chelation.
 * dNTPs bind Mg²⁺ with ~1:1 stoichiometry at physiological pH.
 * Uses the quadratic formula from Owczarzy 2008 supplementary methods.
 */
function freeMg(totalMg: number, totalDNTP: number): number {
	if (totalDNTP === 0) return totalMg;
	// Kd for Mg²⁺–dNTP ≈ 0.57 mM (empirical)
	const Kd = 5.7e-4;
	const b = Kd + totalDNTP - totalMg;
	const discriminant = b * b + 4 * Kd * totalMg;
	return (-b + Math.sqrt(discriminant)) / 2;
}

/**
 * Apply the Owczarzy 2008 Mg²⁺ correction.
 * Uses different equations depending on the [Mg²⁺]/[mono] ratio (R):
 *   R < 0.22 → use 2004 mono-only formula
 *   0.22 ≤ R < 6.0 → use mixed ionic strength equation
 *   R ≥ 6.0 → use pure Mg²⁺ equation
 *
 * Eq. 16 from the paper; in the mixed regime (0.22 ≤ R < 6.0) the same equation
 * is used but with monovalent-dependent correction factors applied to a, d, and g.
 */
function applyMg2008(
	tm1M: number,
	monoConc: number,
	totalMg: number,
	totalDNTP: number,
	oligoLen: number,
	fGC: number,
): number {
	const mg = freeMg(totalMg, totalDNTP);

	if (mg <= 0) {
		// No free Mg²⁺ → fall back to monovalent correction
		return applyMono2004(tm1M, monoConc, fGC);
	}

	const sqrtMg = Math.sqrt(mg);
	const lnMg = Math.log(mg);
	const N = oligoLen; // primer length
	const ratio = sqrtMg / (monoConc > 0 ? monoConc : 1e-9); // [Mg²⁺]^0.5 / [mono]

	let invTm: number;

	if (ratio < 0.22) {
		// Monovalent dominates — use 2004 formula
		return applyMono2004(tm1M, monoConc, fGC);
	} else {
		// Base coefficients from Owczarzy 2008 Table 2 (Eq. 16, pure Mg²⁺ regime).
		// Source: primer3 oligotm.c; Owczarzy et al. 2008 Biochemistry 47:5336–5353.
		let a = 3.92e-5;
		const b = -9.11e-6;
		const c = 6.26e-5;
		let d = 1.42e-5;
		const e = -4.82e-4;
		const f = 5.25e-4;
		let g = 8.31e-5;

		if (ratio < 6.0) {
			// Mixed regime (0.22 ≤ R < 6.0): adjust a, d, g with monovalent-dependent
			// correction factors (Owczarzy 2008 Eq. 16 mixed-regime modification).
			// b, c, e, f are unchanged.
			const lnMono = Math.log(monoConc);
			a = 3.92e-5 * (0.843 - 0.352 * Math.sqrt(monoConc) * lnMono);
			d = 1.42e-5 * (1.279 - 4.03e-3 * lnMono - 8.03e-3 * lnMono * lnMono);
			g = 8.31e-5 * (0.486 - 0.258 * lnMono + 5.25e-3 * lnMono * lnMono * lnMono);
		}

		invTm = 1 / (tm1M + 273.15)
			+ a
			+ b * lnMg
			+ fGC * (c + d * lnMg)
			+ (1 / (2 * (N - 1))) * (e + f * lnMg + g * lnMg * lnMg);
	}

	return 1 / invTm - 273.15;
}

/**
 * Apply salt correction to a Tm calculated at 1M NaCl.
 *
 * @param tm1M   Tm at 1M NaCl in °C (from NN calculation)
 * @param fGC    GC fraction of the oligo [0, 1]
 * @param opts   Salt concentrations and model selection
 * @returns      Corrected Tm in °C
 */
export function applySaltCorrection(
	tm1M: number,
	fGC: number,
	opts: SaltOptions = {},
): number {
	const model = opts.model ?? "owczarzy_2008";
	const mono = opts.monoConc ?? 0.05;
	const mg = opts.mgConc ?? 0.002;
	const dntp = opts.dntpConc ?? 0.0008;
	const N = opts.oligoLen ?? 20;

	if (model === "owczarzy_2004") {
		// 2004: monovalent correction only (primer3 behavior)
		// Simple Mg²⁺ correction: sqrt(Mg) treated as additional mono contribution
		const effectiveMono = mono > 0 ? mono : Math.sqrt(mg) * 120; // rough equivalence
		return applyMono2004(tm1M, effectiveMono, fGC);
	}

	// 2008 model
	return applyMg2008(tm1M, mono, mg, dntp, N, fGC);
}
