/**
 * Template accessibility scoring for primer binding sites.
 *
 * During PCR annealing, the denatured single-stranded template can refold into
 * hairpin structures. If the primer binding site lies within a stem arm of such
 * a structure, the primer is sterically blocked regardless of its own Tm.
 *
 * This module estimates binding-site accessibility using a two-state Boltzmann
 * model: the site is either single-stranded (accessible) or in a hairpin stem
 * (inaccessible). The equilibrium probability of the accessible state is:
 *
 *   P_accessible = 1 / (1 + exp(−ΔG_best / RT_anneal))
 *
 * where ΔG_best is the most stable hairpin whose stem overlaps the primer site,
 * evaluated at the annealing temperature. If no overlapping hairpin is found,
 * the site is fully accessible (P = 1.0).
 *
 * ΔG is computed from SantaLucia 1998 NN stacking parameters (temperature-
 * corrected via ΔG(T) = ΔH − T·ΔS) and Turner 2004 loop penalties scaled to T.
 *
 * Limitation: the scan runs on the top strand for both forward and reverse
 * primers. For the forward primer, the bottom strand is the true binding target;
 * for most sequences the hairpin propensity is similar on both strands, but
 * highly asymmetric sequences may show small discrepancies.
 */

import { calcNNThermo } from "./nearest-neighbor.js";

const R_KCAL = 1.987e-3; // gas constant, kcal/(mol·K)

// Turner 2004 tabulated loop initiation penalties at 37 °C (kcal/mol)
const LOOP_DG_37: Record<number, number> = {
	3: 5.4, 4: 4.1, 5: 4.4, 6: 4.7, 7: 5.0, 8: 5.1, 9: 5.2,
};

/**
 * Hairpin loop penalty at temperature T (Kelvin).
 * Small loops: scale tabulated ΔG°(37°C) by T/310.15 (entropic approximation).
 * Large loops: use the ViennaRNA extrapolation formula, already linear in T.
 */
function loopPenaltyAtT(loopSize: number, T: number): number {
	if (loopSize < 3) return Infinity;
	if (loopSize <= 9) {
		return (LOOP_DG_37[loopSize] ?? 5.2) * (T / 310.15);
	}
	// Entropy-dominated large loops: ΔG ∝ T·ln(size)
	return 5.2 * (T / 310.15) + 1.75 * R_KCAL * T * Math.log(loopSize / 9);
}

const WC: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };

function isWC(a: string, b: string): boolean {
	return WC[a] === b;
}

/**
 * Find the most stable hairpin ΔG (kcal/mol) at temperature T whose
 * stem arms overlap the site [siteStart, siteEnd) within the window string.
 *
 * Only stem-overlap is counted as "blocking"; primer bases that fall in the
 * loop are unpaired and do not impede binding.
 *
 * Returns 0 if no thermodynamically stable overlapping stem is found.
 */
function bestOverlappingHairpinDG(
	win: string,
	siteStart: number,
	siteEnd: number,
	T: number,
): number {
	const n = win.length;
	let bestDG = 0;

	// Exhaustive O(n³) scan — identical to calcHairpinDG but restricted to
	// structures whose left arm [i, i+stemLen) or right arm [j, j+stemLen)
	// overlaps the primer site [siteStart, siteEnd).
	for (let i = 0; i < n - 6; i++) {
		for (let stemLen = 3; i + stemLen * 2 + 3 <= n; stemLen++) {
			const loopStart = i + stemLen;
			for (let j = loopStart + 3; j + stemLen <= n; j++) {
				// Left arm: [i, i+stemLen)  Right arm: [j, j+stemLen)
				const overlapLeft = siteStart < i + stemLen && siteEnd > i;
				const overlapRight = siteStart < j + stemLen && siteEnd > j;
				if (!overlapLeft && !overlapRight) continue;

				// Check WC complementarity: s[i+k] pairs with s[j+stemLen-1-k]
				let ok = true;
				for (let k = 0; k < stemLen; k++) {
					if (!isWC(win[i + k]!, win[j + stemLen - 1 - k]!)) {
						ok = false;
						break;
					}
				}
				if (!ok) continue;

				const loopSize = j - loopStart;
				const { dH, dS } = calcNNThermo(win.slice(i, i + stemLen));
				const stemDG = dH - T * (dS / 1000); // ΔG at annealing T, not 37°C
				const totalDG = stemDG + loopPenaltyAtT(loopSize, T);
				if (totalDG < bestDG) bestDG = totalDG;
			}
		}
	}

	return bestDG;
}

export interface AccessibilityOpts {
	/** Annealing temperature in °C. Default: 55 */
	annealTempC?: number;
	/** Extra flanking bases included in the hairpin scan window. Default: 25 */
	windowExtra?: number;
}

/**
 * Estimate the probability that a primer binding site is single-stranded
 * (accessible to primer annealing) at the annealing temperature.
 *
 * @param template    Full template sequence (5'→3', top strand)
 * @param primerStart 0-indexed start of the primer binding site on the template
 * @param primerLen   Length of the primer in bases
 * @param opts
 * @returns Accessibility ∈ [0, 1].
 *          1.0 = no stable structure overlaps the site (ideal).
 *          ~0.5 = ΔG_best ≈ 0 at anneal T (marginal).
 *          < 0.2 = site is likely structured; expect reduced efficiency.
 */
export function calcAccessibility(
	template: string,
	primerStart: number,
	primerLen: number,
	opts: AccessibilityOpts = {},
): number {
	const T = (opts.annealTempC ?? 55) + 273.15; // °C → K
	const extra = opts.windowExtra ?? 25;
	const RT = R_KCAL * T;

	const winStart = Math.max(0, primerStart - extra);
	const winEnd = Math.min(template.length, primerStart + primerLen + extra);
	const win = template.slice(winStart, winEnd).toUpperCase();

	const siteStart = primerStart - winStart;
	const siteEnd = siteStart + primerLen;

	const bestDG = bestOverlappingHairpinDG(win, siteStart, siteEnd, T);

	if (bestDG >= 0) return 1.0;

	// Two-state equilibrium: folded ⇌ unfolded
	// K_fold = exp(−ΔG_fold / RT);  ΔG_fold = bestDG < 0  →  K_fold > 1
	// P_unfolded = 1 / (1 + K_fold)
	const kFold = Math.exp(-bestDG / RT);
	return 1.0 / (1.0 + kFold);
}

/**
 * Pre-compute template accessibility for every primer start position.
 * Call this once per template in designPCR rather than per candidate.
 *
 * @param template    Full template sequence
 * @param primerLen   Representative primer length for the scan (default 20)
 * @param opts
 * @returns Float32Array of length template.length.
 *          profile[i] = accessibility score for a primer starting at position i.
 */
export function calcAccessibilityProfile(
	template: string,
	primerLen = 20,
	opts: AccessibilityOpts = {},
): Float32Array {
	const n = template.length;
	const profile = new Float32Array(n).fill(1.0);
	for (let i = 0; i <= n - primerLen; i++) {
		profile[i] = calcAccessibility(template, i, primerLen, opts);
	}
	// Trailing positions — use last valid score
	for (let i = n - primerLen + 1; i < n; i++) {
		profile[i] = profile[n - primerLen] ?? 1.0;
	}
	return profile;
}
