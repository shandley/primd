/**
 * Secondary structure prediction for short oligonucleotides (<60 bp).
 *
 * Uses a simplified thermodynamic folding model based on Turner 2004
 * nearest-neighbor parameters for hairpin loops and internal stack stability.
 * This is substantially more accurate than primer3's count-based approach
 * for predicting whether a primer will form self-inhibiting structures.
 *
 * For the purposes of primer design (oligos 18–35 bp), the dominant
 * secondary structures are:
 *   1. Hairpin loops (intramolecular)
 *   2. Self-dimers (intermolecular, especially 3′ end complementarity)
 *   3. Primer-pair dimers (heterodimers, 3′ end complementarity)
 *
 * Full RNA-folding algorithms (Zuker mfold, ViennaRNA) are O(n³) and
 * unnecessary for oligos this short. We use a sliding stem-loop scan.
 *
 * ΔG values in kcal/mol. Negative = thermodynamically favorable (bad for primers).
 */

import { calcNNThermo } from "./nearest-neighbor.js";

// Turner 2004 hairpin loop initiation penalties (kcal/mol)
// For loops of size 3–9+; size 3 (triloop) is special-cased
const HAIRPIN_LOOP_PENALTY: Record<number, number> = {
	3: 5.4,
	4: 4.1,
	5: 4.4,
	6: 4.7,
	7: 5.0,
	8: 5.1,
	9: 5.2,
};

function hairpinLoopPenalty(loopSize: number): number {
	if (loopSize < 3) return Infinity;
	if (loopSize <= 9) return HAIRPIN_LOOP_PENALTY[loopSize] ?? 5.2;
	// For larger loops: 5.2 + 1.75 * R * T * ln(loopSize/9)
	return 5.2 + 1.75 * 1.987e-3 * 310.15 * Math.log(loopSize / 9);
}

const COMP: Record<string, string> = { A: "T", T: "A", G: "C", C: "G" };

function complement(b: string): string {
	return COMP[b] ?? "N";
}

function isWC(a: string, b: string): boolean {
	return complement(a) === b;
}

/**
 * Compute the ΔG of the most stable hairpin in a short oligo.
 * Scans all possible stem-loop configurations.
 *
 * @returns ΔG in kcal/mol (negative = stable hairpin present = bad)
 */
export function calcHairpinDG(seq: string): number {
	const s = seq.toUpperCase();
	const n = s.length;
	let bestDG = 0; // 0 = no structure

	// Minimum stem length: 3 bp. Minimum loop: 3 nt.
	for (let stemStart = 0; stemStart < n - 6; stemStart++) {
		for (let stemLen = 3; stemLen <= Math.floor((n - stemStart) / 2) - 1; stemLen++) {
			const loopStart = stemStart + stemLen;
			const stemEnd = n - stemStart - stemLen; // exclusive end of second stem arm

			// Check that the two stem halves are complementary
			let stemOk = true;
			for (let k = 0; k < stemLen; k++) {
				const b1 = s[stemStart + k];
				const b2 = s[n - 1 - stemStart - k];
				if (!b1 || !b2 || !isWC(b1, b2)) {
					stemOk = false;
					break;
				}
			}
			if (!stemOk) continue;

			const loopSize = stemEnd - loopStart;
			if (loopSize < 3) continue;

			// Compute stem ΔG using NN parameters for the stem sequence
			const stemSeq = s.slice(stemStart, stemStart + stemLen);
			const { dH, dS } = calcNNThermo(stemSeq);
			const stemDG = dH - 310.15 * (dS / 1000); // ΔG at 37°C

			// Total ΔG = stem stability + loop penalty
			const loopDG = hairpinLoopPenalty(loopSize);
			const totalDG = stemDG + loopDG;

			if (totalDG < bestDG) bestDG = totalDG;
		}
	}

	return bestDG;
}

/**
 * Compute the ΔG of the most stable 3′-end self-dimer.
 * Scans the last `scanLen` bases of the 3′ end against the full sequence.
 * This is the critical check for primer extension artifacts.
 *
 * @param seq       Primer sequence 5′→3′
 * @param scanLen   How many bases from the 3′ end to check. Default: 10
 * @returns ΔG in kcal/mol (negative = stable 3′ dimer)
 */
export function calcSelfDimerDG(seq: string, scanLen = 10): number {
	const s = seq.toUpperCase();
	const n = s.length;
	const threeEnd = s.slice(Math.max(0, n - scanLen));
	return calcHeteroDimerDG(threeEnd, s);
}

/**
 * Compute the ΔG of the most stable 3′-end heterodimer between two primers.
 * Aligns the 3′ end of seq1 against seq2 and finds the best complementary run.
 *
 * @returns ΔG in kcal/mol (negative = stable dimer)
 */
export function calcHeteroDimerDG(seq1: string, seq2: string, scanLen = 10): number {
	const s1 = seq1.toUpperCase();
	// Reverse seq2 so we can slide the 3' end of s1 along it
	const s2rc = seq2
		.toUpperCase()
		.split("")
		.reverse()
		.map((b) => complement(b))
		.join("");

	const probe = s1.slice(Math.max(0, s1.length - scanLen));
	let bestDG = 0;

	// Slide probe (3' end of seq1) along s2rc
	for (let offset = -(probe.length - 1); offset < s2rc.length; offset++) {
		let run = "";
		for (let i = 0; i < probe.length; i++) {
			const j = offset + i;
			const b1 = probe[i];
			const b2 = j >= 0 && j < s2rc.length ? s2rc[j] : undefined;
			if (b1 && b2 && isWC(b1, b2)) {
				run += b1;
			} else {
				// Evaluate accumulated run
				if (run.length >= 2) {
					const { dH, dS } = calcNNThermo(run);
					const dG = dH - 310.15 * (dS / 1000);
					if (dG < bestDG) bestDG = dG;
				}
				run = "";
			}
		}
		if (run.length >= 2) {
			const { dH, dS } = calcNNThermo(run);
			const dG = dH - 310.15 * (dS / 1000);
			if (dG < bestDG) bestDG = dG;
		}
	}

	return bestDG;
}
