import { calcGC } from "../thermodynamics/index.js";
import { calcHairpinDG } from "../thermodynamics/secondary-structure.js";
import type { PrimerPair, QPCROptions, QPCRResult } from "../types.js";
import { designPCR } from "./pcr.js";

// qPCR-specific defaults — tighter amplicon size, same Tm target
const QPCR_PRODUCT_RANGE: [number, number] = [70, 200];

/**
 * Amplicon Tm estimate using the Marmur-Schildkraut-Doty formula for long duplexes.
 * More appropriate for PCR products than the oligonucleotide NN model.
 *
 * Tm = 81.5 + 16.6 × log₁₀([Na⁺]) + 0.41 × %GC − 675 / length
 *
 * Reference: Marmur J & Doty P (1962) J Mol Biol 5:109–118
 */
function calcAmpliconTm(seq: string, monoConc: number): number {
	const len = seq.length;
	const gcPct = calcGC(seq) * 100;
	return 81.5 + 16.6 * Math.log10(monoConc) + 0.41 * gcPct - 675 / len;
}

/**
 * Efficiency score [0,1] based on amplicon properties relevant to qPCR.
 *
 * Scoring components (equal weight):
 *   - Size: sweet spot 80–150 bp; penalties below 70 or above 180
 *   - GC%: optimal 45–60%; penalties outside 35–70%
 *   - Secondary structure: penalises hairpins with ΔG < −3 kcal/mol
 *
 * A score ≥ 0.8 indicates a high-quality qPCR amplicon.
 */
function calcEfficiencyScore(seq: string, ampliconDG: number): number {
	const len = seq.length;
	const gc = calcGC(seq);

	// Size score — optimal 80–150 bp
	let sizeScore: number;
	if (len >= 80 && len <= 150) {
		sizeScore = 1.0;
	} else if (len < 80) {
		sizeScore = Math.max(0.2, 1.0 - (80 - len) / 60);
	} else {
		sizeScore = Math.max(0.1, 1.0 - (len - 150) / 100);
	}

	// GC score — optimal 45–60%
	let gcScore: number;
	if (gc >= 0.45 && gc <= 0.60) {
		gcScore = 1.0;
	} else if (gc >= 0.35 && gc < 0.45) {
		gcScore = 0.5 + (gc - 0.35) * 5; // linear 0.5→1.0 from 35→45%
	} else if (gc > 0.60 && gc <= 0.70) {
		gcScore = 0.5 + (0.70 - gc) * 5; // linear 0.5→1.0 from 70→60%
	} else {
		gcScore = 0.2; // <35% or >70% — low efficiency expected
	}

	// Structure score — hairpin ΔG below −3 kcal/mol progressively penalised
	let structureScore: number;
	if (ampliconDG >= -3.0) {
		structureScore = 1.0;
	} else if (ampliconDG <= -12.0) {
		structureScore = 0.1;
	} else {
		// linear from 1.0 at −3 to 0.1 at −12
		structureScore = 1.0 - ((-ampliconDG - 3.0) / 9.0) * 0.9;
	}

	return Math.max(0, Math.min(1, (sizeScore + gcScore + structureScore) / 3));
}

type QPCRPair = QPCRResult["pairs"][number];

function augmentPair(
	pair: PrimerPair,
	template: string,
	monoConc: number,
	scoreStructure: boolean,
): QPCRPair {
	const ampliconSeq = template.slice(pair.fwd.start, pair.rev.end);
	const ampliconTm = calcAmpliconTm(ampliconSeq, monoConc);
	// calcHairpinDG on the amplicon gives the most stable intramolecular fold —
	// a strong hairpin here could block polymerase extension.
	const ampliconDG = scoreStructure ? calcHairpinDG(ampliconSeq) : 0;
	const efficiencyScore = calcEfficiencyScore(ampliconSeq, ampliconDG);
	return { ...pair, ampliconTm, ampliconDG, efficiencyScore };
}

/**
 * Design primer pairs optimised for quantitative PCR (qPCR / RT-qPCR).
 *
 * Extends designPCR with:
 *   - Smaller default amplicon range [70, 200] bp (optimal for real-time detection)
 *   - Amplicon Tm (Marmur-Schildkraut formula, appropriate for PCR products)
 *   - Amplicon secondary structure ΔG (hairpin potential at binding site)
 *   - Efficiency score [0,1] — composite of size, GC%, and structure penalties
 *
 * Pairs are re-ranked by efficiencyScore descending, with penalty as tiebreaker.
 *
 * @param template   Full template sequence (5′→3′, top strand)
 * @param regionStart  0-indexed start of target region (primers placed just outside)
 * @param regionEnd    0-indexed end of target region (exclusive)
 * @param opts       Design parameters (all PCROptions apply; productSizeRange defaults to [70,200])
 */
export function designQPCR(
	template: string,
	regionStart: number,
	regionEnd: number,
	opts: QPCROptions = {},
): QPCRResult {
	const scoreStructure = opts.scoreAmpliconStructure ?? true;
	const monoConc = opts.monoConc ?? 0.05; // 50 mM Na⁺ for Tm formula

	// Run PCR design with qPCR defaults
	const pcrResult = designPCR(template, regionStart, regionEnd, {
		productSizeRange: QPCR_PRODUCT_RANGE,
		...opts,
	});

	if (pcrResult.pairs.length === 0) {
		const result: QPCRResult = { pairs: [] };
		if (pcrResult.warning) result.warning = pcrResult.warning;
		return result;
	}

	// Augment each pair with qPCR-specific metrics
	const seq = template.toUpperCase();
	const augmented = pcrResult.pairs.map((pair) =>
		augmentPair(pair, seq, monoConc, scoreStructure),
	);

	// Re-rank: highest efficiencyScore first; ties broken by PCR penalty (lower = better)
	augmented.sort((a, b) =>
		b.efficiencyScore !== a.efficiencyScore
			? b.efficiencyScore - a.efficiencyScore
			: a.penalty - b.penalty,
	);

	const result: QPCRResult = { pairs: augmented };
	if (pcrResult.warning) result.warning = pcrResult.warning;
	return result;
}
