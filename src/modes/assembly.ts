/**
 * designAssembly — Gibson Assembly and Golden Gate primer design.
 *
 * Primers consist of a 5′ tail (overlap/RE site) + 3′ annealing region:
 *
 *   Gibson:       [←overlap from template upstream→][annealing]
 *   Golden Gate:  [AA][GGTCTC][A][4-nt overhang from template][annealing]
 *
 * The annealing regions are selected by the same thermodynamic optimisation as
 * designPCR. Overlaps are derived from the flanking template sequence so that
 * adjacent fragments (or the vector backbone) carry the complementary ends.
 *
 * For circular plasmids spanning the origin, rotate the sequence before calling
 * so the target region is contiguous — same pattern as the Ori Web Worker.
 */

import { reverseComplement } from "../thermodynamics/index.js";
import type { AssemblyPrimerOptions, AssemblyPrimerPair, AssemblyResult, PCROptions } from "../types.js";
import { designPCR } from "./pcr.js";

const DEFAULT_GIBSON_OVERLAP = 20;
const DEFAULT_ANNEALING_LEN_RANGE: [number, number] = [18, 25];
const DEFAULT_GG_ENZYME_SITE = "GGTCTC"; // BsaI recognition sequence
const GG_LEADER = "AA"; // extra bases upstream of RE site for efficient cutting
const GG_SPACER = "A"; // 1 nt between RE site and overhang (BsaI cuts 1 nt after site)
const GG_OVERHANG_LEN = 4;

/**
 * Design primers for overlap-based cloning (Gibson Assembly or Golden Gate).
 *
 * @param template    Full template sequence (5′→3′, top strand). Should include
 *                    flanking sequence beyond the amplified region — at least
 *                    `gibsonOverlap` bp on each side (default 20 bp).
 * @param regionStart 0-indexed start of the insert to amplify (inclusive)
 * @param regionEnd   0-indexed end of the insert (exclusive)
 * @param options     Design parameters; method defaults to "gibson"
 */
export function designAssembly(
	template: string,
	regionStart: number,
	regionEnd: number,
	options?: Partial<AssemblyPrimerOptions>,
): AssemblyResult {
	const method = options?.method ?? "gibson";
	const gibsonOverlap = options?.gibsonOverlap ?? DEFAULT_GIBSON_OVERLAP;
	const annealingLenRange = options?.annealingLenRange ?? DEFAULT_ANNEALING_LEN_RANGE;
	const enzymeSite = options?.ggEnzymeSite ?? DEFAULT_GG_ENZYME_SITE;

	const seq = template.toUpperCase();
	const regionLen = regionEnd - regionStart;

	// Pass thermo options through; only set properties that are actually defined
	// (exactOptionalPropertyTypes forbids assigning T|undefined to optional T).
	const pcrOpts: PCROptions = {
		primerLenRange: annealingLenRange,
		productSizeRange: [Math.max(100, regionLen + 10), regionLen + 200],
		numReturn: 5,
	};
	if (options?.nnModel !== undefined) pcrOpts.nnModel = options.nnModel;
	if (options?.saltModel !== undefined) pcrOpts.saltModel = options.saltModel;
	if (options?.oligoConc !== undefined) pcrOpts.oligoConc = options.oligoConc;
	if (options?.monoConc !== undefined) pcrOpts.monoConc = options.monoConc;
	if (options?.mgConc !== undefined) pcrOpts.mgConc = options.mgConc;
	if (options?.dntpConc !== undefined) pcrOpts.dntpConc = options.dntpConc;

	const pcrResult = designPCR(seq, regionStart, regionEnd, pcrOpts);

	if (pcrResult.pairs.length === 0) {
		return {
			pairs: [],
			warning: pcrResult.warning ?? "No valid annealing regions found in the specified region.",
		};
	}

	let warning: string | undefined;
	const pairs: AssemblyPrimerPair[] = [];

	for (const pair of pcrResult.pairs) {
		let fwdTail: string;
		let revTail: string;

		if (method === "gibson") {
			// Fwd tail: last `gibsonOverlap` bp of the template upstream of the amplicon start.
			const upstreamStart = pair.fwd.start - gibsonOverlap;
			if (upstreamStart < 0) {
				if (!warning) {
					warning =
						`Only ${pair.fwd.start} bp of upstream sequence available; ` +
						`full ${gibsonOverlap}-bp Gibson overlap not possible on fwd primer. ` +
						"Provide a longer template or reduce gibsonOverlap.";
				}
				fwdTail = seq.slice(0, pair.fwd.start);
			} else {
				fwdTail = seq.slice(upstreamStart, pair.fwd.start);
			}

			// Rev tail: RC of the first `gibsonOverlap` bp downstream of the amplicon end.
			const downstreamEnd = pair.rev.end + gibsonOverlap;
			if (downstreamEnd > seq.length) {
				if (!warning) {
					warning =
						`Only ${seq.length - pair.rev.end} bp of downstream sequence available; ` +
						`full ${gibsonOverlap}-bp Gibson overlap not possible on rev primer. ` +
						"Provide a longer template or reduce gibsonOverlap.";
				}
				revTail = reverseComplement(seq.slice(pair.rev.end));
			} else {
				revTail = reverseComplement(seq.slice(pair.rev.end, downstreamEnd));
			}
		} else {
			// Golden Gate — use natural flanking sequence as the 4-nt junction overhang.
			// Fwd: the 4 nt immediately upstream of the annealing start.
			const fwdOverhangStart = Math.max(0, pair.fwd.start - GG_OVERHANG_LEN);
			const fwdOverhang = seq.slice(fwdOverhangStart, pair.fwd.start);
			// Rev: RC of the 4 nt immediately downstream of the annealing end.
			const revOverhang = reverseComplement(
				seq.slice(pair.rev.end, Math.min(seq.length, pair.rev.end + GG_OVERHANG_LEN)),
			);

			// Primer structure: [leader][RE site][spacer][overhang][annealing]
			const enzymeSiteRC = reverseComplement(enzymeSite);
			fwdTail = `${GG_LEADER}${enzymeSite}${GG_SPACER}${fwdOverhang}`;
			revTail = `${GG_LEADER}${enzymeSiteRC}${GG_SPACER}${revOverhang}`;
		}

		pairs.push({
			fwd: {
				...pair.fwd,
				tail: fwdTail,
				fullSeq: fwdTail + pair.fwd.seq,
			},
			rev: {
				...pair.rev,
				tail: revTail,
				fullSeq: revTail + pair.rev.seq,
			},
			annealingTm: (pair.fwd.tm + pair.rev.tm) / 2,
			productSize: pair.productSize,
		});
	}

	return { pairs, ...(warning ? { warning } : {}) };
}
