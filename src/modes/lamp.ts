/**
 * LAMP primer design (Notomi et al. 2000 geometry).
 *
 * Six-region layout on the template top strand (5'→3'):
 *   [F3]─gap─[F2]─gap─[F1]─[LoopF region]─target─[LoopB region]─[B1c]─gap─[B2c]─gap─[B3c]
 *
 * Primer sequences derived from the top strand:
 *   F3   = top[f3Start:f3End]
 *   FIP  = revcomp(top[f1Start:f1End]) + top[f2Start:f2End]   (F1c at 5′, F2 at 3′)
 *   B3   = revcomp(top[b3cStart:b3cEnd])
 *   BIP  = top[b1cStart:b1cEnd] + revcomp(top[b2cStart:b2cEnd])   (B1c at 5′, B2 at 3′)
 *   LoopF = revcomp(top[lfStart:lfEnd])   (anneals to bottom strand, same direction as FIP)
 *   LoopB = top[lbStart:lbEnd]            (anneals to top strand, same direction as BIP)
 *
 * Inner amplicon (F2_start → B2c_end): 120–200 bp.
 */

import { calcTm, calcGC, reverseComplement } from "../thermodynamics/index.js";
import type { TmOptions } from "../thermodynamics/index.js";
import { calcHairpinDG, calcSelfDimerDG } from "../thermodynamics/secondary-structure.js";
import type {
	LAMPOptions,
	LAMPResult,
	LAMPPrimerSet,
	LAMPInnerPrimer,
	PrimerCandidate,
	ThermoOptions,
} from "../types.js";
import { buildTmOpts, maxPolyRun } from "./utils.js";

// ── Geometry constants ────────────────────────────────────────────────────────

const INNER_AMP_SIZE: [number, number] = [120, 200]; // F2_start → B2c_end

const F2_LEN: [number, number] = [20, 22];
const F1_LEN: [number, number] = [18, 22];
const F3_LEN: [number, number] = [18, 22];
const B1C_LEN: [number, number] = [18, 22];
const B2C_LEN: [number, number] = [20, 22];
const B3_LEN: [number, number] = [18, 22];
const LOOP_LEN: [number, number] = [17, 21];

const F3_F2_GAP: [number, number] = [0, 20]; // F3_end → F2_start
const F2_F1_GAP: [number, number] = [0, 10]; // F2_end → F1_start
const B1C_B2C_GAP: [number, number] = [0, 10]; // B1c_end → B2c_start
const B2C_B3_GAP: [number, number] = [0, 20]; // B2c_end → B3c_start

// Per Notomi 2000 / Primer Explorer v5: outer (F3/B3) must be cooler than
// inner parts so that at inner annealing temperature the outer primers are
// fully displaced. Violating this prevents efficient LAMP amplification.
const DEFAULT_OUTER_TM: [number, number] = [58, 62];
const DEFAULT_INNER_TM: [number, number] = [63, 68];
const DEFAULT_LOOP_TM: [number, number] = [60, 65];

// ── Helpers ───────────────────────────────────────────────────────────────────

function makeCandidate(
	primerSeq: string,
	templateStart: number,
	templateEnd: number,
	tmOpts: TmOptions,
): PrimerCandidate {
	const tm = calcTm(primerSeq, tmOpts).tm;
	const gc = calcGC(primerSeq);
	const gcClamp =
		primerSeq[primerSeq.length - 1] === "G" ||
		primerSeq[primerSeq.length - 1] === "C";
	return {
		seq: primerSeq,
		start: templateStart,
		end: templateEnd,
		len: primerSeq.length,
		tm,
		gc,
		gcClamp,
		hairpinDG: calcHairpinDG(primerSeq),
		selfDimerDG: calcSelfDimerDG(primerSeq),
		polyRun: maxPolyRun(primerSeq),
		offTarget: 0,
		templateAccessibility: 1.0, // not yet computed for LAMP primers
		penalty: 0,
	};
}

function gcOk(seq: string, min = 0.35, max = 0.70): boolean {
	const gc = calcGC(seq);
	return gc >= min && gc <= max;
}

function tmOk(tm: number, range: [number, number], slack = 5): boolean {
	return tm >= range[0] - slack && tm <= range[1] + slack;
}

function tmCenter(range: [number, number]): number {
	return (range[0] + range[1]) / 2;
}

// ── Main function ─────────────────────────────────────────────────────────────

/**
 * Design a LAMP primer set flanking a target region.
 *
 * @param template    Full template sequence (5′→3′, top strand)
 * @param regionStart 0-indexed start of the region to amplify/detect
 * @param regionEnd   0-indexed end of the region (exclusive)
 * @param opts        Design parameters
 */
export function designLAMP(
	template: string,
	regionStart: number,
	regionEnd: number,
	opts: LAMPOptions = {},
): LAMPResult {
	const seq = template.toUpperCase();
	const outerTm = opts.outerTmRange ?? DEFAULT_OUTER_TM;
	const innerTm = opts.innerTmRange ?? DEFAULT_INNER_TM;
	const loopTm = opts.loopTmRange ?? DEFAULT_LOOP_TM;
	const designLoops = opts.designLoops ?? true;
	const numReturn = opts.numReturn ?? 3;

	const tmOpts = buildTmOpts(opts);

	const targetLen = regionEnd - regionStart;
	if (targetLen > INNER_AMP_SIZE[1]) {
		return {
			sets: [],
			warning:
				`Target region (${targetLen} bp) is too large for LAMP. ` +
				`The inner amplicon (F2→B2c) is capped at ${INNER_AMP_SIZE[1]} bp; ` +
				`select a target region ≤ ${INNER_AMP_SIZE[1]} bp.`,
		};
	}

	// ── Step 1: Enumerate (F2, F1) forward-inner pairs ───────────────────────
	// F2 starts before regionStart; F1 immediately follows F2.
	// Both F2 and F1 must be upstream (5') of the target region core.

	interface FwdPair {
		f2Seq: string; f2Start: number; f2End: number; f2Tm: number;
		f1Seq: string; f1Start: number; f1End: number; f1Tm: number;
		fipSeq: string;
	}

	const fwdPairs: FwdPair[] = [];
	const f2SearchStart = Math.max(0, regionStart - 100);
	const f2SearchEnd = regionStart + 5;

	for (let f2Start = f2SearchStart; f2Start < f2SearchEnd; f2Start++) {
		for (let f2Len = F2_LEN[0]; f2Len <= F2_LEN[1]; f2Len++) {
			const f2End = f2Start + f2Len;
			if (f2End > seq.length) continue;
			const f2Seq = seq.slice(f2Start, f2End);
			if (!gcOk(f2Seq)) continue;
			const f2Tm = calcTm(f2Seq, tmOpts).tm;
			if (!tmOk(f2Tm, innerTm)) continue;

			for (let f1Gap = F2_F1_GAP[0]; f1Gap <= F2_F1_GAP[1]; f1Gap++) {
				for (let f1Len = F1_LEN[0]; f1Len <= F1_LEN[1]; f1Len++) {
					const f1Start = f2End + f1Gap;
					const f1End = f1Start + f1Len;
					if (f1End > seq.length) continue;
					// F1 region must stay well before the target core
					if (f1End > regionStart + 30) continue;

					const f1TopSeq = seq.slice(f1Start, f1End);
					if (!gcOk(f1TopSeq)) continue;
					const f1Tm = calcTm(f1TopSeq, tmOpts).tm; // Tm of F1c = Tm of F1
					if (!tmOk(f1Tm, innerTm)) continue;

					// FIP = revcomp(F1) + F2  (F1c at 5′, F2 at 3′)
					fwdPairs.push({
						f2Seq, f2Start, f2End, f2Tm,
						f1Seq: f1TopSeq, f1Start, f1End, f1Tm,
						fipSeq: reverseComplement(f1TopSeq) + f2Seq,
					});
				}
			}
		}
	}

	if (fwdPairs.length === 0) {
		return {
			sets: [],
			warning: "No valid forward inner primers (FIP components) found. Try widening the region or relaxing Tm constraints.",
		};
	}

	// Keep top 60 by Tm proximity to target
	const innerTmMid = tmCenter(innerTm);
	fwdPairs.sort((a, b) => {
		const da = Math.abs((a.f2Tm + a.f1Tm) / 2 - innerTmMid);
		const db = Math.abs((b.f2Tm + b.f1Tm) / 2 - innerTmMid);
		return da - db;
	});
	const topFwd = fwdPairs.slice(0, 60);

	// ── Step 2: Enumerate (B1c, B2c) reverse-inner pairs ─────────────────────
	// B2c ends at or after regionEnd; B1c precedes B2c (toward F side).

	interface RevPair {
		b2cTopSeq: string; b2cStart: number; b2cEnd: number; b2Tm: number;
		b1cSeq: string; b1cStart: number; b1cEnd: number; b1cTm: number;
		bipSeq: string;
	}

	const revPairs: RevPair[] = [];
	const b2cEndSearchStart = regionEnd - 5;
	const b2cEndSearchEnd = Math.min(seq.length, regionEnd + 100);

	for (let b2cEnd = b2cEndSearchStart; b2cEnd <= b2cEndSearchEnd; b2cEnd++) {
		for (let b2cLen = B2C_LEN[0]; b2cLen <= B2C_LEN[1]; b2cLen++) {
			const b2cStart = b2cEnd - b2cLen;
			if (b2cStart < 0) continue;
			const b2cTopSeq = seq.slice(b2cStart, b2cEnd);
			if (!gcOk(b2cTopSeq)) continue;
			const b2Tm = calcTm(b2cTopSeq, tmOpts).tm; // Tm of B2 = Tm of B2c region
			if (!tmOk(b2Tm, innerTm)) continue;

			for (let b1cGap = B1C_B2C_GAP[0]; b1cGap <= B1C_B2C_GAP[1]; b1cGap++) {
				for (let b1cLen = B1C_LEN[0]; b1cLen <= B1C_LEN[1]; b1cLen++) {
					const b1cEnd = b2cStart - b1cGap;
					const b1cStart = b1cEnd - b1cLen;
					if (b1cStart < 0) continue;
					// B1c must stay well after the target core
					if (b1cStart < regionEnd - 30) continue;

					const b1cSeq = seq.slice(b1cStart, b1cEnd);
					if (!gcOk(b1cSeq)) continue;
					const b1cTm = calcTm(b1cSeq, tmOpts).tm;
					if (!tmOk(b1cTm, innerTm)) continue;

					// BIP = B1c + revcomp(B2c)  (B1c at 5′, B2 at 3′)
					revPairs.push({
						b2cTopSeq, b2cStart, b2cEnd, b2Tm,
						b1cSeq, b1cStart, b1cEnd, b1cTm,
						bipSeq: b1cSeq + reverseComplement(b2cTopSeq),
					});
				}
			}
		}
	}

	if (revPairs.length === 0) {
		return {
			sets: [],
			warning: "No valid reverse inner primers (BIP components) found. Try widening the region or relaxing Tm constraints.",
		};
	}

	revPairs.sort((a, b) => {
		const da = Math.abs((a.b1cTm + a.b2Tm) / 2 - innerTmMid);
		const db = Math.abs((b.b1cTm + b.b2Tm) / 2 - innerTmMid);
		return da - db;
	});
	const topRev = revPairs.slice(0, 60);

	// ── Step 3: Pair inner primers by inner-amplicon size ─────────────────────

	interface InnerQuad {
		fwd: FwdPair;
		rev: RevPair;
		innerAmpSize: number;
		targetCore: number; // distance F1_end → B1c_start
	}

	const innerQuads: InnerQuad[] = [];

	for (const fwd of topFwd) {
		for (const rev of topRev) {
			if (fwd.f1End >= rev.b1cStart) continue; // F1 must precede B1c
			const targetCore = rev.b1cStart - fwd.f1End;
			if (targetCore < 10) continue; // need at least a token gap
			const innerAmpSize = rev.b2cEnd - fwd.f2Start;
			if (innerAmpSize < INNER_AMP_SIZE[0] || innerAmpSize > INNER_AMP_SIZE[1]) continue;
			// The target region of interest must lie within F2_start..B2c_end
			if (fwd.f2Start > regionStart) continue;
			if (rev.b2cEnd < regionEnd) continue;
			innerQuads.push({ fwd, rev, innerAmpSize, targetCore });
		}
	}

	if (innerQuads.length === 0) {
		return {
			sets: [],
			warning:
				"No valid inner primer combinations meet the geometric constraints. " +
				"The region may be too large for LAMP (inner amplicon max 200 bp).",
		};
	}

	// Score quads by Tm quality; keep top 30
	innerQuads.sort((a, b) => {
		const aTmDev =
			Math.abs(a.fwd.f2Tm - innerTmMid) +
			Math.abs(a.fwd.f1Tm - innerTmMid) +
			Math.abs(a.rev.b1cTm - innerTmMid) +
			Math.abs(a.rev.b2Tm - innerTmMid);
		const bTmDev =
			Math.abs(b.fwd.f2Tm - innerTmMid) +
			Math.abs(b.fwd.f1Tm - innerTmMid) +
			Math.abs(b.rev.b1cTm - innerTmMid) +
			Math.abs(b.rev.b2Tm - innerTmMid);
		return aTmDev - bTmDev;
	});
	const topQuads = innerQuads.slice(0, 30);

	// ── Step 4: Find F3 and B3 for each quad ─────────────────────────────────

	const outerTmMid = tmCenter(outerTm);
	const sets: LAMPPrimerSet[] = [];

	for (const quad of topQuads) {
		// F3: upstream of F2, forward direction
		const f3EndMax = quad.fwd.f2Start - F3_F2_GAP[0];
		const f3StartMin = Math.max(0, quad.fwd.f2Start - F3_F2_GAP[1] - F3_LEN[1]);

		let bestF3: PrimerCandidate | null = null;
		let bestF3Dev = Infinity;

		for (let f3Start = f3StartMin; f3Start < f3EndMax; f3Start++) {
			for (let f3Len = F3_LEN[0]; f3Len <= F3_LEN[1]; f3Len++) {
				const f3End = f3Start + f3Len;
				const gap = quad.fwd.f2Start - f3End;
				if (gap < F3_F2_GAP[0] || gap > F3_F2_GAP[1]) continue;
				if (f3End > seq.length) continue;
				const f3Seq = seq.slice(f3Start, f3End);
				if (!gcOk(f3Seq)) continue;
				const f3Tm = calcTm(f3Seq, tmOpts).tm;
				if (!tmOk(f3Tm, outerTm)) continue;
				const dev = Math.abs(f3Tm - outerTmMid);
				if (dev < bestF3Dev) {
					bestF3Dev = dev;
					bestF3 = makeCandidate(f3Seq, f3Start, f3End, tmOpts);
				}
			}
		}

		if (!bestF3) continue;

		// B3: downstream of B2c, reverse direction (RC of top strand)
		const b3StartMin = quad.rev.b2cEnd + B2C_B3_GAP[0];
		const b3StartMax = Math.min(seq.length - B3_LEN[0], quad.rev.b2cEnd + B2C_B3_GAP[1]);

		let bestB3: PrimerCandidate | null = null;
		let bestB3Dev = Infinity;

		for (let b3cStart = b3StartMin; b3cStart <= b3StartMax; b3cStart++) {
			const gap = b3cStart - quad.rev.b2cEnd;
			if (gap < B2C_B3_GAP[0] || gap > B2C_B3_GAP[1]) continue;
			for (let b3Len = B3_LEN[0]; b3Len <= B3_LEN[1]; b3Len++) {
				const b3cEnd = b3cStart + b3Len;
				if (b3cEnd > seq.length) continue;
				const b3cTopSeq = seq.slice(b3cStart, b3cEnd);
				const b3Seq = reverseComplement(b3cTopSeq);
				if (!gcOk(b3Seq)) continue;
				const b3Tm = calcTm(b3Seq, tmOpts).tm;
				if (!tmOk(b3Tm, outerTm)) continue;
				const dev = Math.abs(b3Tm - outerTmMid);
				if (dev < bestB3Dev) {
					bestB3Dev = dev;
					bestB3 = makeCandidate(b3Seq, b3cStart, b3cEnd, tmOpts);
				}
			}
		}

		if (!bestB3) continue;

		// ── Step 5: Optional loop primers ──────────────────────────────────────
		// Target core = [F1_end, B1c_start]. Split into LoopF (first half) and LoopB (second half).
		let loopF: PrimerCandidate | null = null;
		let loopB: PrimerCandidate | null = null;

		if (designLoops && quad.targetCore >= LOOP_LEN[0] * 2 + 10) {
			const coreStart = quad.fwd.f1End;
			const coreEnd = quad.rev.b1cStart;
			const coreMid = Math.floor((coreStart + coreEnd) / 2);

			// LoopF: in [coreStart, coreMid], primer = revcomp(top strand)
			const loopTmMid = tmCenter(loopTm);
			let bestLFDev = Infinity;

			for (let lfStart = coreStart; lfStart < coreMid - LOOP_LEN[0]; lfStart++) {
				for (let lfLen = LOOP_LEN[0]; lfLen <= LOOP_LEN[1]; lfLen++) {
					const lfEnd = lfStart + lfLen;
					if (lfEnd > coreMid) continue;
					const lfTopSeq = seq.slice(lfStart, lfEnd);
					const loopFSeq = reverseComplement(lfTopSeq);
					if (!gcOk(loopFSeq)) continue;
					const lfTm = calcTm(loopFSeq, tmOpts).tm;
					if (!tmOk(lfTm, loopTm, 3)) continue;
					const dev = Math.abs(lfTm - loopTmMid);
					if (dev < bestLFDev) {
						bestLFDev = dev;
						loopF = makeCandidate(loopFSeq, lfStart, lfEnd, tmOpts);
					}
				}
			}

			// LoopB: in [coreMid, coreEnd], primer = top strand sequence
			let bestLBDev = Infinity;

			for (let lbStart = coreMid; lbStart < coreEnd - LOOP_LEN[0]; lbStart++) {
				for (let lbLen = LOOP_LEN[0]; lbLen <= LOOP_LEN[1]; lbLen++) {
					const lbEnd = lbStart + lbLen;
					if (lbEnd > coreEnd) continue;
					const loopBSeq = seq.slice(lbStart, lbEnd);
					if (!gcOk(loopBSeq)) continue;
					const lbTm = calcTm(loopBSeq, tmOpts).tm;
					if (!tmOk(lbTm, loopTm, 3)) continue;
					const dev = Math.abs(lbTm - loopTmMid);
					if (dev < bestLBDev) {
						bestLBDev = dev;
						loopB = makeCandidate(loopBSeq, lbStart, lbEnd, tmOpts);
					}
				}
			}
		}

		// ── Step 6: Assemble primer set and compute penalty ─────────────────────
		const f2Cand = makeCandidate(quad.fwd.f2Seq, quad.fwd.f2Start, quad.fwd.f2End, tmOpts);
		// F1c: primer sequence is revcomp of F1 region; use F1 template coords
		const f1cCand = makeCandidate(
			reverseComplement(quad.fwd.f1Seq),
			quad.fwd.f1Start,
			quad.fwd.f1End,
			tmOpts,
		);
		const b1cCand = makeCandidate(quad.rev.b1cSeq, quad.rev.b1cStart, quad.rev.b1cEnd, tmOpts);
		const b2Cand = makeCandidate(
			reverseComplement(quad.rev.b2cTopSeq),
			quad.rev.b2cStart,
			quad.rev.b2cEnd,
			tmOpts,
		);

		const fip: LAMPInnerPrimer = {
			seq: quad.fwd.fipSeq,
			part1: f1cCand,
			part2: f2Cand,
			tm1: quad.fwd.f1Tm,
			tm2: quad.fwd.f2Tm,
		};

		const bip: LAMPInnerPrimer = {
			seq: quad.rev.bipSeq,
			part1: b1cCand,
			part2: b2Cand,
			tm1: quad.rev.b1cTm,
			tm2: quad.rev.b2Tm,
		};

		const penalty =
			Math.pow(bestF3.tm - outerTmMid, 2) * 0.5 +
			Math.pow(bestB3.tm - outerTmMid, 2) * 0.5 +
			Math.pow(quad.fwd.f2Tm - innerTmMid, 2) * 0.3 +
			Math.pow(quad.fwd.f1Tm - innerTmMid, 2) * 0.3 +
			Math.pow(quad.rev.b1cTm - innerTmMid, 2) * 0.3 +
			Math.pow(quad.rev.b2Tm - innerTmMid, 2) * 0.3 +
			(bestF3.hairpinDG < -2 ? (-bestF3.hairpinDG - 2) * 3 : 0) +
			(bestB3.hairpinDG < -2 ? (-bestB3.hairpinDG - 2) * 3 : 0) +
			(bestF3.gcClamp ? 0 : 2) +
			(bestB3.gcClamp ? 0 : 2) +
			(loopF ? 0 : 5) + (loopB ? 0 : 5); // reward having loop primers

		sets.push({
			F3: bestF3,
			B3: bestB3,
			FIP: fip,
			BIP: bip,
			LoopF: loopF,
			LoopB: loopB,
			penalty,
		});
	}

	if (sets.length === 0) {
		return {
			sets: [],
			warning:
				"No valid outer primers (F3/B3) found within the required gap distances. " +
				"Try providing a longer template sequence upstream/downstream of the region.",
		};
	}

	sets.sort((a, b) => a.penalty - b.penalty);
	return { sets: sets.slice(0, numReturn) };
}
