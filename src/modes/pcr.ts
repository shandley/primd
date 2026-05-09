import { calcTm, calcGC, reverseComplement } from "../thermodynamics/index.js";
import type { TmOptions } from "../thermodynamics/index.js";
import { calcHairpinDG, calcSelfDimerDG, calcHeteroDimerDG } from "../thermodynamics/secondary-structure.js";
import { calcAccessibilityProfile } from "../thermodynamics/accessibility.js";
import type { PCROptions, PCRResult, PrimerCandidate, PrimerPair, ThermoOptions } from "../types.js";
import { buildTmOpts, maxPolyRun } from "./utils.js";

const DEFAULT_PRODUCT_RANGE: [number, number] = [100, 1000];
const DEFAULT_PRIMER_LEN: [number, number] = [18, 27];
const DEFAULT_TM_TARGET = 60;
const DEFAULT_MAX_TM_DIFF = 3;
const DEFAULT_MAX_HAIRPIN_DG = -2.0;
const DEFAULT_MAX_SELF_DIMER_DG = -5.0;
const DEFAULT_MAX_HETERO_DIMER_DG = -5.0;
const DEFAULT_GC_RANGE: [number, number] = [0.40, 0.65];
const DEFAULT_MAX_POLY_RUN = 4;
const DEFAULT_NUM_RETURN = 5;

interface ScoringOpts {
	tmTarget: number;
	gcRange: [number, number];
	maxHairpinDG: number;
	maxSelfDimerDG: number;
}

function makePrimerCandidate(
	seq: string,
	start: number,
	end: number,
	tmOpts: TmOptions,
	primerOpts: ScoringOpts,
	templateAccessibility: number,
): PrimerCandidate {
	const tmResult = calcTm(seq, tmOpts);
	const gc = calcGC(seq);
	const gcClamp = seq[seq.length - 1] === "G" || seq[seq.length - 1] === "C";
	const hairpinDG = calcHairpinDG(seq);
	const selfDimerDG = calcSelfDimerDG(seq);
	const poly = maxPolyRun(seq);

	// Penalty: quadratic distance from Tm target + weighted structure penalties
	const tmPenalty = Math.pow(tmResult.tm - primerOpts.tmTarget, 2) * 0.5;

	// Soft GC zone: 5 percentage points inside each hard limit
	// For default gcRange=[0.40,0.65]: soft range is [0.45,0.60] — preserves existing behavior
	const GC_SOFT_MARGIN = 0.05;
	const gcSoftMin = primerOpts.gcRange[0] + GC_SOFT_MARGIN;
	const gcSoftMax = primerOpts.gcRange[1] - GC_SOFT_MARGIN;
	const gcPenalty = gc < gcSoftMin || gc > gcSoftMax ? 5 : 0;

	// Penalty kickoffs at 50% of each filter threshold (soft gradient before hard cutoff)
	// For default maxHairpinDG=-2.0: kickoff at -1.0 — preserves existing behavior
	// For default maxSelfDimerDG=-5.0: kickoff at -3.0 (0.6 ratio) — preserves existing behavior
	const hairpinKickoff = primerOpts.maxHairpinDG * 0.5;
	const dimerKickoff = primerOpts.maxSelfDimerDG * 0.6;
	const hairpinPenalty = Math.max(0, -hairpinDG - hairpinKickoff) * 3;
	const dimerPenalty = Math.max(0, -selfDimerDG - dimerKickoff) * 4;
	const polyPenalty = Math.max(0, poly - 3) * 8;
	const clampPenalty = gcClamp ? 0 : 3;
	// Accessibility penalty: quadratic, starts when site is <50% likely single-stranded
	const accessPenalty = Math.pow(Math.max(0, 0.5 - templateAccessibility), 2) * 80;

	const penalty = tmPenalty + gcPenalty + hairpinPenalty + dimerPenalty + polyPenalty + clampPenalty + accessPenalty;

	return { seq, start, end, len: seq.length, tm: tmResult.tm, gc, gcClamp, hairpinDG, selfDimerDG, polyRun: poly, offTarget: 0, templateAccessibility, penalty };
}

/**
 * Design PCR primer pairs flanking a target region.
 *
 * @param template  Full template sequence (5′→3′, top strand)
 * @param regionStart  0-indexed start of target region (primers placed just outside)
 * @param regionEnd    0-indexed end of target region (exclusive)
 * @param opts      Design parameters
 */
export function designPCR(
	template: string,
	regionStart: number,
	regionEnd: number,
	opts: PCROptions = {},
): PCRResult {
	const seq = template.toUpperCase();
	const productRange = opts.productSizeRange ?? DEFAULT_PRODUCT_RANGE;
	const primerLen = opts.primerLenRange ?? DEFAULT_PRIMER_LEN;
	const tmTarget = opts.tmTarget ?? DEFAULT_TM_TARGET;
	const maxTmDiff = opts.maxTmDiff ?? DEFAULT_MAX_TM_DIFF;
	const maxHairpinDG = opts.maxHairpinDG ?? DEFAULT_MAX_HAIRPIN_DG;
	const maxSelfDimerDG = opts.maxSelfDimerDG ?? DEFAULT_MAX_SELF_DIMER_DG;
	const maxHeteroDimerDG = opts.maxHeteroDimerDG ?? DEFAULT_MAX_HETERO_DIMER_DG;
	const gcRange = opts.gcRange ?? DEFAULT_GC_RANGE;
	const maxPoly = opts.maxPolyRun ?? DEFAULT_MAX_POLY_RUN;
	const numReturn = opts.numReturn ?? DEFAULT_NUM_RETURN;

	const annealTempC = opts.annealTempC ?? (tmTarget - 5);
	const minAccess = opts.minTemplateAccessibility ?? 0.05;

	const tmOpts: TmOptions = buildTmOpts(opts);

	// Pre-compute accessibility profiles once — O(n·W³) total, not per candidate.
	// Forward primers anneal to the bottom strand template; reverse primers anneal to the top strand.
	// Each profile is computed on the strand the primer physically anneals to.
	const avgPrimerLen = Math.round((primerLen[0] + primerLen[1]) / 2);
	const fwdAccessProfile = calcAccessibilityProfile(reverseComplement(seq), avgPrimerLen, { annealTempC });
	const revAccessProfile = calcAccessibilityProfile(seq, avgPrimerLen, { annealTempC });

	const scoringOpts: ScoringOpts = { tmTarget, gcRange, maxHairpinDG, maxSelfDimerDG };

	// Generate forward candidates: left of the target region
	const fwdCandidates: PrimerCandidate[] = [];
	const fwdSearchStart = Math.max(0, regionStart - productRange[1]);
	const fwdSearchEnd = regionStart;

	for (let start = fwdSearchStart; start < fwdSearchEnd; start++) {
		for (let len = primerLen[0]; len <= primerLen[1]; len++) {
			const end = start + len;
			if (end > seq.length || end > regionEnd) continue;
			const primerSeq = seq.slice(start, end);
			const gc = calcGC(primerSeq);
			if (gc < gcRange[0] || gc > gcRange[1]) continue;
			if (maxPolyRun(primerSeq) > maxPoly) continue;
			// Forward primer anneals to bottom strand; RC coordinate = seq.length - end
			const access = fwdAccessProfile[seq.length - end] ?? 1.0;
			if (access < minAccess) continue;
			const candidate = makePrimerCandidate(primerSeq, start, end, tmOpts, scoringOpts, access);
			if (Math.abs(candidate.tm - tmTarget) > 15) continue;
			if (candidate.hairpinDG < maxHairpinDG) continue;
			if (candidate.selfDimerDG < maxSelfDimerDG) continue;
			fwdCandidates.push(candidate);
		}
	}

	// Generate reverse candidates: right of the target region
	const revCandidates: PrimerCandidate[] = [];
	const revSearchStart = regionEnd;
	const revSearchEnd = Math.min(seq.length, regionEnd + productRange[1]);

	for (let start = revSearchStart; start < revSearchEnd; start++) {
		for (let len = primerLen[0]; len <= primerLen[1]; len++) {
			const end = start + len;
			if (end > seq.length) continue;
			const topSeq = seq.slice(start, end);
			const primerSeq = reverseComplement(topSeq); // rev primer is RC of top strand
			const gc = calcGC(primerSeq);
			if (gc < gcRange[0] || gc > gcRange[1]) continue;
			if (maxPolyRun(primerSeq) > maxPoly) continue;
			// Reverse primer anneals to top strand; top-strand coordinate = start
			const access = revAccessProfile[start] ?? 1.0;
			if (access < minAccess) continue;
			const candidate = makePrimerCandidate(primerSeq, start, end, tmOpts, scoringOpts, access);
			if (Math.abs(candidate.tm - tmTarget) > 15) continue;
			if (candidate.hairpinDG < maxHairpinDG) continue;
			if (candidate.selfDimerDG < maxSelfDimerDG) continue;
			revCandidates.push(candidate);
		}
	}

	if (fwdCandidates.length === 0 || revCandidates.length === 0) {
		return {
			pairs: [],
			warning:
				`No valid primers found. Try widening the product size range, relaxing Tm ` +
				`constraints, or selecting a longer region (≥ 200 bp).`,
		};
	}

	// Sort individually by penalty
	fwdCandidates.sort((a, b) => a.penalty - b.penalty);
	revCandidates.sort((a, b) => a.penalty - b.penalty);

	// Pair top candidates, filter on product size and Tm matching
	const pairs: PrimerPair[] = [];
	const topFwd = fwdCandidates.slice(0, 40);
	const topRev = revCandidates.slice(0, 40);

	for (const fwd of topFwd) {
		for (const rev of topRev) {
			const productSize = rev.end - fwd.start;
			if (productSize < productRange[0] || productSize > productRange[1]) continue;
			const tmDiff = Math.abs(fwd.tm - rev.tm);
			if (tmDiff > maxTmDiff) continue;
			const heteroDimerDG = Math.min(
				calcHeteroDimerDG(fwd.seq, rev.seq),
				calcHeteroDimerDG(rev.seq, fwd.seq),
			);
			if (heteroDimerDG < maxHeteroDimerDG) continue;
			// Pair penalty kickoff at 40% of maxHeteroDimerDG threshold
			// For default maxHeteroDimerDG=-5.0: kickoff at -2.0 — preserves existing behavior
			const heteroDimerKickoff = maxHeteroDimerDG * 0.4;
			const pairPenalty =
				fwd.penalty + rev.penalty + tmDiff * 2 + Math.max(0, -heteroDimerDG - heteroDimerKickoff) * 5;
			pairs.push({ fwd, rev, productSize, heteroDimerDG, tmDiff, penalty: pairPenalty });
		}
	}

	if (pairs.length === 0) {
		return {
			pairs: [],
			warning:
				`Candidate primers exist but no compatible pairs found. ` +
				`Try increasing maxTmDiff or widening the product size range.`,
		};
	}

	pairs.sort((a, b) => a.penalty - b.penalty);
	return { pairs: pairs.slice(0, numReturn) };
}
