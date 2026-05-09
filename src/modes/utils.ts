import type { ThermoOptions } from "../types.js";
import type { TmOptions } from "../thermodynamics/index.js";

export function buildTmOpts(opts: ThermoOptions): TmOptions {
	const out: TmOptions = {};
	if (opts.nnModel !== undefined) out.nnModel = opts.nnModel;
	if (opts.saltModel !== undefined) out.saltModel = opts.saltModel;
	if (opts.oligoConc !== undefined) out.oligoConc = opts.oligoConc;
	if (opts.monoConc !== undefined) out.monoConc = opts.monoConc;
	if (opts.mgConc !== undefined) out.mgConc = opts.mgConc;
	if (opts.dntpConc !== undefined) out.dntpConc = opts.dntpConc;
	return out;
}

export function maxPolyRun(seq: string): number {
	let max = 1;
	let cur = 1;
	for (let i = 1; i < seq.length; i++) {
		cur = seq[i] === seq[i - 1] ? cur + 1 : 1;
		if (cur > max) max = cur;
	}
	return max;
}
