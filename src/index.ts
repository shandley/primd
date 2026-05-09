/**
 * primd — Modern, web-native primer design library
 *
 * Modules:
 *   thermodynamics  calcTm, calcGC, reverseComplement, calcHairpinDG, etc.
 *   modes/pcr       designPCR — standard PCR primer pairs
 *   modes/lamp      designLAMP — LAMP primer sets (coming soon)
 *   modes/assembly  designAssembly — Gibson/Golden Gate primers (coming soon)
 *   modes/qpcr      designQPCR — qPCR primer pairs with efficiency scoring (coming soon)
 */

// Thermodynamics
export { calcTm, calcGC, reverseComplement } from "./thermodynamics/index.js";
export { calcNNThermo } from "./thermodynamics/nearest-neighbor.js";
export { applySaltCorrection } from "./thermodynamics/salt-correction.js";
export { calcHairpinDG, calcSelfDimerDG, calcHeteroDimerDG } from "./thermodynamics/secondary-structure.js";

// Design modes
export { designPCR } from "./modes/pcr.js";

// Types
export type {
	NNModel,
	SaltModel,
	ThermoOptions,
	ThermoResult,
	SecondaryStructure,
	PrimerCandidate,
	PrimerPair,
	PCROptions,
	PCRResult,
	QPCROptions,
	QPCRResult,
	LAMPPrimerSet,
	LAMPInnerPrimer,
	LAMPOptions,
	LAMPResult,
	AssemblyMethod,
	AssemblyPrimerOptions,
	AssemblyPrimerPair,
	AssemblyResult,
} from "./types.js";
