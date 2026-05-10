/**
 * primd — Modern, web-native primer design library
 *
 * Modules:
 *   thermodynamics  calcTm, calcGC, reverseComplement, calcHairpinDG, etc. (SantaLucia 1998 + Owczarzy 2004/2008)
 *   modes/pcr       designPCR — standard PCR primer pairs
 *   modes/lamp      designLAMP — LAMP primer sets (F3/FIP/B3/BIP + optional loop primers)
 *   modes/assembly  designAssembly — Gibson Assembly and Golden Gate primers
 *   modes/qpcr      designQPCR — qPCR primer pairs with efficiency scoring
 */

// Thermodynamics
export { calcTm, calcGC, reverseComplement } from "./thermodynamics/index.js";
export { calcNNThermo } from "./thermodynamics/nearest-neighbor.js";
export { applySaltCorrection } from "./thermodynamics/salt-correction.js";
export { calcHairpinDG, calcSelfDimerDG, calcHeteroDimerDG } from "./thermodynamics/secondary-structure.js";
export { calcAccessibility, calcAccessibilityProfile } from "./thermodynamics/accessibility.js";
export type { AccessibilityOpts } from "./thermodynamics/accessibility.js";

// Design modes
export { designPCR } from "./modes/pcr.js";
export { designQPCR } from "./modes/qpcr.js";
export { designLAMP } from "./modes/lamp.js";
export { designAssembly } from "./modes/assembly.js";

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
