// ── Thermodynamics ────────────────────────────────────────────────────────────

/**
 * Which nearest-neighbor parameter set to use for Tm calculation.
 * - "santa_lucia_1998": classic unified parameters, industry standard
 * - "dna24": updated parameters derived from 27,732 hairpin measurements
 *   (Wayment-Steele et al. 2024, Nature Communications)
 */
export type NNModel = "santa_lucia_1998" | "dna24";

/**
 * Which salt correction formula to apply.
 * - "owczarzy_2004": monovalent + simple Mg²⁺ (used by primer3)
 * - "owczarzy_2008": improved Mg²⁺ accounting for free vs. total Mg²⁺ (recommended)
 */
export type SaltModel = "owczarzy_2004" | "owczarzy_2008";

export interface ThermoOptions {
	/** Nearest-neighbor parameter set. Default: "santa_lucia_1998" */
	nnModel?: NNModel;
	/** Salt correction formula. Default: "owczarzy_2008" */
	saltModel?: SaltModel;
	/** Total oligo (primer) concentration in mol/L. Default: 250e-9 (250 nM) */
	oligoConc?: number;
	/** Monovalent cation concentration in mol/L. Default: 0.05 (50 mM) */
	monoConc?: number;
	/** Total Mg²⁺ concentration in mol/L. Default: 0.002 (2 mM, typical PCR) */
	mgConc?: number;
	/** Total dNTP concentration in mol/L (chelates free Mg²⁺). Default: 0.0008 (0.8 mM) */
	dntpConc?: number;
}

export interface ThermoResult {
	/** Melting temperature in °C */
	tm: number;
	/** ΔH in kcal/mol */
	dH: number;
	/** ΔS in cal/mol/K */
	dS: number;
	/** ΔG at 37°C in kcal/mol */
	dG37: number;
}

export interface SecondaryStructure {
	/** Most stable self-fold ΔG in kcal/mol (negative = more stable) */
	dG: number;
	/** Position of the fold (0-indexed start of the stem) */
	pos: number;
	/** Length of the stem in bp */
	stemLen: number;
}

// ── Primer candidates ─────────────────────────────────────────────────────────

export interface PrimerCandidate {
	/** Primer sequence, 5'→3' */
	seq: string;
	/** 0-indexed start in the template (on the strand the primer binds) */
	start: number;
	/** 0-indexed end (exclusive) in the template */
	end: number;
	/** Length in bases */
	len: number;
	/** Melting temperature in °C */
	tm: number;
	/** GC fraction [0,1] */
	gc: number;
	/** Whether the 3' end has a G or C (GC clamp) */
	gcClamp: boolean;
	/** Most stable hairpin ΔG (kcal/mol); less negative = better */
	hairpinDG: number;
	/** Most stable self-dimer ΔG at the 3' end (kcal/mol) */
	selfDimerDG: number;
	/** Longest poly-nucleotide run length */
	polyRun: number;
	/** Number of off-target binding sites within the template */
	offTarget: number;
	/** Overall penalty score (lower = better) */
	penalty: number;
}

export interface PrimerPair {
	fwd: PrimerCandidate;
	rev: PrimerCandidate;
	/** Predicted product size in bp */
	productSize: number;
	/** 3' heterodimer ΔG (kcal/mol) */
	heteroDimerDG: number;
	/** |Tm_fwd - Tm_rev| in °C */
	tmDiff: number;
	/** Pair penalty (lower = better) */
	penalty: number;
}

// ── PCR mode ──────────────────────────────────────────────────────────────────

export interface PCROptions extends ThermoOptions {
	/** Target amplicon size range [min, max] bp. Default: [100, 1000] */
	productSizeRange?: [number, number];
	/** Primer length range [min, max] bp. Default: [18, 27] */
	primerLenRange?: [number, number];
	/** Target Tm in °C. Default: 60 */
	tmTarget?: number;
	/** Maximum |Tm_fwd - Tm_rev| in °C. Default: 3 */
	maxTmDiff?: number;
	/** Maximum hairpin ΔG in kcal/mol (less negative = more permissive). Default: -2.0 */
	maxHairpinDG?: number;
	/** Maximum self-dimer ΔG at 3' end. Default: -5.0 */
	maxSelfDimerDG?: number;
	/** Maximum 3' heterodimer ΔG. Default: -5.0 */
	maxHeteroDimerDG?: number;
	/** GC content range [min, max] as fractions. Default: [0.40, 0.65] */
	gcRange?: [number, number];
	/** Maximum poly-nucleotide run. Default: 4 */
	maxPolyRun?: number;
	/** Number of top pairs to return. Default: 5 */
	numReturn?: number;
}

export interface PCRResult {
	pairs: PrimerPair[];
	/** Human-readable explanation if no pairs were found */
	warning?: string;
}

// ── qPCR mode ─────────────────────────────────────────────────────────────────

export interface QPCROptions extends PCROptions {
	/** Target amplicon size range. Default: [70, 200] (smaller for qPCR) */
	productSizeRange?: [number, number];
	/** Whether to score amplicon secondary structure. Default: true */
	scoreAmpliconStructure?: boolean;
}

export interface QPCRResult extends PCRResult {
	pairs: (PrimerPair & {
		/** Predicted amplicon Tm in °C */
		ampliconTm: number;
		/** Amplicon secondary structure ΔG (kcal/mol) — more negative = less efficient */
		ampliconDG: number;
		/** Predicted amplification efficiency score [0,1] */
		efficiencyScore: number;
	})[];
}

// ── LAMP mode ─────────────────────────────────────────────────────────────────

/**
 * Standard LAMP primer set (Notomi et al. 2000 geometry).
 *
 *   F3 →        B3 ←
 *   FIP (F1c + F2) →     BIP (B1c + B2) ←
 *   LoopF →  LoopB ←   (optional)
 *
 * Coordinates are 0-indexed positions on the forward (top) strand.
 */
export interface LAMPPrimerSet {
	F3: PrimerCandidate;
	B3: PrimerCandidate;
	FIP: LAMPInnerPrimer;
	BIP: LAMPInnerPrimer;
	LoopF: PrimerCandidate | null;
	LoopB: PrimerCandidate | null;
	/** Overall set penalty (lower = better) */
	penalty: number;
}

export interface LAMPInnerPrimer {
	/** Full FIP/BIP sequence: 5'→3' (sense + antisense portions joined) */
	seq: string;
	/** The F1c or B1c component (reverse complement portion) */
	part1: PrimerCandidate;
	/** The F2 or B2 component (forward-binding portion) */
	part2: PrimerCandidate;
	/** Tm of part1 */
	tm1: number;
	/** Tm of part2 */
	tm2: number;
}

export interface LAMPOptions extends ThermoOptions {
	/** Target Tm for outer primers (F3/B3). Default: 60–65°C */
	outerTmRange?: [number, number];
	/** Target Tm for inner primers (FIP/BIP parts). Default: 60–65°C */
	innerTmRange?: [number, number];
	/** Target Tm for loop primers. Default: 60–65°C */
	loopTmRange?: [number, number];
	/** Whether to design loop primers. Default: true */
	designLoops?: boolean;
	/** Number of top sets to return. Default: 3 */
	numReturn?: number;
}

export interface LAMPResult {
	sets: LAMPPrimerSet[];
	warning?: string;
}

// ── Assembly primer mode ──────────────────────────────────────────────────────

export type AssemblyMethod = "gibson" | "golden_gate";

export interface AssemblyPrimerOptions extends ThermoOptions {
	/** Assembly method determines overlap and tail design. */
	method: AssemblyMethod;
	/** Overlap length for Gibson Assembly in bp. Default: 20 */
	gibsonOverlap?: number;
	/** Restriction enzyme site to include for Golden Gate (e.g., "GGTCTC" for BsaI) */
	ggEnzymeSite?: string;
	/** Primer length for the annealing region only. Default: [18, 25] */
	annealingLenRange?: [number, number];
}

export interface AssemblyPrimerPair {
	fwd: PrimerCandidate & { tail: string; fullSeq: string };
	rev: PrimerCandidate & { tail: string; fullSeq: string };
	/** Annealing Tm of the 3' annealing portion alone */
	annealingTm: number;
	productSize: number;
}

export interface AssemblyResult {
	pairs: AssemblyPrimerPair[];
	warning?: string;
}
