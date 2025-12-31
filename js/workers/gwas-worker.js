/**
 * GWAS Web Worker
 * Handles heavy computations in background thread
 * Prevents UI blocking during data generation and processing
 */

let wasmModule = null;
let gwasEngine = null;
let prsEngine = null;

// Initialize WASM if available
async function initWasm() {
    try {
        const wasmUrl = new URL('../../wasm/pkg/gwas_wasm.js', self.location.href);
        wasmModule = await import(wasmUrl.href);
        await wasmModule.default();
        return true;
    } catch (e) {
        console.log('WASM not available, using pure JS fallback');
        return false;
    }
}

// Message handler
self.onmessage = async function(e) {
    const { type, data, id } = e.data;

    try {
        let result;

        switch (type) {
            case 'init':
                const hasWasm = await initWasm();
                result = { initialized: true, wasmAvailable: hasWasm };
                break;

            case 'generateGWAS':
                result = await generateGWASData(data);
                break;

            case 'calculateStats':
                result = calculateStatistics(data);
                break;

            case 'filterSNPs':
                result = filterSNPs(data.snps, data.threshold);
                break;

            case 'calculateQQ':
                result = calculateQQData(data.pvalues);
                break;

            case 'calculatePRS':
                result = await calculatePRS(data);
                break;

            case 'sortSNPs':
                result = sortSNPs(data.snps, data.sortBy);
                break;

            default:
                throw new Error(`Unknown message type: ${type}`);
        }

        self.postMessage({ id, success: true, result });
    } catch (error) {
        self.postMessage({ id, success: false, error: error.message });
    }
};

/**
 * Generate GWAS data using WASM or JS fallback
 */
async function generateGWASData(options) {
    const startTime = performance.now();
    const {
        numSNPs = 50000,
        numChromosomes = 22,
        significantLoci = 12,
        genomicInflation = 1.02
    } = options;

    let snps;

    if (wasmModule && wasmModule.GWASEngine) {
        // Use WASM engine
        gwasEngine = new wasmModule.GWASEngine(numSNPs, numChromosomes, genomicInflation);
        const rawData = gwasEngine.generate_pvalues(significantLoci);

        // Parse flat array into SNP objects
        snps = [];
        for (let i = 0; i < rawData.length; i += 4) {
            snps.push({
                chr: rawData[i],
                pos: rawData[i + 1],
                pvalue: rawData[i + 2],
                logp: rawData[i + 3]
            });
        }
    } else {
        // Pure JS fallback
        snps = generateGWASJS(numSNPs, numChromosomes, significantLoci, genomicInflation);
    }

    // Add additional fields and sort
    const chrLengths = [
        249, 243, 198, 191, 182, 171, 159, 145, 138, 134,
        135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51
    ];

    // Calculate cumulative positions
    const chrOffsets = {};
    let cumOffset = 0;
    for (let chr = 1; chr <= numChromosomes; chr++) {
        chrOffsets[chr] = cumOffset;
        cumOffset += chrLengths[chr - 1] * 1e6;
    }

    // Add IDs, cumPos, and other fields
    snps.forEach((snp, i) => {
        snp.id = `rs${10000000 + i}`;
        snp.cumPos = chrOffsets[snp.chr] + snp.pos;
        snp.isSignificant = snp.pvalue < 5e-8;
        snp.beta = (Math.random() - 0.5) * (snp.isSignificant ? 0.3 : 0.1);
        snp.maf = 0.05 + Math.random() * 0.4;
    });

    // Sort by chromosome and position
    snps.sort((a, b) => a.chr !== b.chr ? a.chr - b.chr : a.pos - b.pos);

    // Calculate metadata
    const metadata = calculateMetadata(snps, chrLengths, numChromosomes);

    const endTime = performance.now();
    console.log(`Generated ${snps.length} SNPs in ${(endTime - startTime).toFixed(0)}ms`);

    return { snps, metadata };
}

/**
 * Pure JS GWAS generation (fallback)
 */
function generateGWASJS(numSNPs, numChromosomes, significantLoci, genomicInflation) {
    const chrLengths = [
        249, 243, 198, 191, 182, 171, 159, 145, 138, 134,
        135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51
    ];

    // Generate significant loci
    const sigPositions = [];
    for (let i = 0; i < significantLoci; i++) {
        const chr = Math.floor(Math.random() * numChromosomes) + 1;
        const pos = Math.random() * chrLengths[chr - 1] * 1e6;
        const effect = 0.1 + Math.random() * 0.4;
        const radius = 500000 + Math.random() * 1000000;
        sigPositions.push({ chr, pos, effect, radius });
    }

    const snps = [];
    const snpsPerChr = Math.floor(numSNPs / numChromosomes);

    for (let chr = 1; chr <= numChromosomes; chr++) {
        const chrLength = chrLengths[chr - 1] * 1e6;

        for (let i = 0; i < snpsPerChr; i++) {
            const pos = Math.random() * chrLength;

            // Check if near significant locus
            const nearSig = sigPositions.find(l =>
                l.chr === chr && Math.abs(l.pos - pos) < l.radius
            );

            let pvalue;
            if (nearSig) {
                const distance = Math.abs(nearSig.pos - pos);
                const correlation = Math.exp(-distance / (nearSig.radius * 0.3));

                if (distance < 10000) {
                    pvalue = Math.pow(10, -(8 + Math.random() * 4));
                } else {
                    const base = 8 + Math.random() * 2;
                    pvalue = Math.pow(10, -(base * correlation + Math.random()));
                }
            } else {
                const chi2 = generateInflatedChi2(genomicInflation);
                pvalue = chi2ToPvalue(chi2);
            }

            snps.push({
                chr,
                pos,
                pvalue,
                logp: -Math.log10(pvalue)
            });
        }
    }

    return snps;
}

function generateInflatedChi2(inflation) {
    const u1 = Math.random();
    const u2 = Math.random();
    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    return z * z * inflation;
}

function chi2ToPvalue(chi2) {
    const z = Math.sqrt(chi2);
    return Math.max(2 * (1 - normalCDF(z)), 1e-300);
}

function normalCDF(x) {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x) / Math.SQRT2;

    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.exp(-x*x);

    return 0.5 * (1 + sign * y);
}

/**
 * Calculate metadata for GWAS results
 */
function calculateMetadata(snps, chrLengths, numChromosomes) {
    // Calculate chromosome info
    const chromosomeInfo = [];
    let cumOffset = 0;
    for (let chr = 1; chr <= numChromosomes; chr++) {
        const length = chrLengths[chr - 1] * 1e6;
        chromosomeInfo.push({
            chr,
            start: cumOffset,
            end: cumOffset + length,
            center: cumOffset + length / 2
        });
        cumOffset += length;
    }

    // Calculate lambda (genomic inflation)
    const pvalues = snps.map(s => s.pvalue);
    const sortedPvals = [...pvalues].sort((a, b) => a - b);
    const medianPval = sortedPvals[Math.floor(sortedPvals.length / 2)];

    const observedChi2 = pvalueToChiSquared(medianPval);
    const expectedChi2 = pvalueToChiSquared(0.5);
    const lambda = observedChi2 / expectedChi2;

    return {
        totalSNPs: snps.length,
        significantSNPs: snps.filter(s => s.pvalue < 5e-8).length,
        suggestiveSNPs: snps.filter(s => s.pvalue < 1e-5 && s.pvalue >= 5e-8).length,
        lambda,
        maxLogP: Math.max(...snps.map(s => s.logp)),
        minPvalue: Math.min(...pvalues),
        chromosomeInfo
    };
}

function pvalueToChiSquared(p) {
    const z = inverseNormalCDF(1 - p/2);
    return z * z;
}

function inverseNormalCDF(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p === 0.5) return 0;

    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];

    const pLow = 0.02425;

    if (p < pLow) {
        const q = Math.sqrt(-2 * Math.log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= 1 - pLow) {
        const q = p - 0.5;
        const r = q * q;
        const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                   1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
        const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                   6.680131188771972e+01, -1.328068155288572e+01];
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
        const q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
}

/**
 * Calculate QQ plot data
 */
function calculateQQData(pvalues) {
    const sorted = [...pvalues].sort((a, b) => a - b);
    const n = sorted.length;

    const expected = [];
    const observed = [];

    // Sample points for large datasets
    const step = n > 10000 ? Math.floor(n / 5000) : 1;

    for (let i = 0; i < n; i += step) {
        const exp = (i + 0.5) / n;
        expected.push(-Math.log10(exp));
        observed.push(-Math.log10(sorted[i]));
    }

    return { expected, observed };
}

/**
 * Calculate PRS
 */
async function calculatePRS(data) {
    const { genotypes, weights } = data;

    if (wasmModule && wasmModule.PRSEngine && !prsEngine) {
        prsEngine = new wasmModule.PRSEngine(weights.length);
    }

    if (prsEngine) {
        const result = prsEngine.calculate(new Float64Array(genotypes));
        return {
            score: result[0],
            standardizedScore: result[1],
            percentile: result[2]
        };
    }

    // JS fallback
    const score = genotypes.reduce((sum, g, i) => sum + g * weights[i], 0);
    return { score, standardizedScore: 0, percentile: 50 };
}

/**
 * Filter SNPs by p-value threshold
 */
function filterSNPs(snps, threshold) {
    return snps.filter(s => s.pvalue < threshold);
}

/**
 * Sort SNPs by specified field
 */
function sortSNPs(snps, sortBy) {
    return [...snps].sort((a, b) => {
        switch (sortBy) {
            case 'pvalue': return a.pvalue - b.pvalue;
            case 'chr': return a.chr !== b.chr ? a.chr - b.chr : a.pos - b.pos;
            case 'logp': return b.logp - a.logp;
            default: return 0;
        }
    });
}

/**
 * Calculate basic statistics
 */
function calculateStatistics(data) {
    const { values } = data;
    const n = values.length;

    const sum = values.reduce((a, b) => a + b, 0);
    const mean = sum / n;
    const variance = values.reduce((s, v) => s + Math.pow(v - mean, 2), 0) / n;

    const sorted = [...values].sort((a, b) => a - b);
    const median = sorted[Math.floor(n / 2)];

    return { mean, median, variance, sd: Math.sqrt(variance), min: sorted[0], max: sorted[n-1] };
}
