/**
 * GWAS Data Simulator
 * Generates realistic GWAS summary statistics for demonstration
 */

class GWASSimulator {
    constructor(options = {}) {
        this.numSNPs = options.numSNPs || 50000;
        this.numChromosomes = options.numChromosomes || 22;
        this.significantLoci = options.significantLoci || 8;
        this.genomicInflation = options.genomicInflation || 1.02;
        
        // Chromosome lengths (in Mb, approximate)
        this.chrLengths = [
            249, 243, 198, 191, 182, 171, 159, 145, 138, 134,
            135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51
        ];
        
        // Sample gene names for annotation
        this.sampleGenes = [
            'APOE', 'BRCA1', 'BRCA2', 'TP53', 'EGFR', 'KRAS', 'MYC', 'PTEN',
            'RB1', 'APC', 'CDKN2A', 'PIK3CA', 'BRAF', 'NRAS', 'CTNNB1', 'SMAD4',
            'VHL', 'NF1', 'NF2', 'TSC1', 'TSC2', 'STK11', 'MLH1', 'MSH2',
            'GATA3', 'FOXA1', 'ESR1', 'AR', 'NOTCH1', 'JAK2', 'BCR', 'ABL1',
            'FLT3', 'NPM1', 'DNMT3A', 'IDH1', 'IDH2', 'TET2', 'ASXL1', 'EZH2',
            'RUNX1', 'CEBPA', 'WT1', 'KIT', 'PDGFRA', 'CSF3R', 'CALR', 'MPL',
            'SH2B3', 'TCF7L2', 'PPARG', 'KCNJ11', 'HNF1A', 'HNF4A', 'GCK', 'INS',
            'LDLR', 'PCSK9', 'APOB', 'LPL', 'CETP', 'LIPC', 'ANGPTL3', 'ANGPTL4'
        ];
    }

    /**
     * Generate complete GWAS dataset
     */
    generate() {
        console.log('Generating GWAS data...');
        const startTime = performance.now();
        
        const snps = [];
        const snpsPerChr = Math.floor(this.numSNPs / this.numChromosomes);
        
        // Determine positions for significant loci
        const significantPositions = this.generateSignificantLoci();
        
        let snpIndex = 0;
        
        for (let chr = 1; chr <= this.numChromosomes; chr++) {
            const chrLength = this.chrLengths[chr - 1] * 1e6; // Convert to bp
            const numSNPsThisChr = chr === this.numChromosomes ? 
                this.numSNPs - snpIndex : snpsPerChr;
            
            for (let i = 0; i < numSNPsThisChr; i++) {
                const position = Math.floor(Math.random() * chrLength);
                const snp = this.generateSNP(chr, position, snpIndex, significantPositions);
                snps.push(snp);
                snpIndex++;
            }
        }
        
        // Sort by chromosome and position
        snps.sort((a, b) => {
            if (a.chr !== b.chr) return a.chr - b.chr;
            return a.pos - b.pos;
        });
        
        // Calculate cumulative positions for Manhattan plot
        this.addCumulativePositions(snps);
        
        const endTime = performance.now();
        console.log(`Generated ${snps.length} SNPs in ${(endTime - startTime).toFixed(0)}ms`);
        
        return {
            snps,
            metadata: this.calculateMetadata(snps),
            significantLoci: this.extractSignificantLoci(snps)
        };
    }

    /**
     * Generate positions for significant loci
     */
    generateSignificantLoci() {
        const positions = [];
        const usedChromosomes = new Set();
        
        for (let i = 0; i < this.significantLoci; i++) {
            let chr;
            // Try to spread significant loci across chromosomes
            do {
                chr = Math.floor(Math.random() * this.numChromosomes) + 1;
            } while (usedChromosomes.has(chr) && usedChromosomes.size < this.numChromosomes);
            
            usedChromosomes.add(chr);
            
            const chrLength = this.chrLengths[chr - 1] * 1e6;
            const pos = Math.floor(Math.random() * chrLength);
            const effectSize = 0.1 + Math.random() * 0.4; // Effect size between 0.1 and 0.5
            const radius = 500000 + Math.random() * 1000000; // LD region 0.5-1.5 Mb
            
            positions.push({ chr, pos, effectSize, radius });
        }
        
        return positions;
    }

    /**
     * Generate a single SNP
     */
    generateSNP(chr, pos, index, significantPositions) {
        // Check if near a significant locus
        const nearbyLocus = significantPositions.find(locus => 
            locus.chr === chr && Math.abs(locus.pos - pos) < locus.radius
        );
        
        let pvalue, beta, maf;
        
        if (nearbyLocus) {
            // Generate correlated p-value based on distance from lead SNP
            const distance = Math.abs(nearbyLocus.pos - pos);
            const correlation = Math.exp(-distance / (nearbyLocus.radius * 0.3));
            
            // Lead SNP has very low p-value, others are correlated
            if (distance < 10000) {
                // Lead SNP
                pvalue = Math.pow(10, -(8 + Math.random() * 4)); // 1e-8 to 1e-12
                beta = nearbyLocus.effectSize * (Math.random() > 0.5 ? 1 : -1);
            } else {
                // LD SNPs
                const basePval = -Math.log10(Math.pow(10, -(8 + Math.random() * 2)));
                pvalue = Math.pow(10, -(basePval * correlation + Math.random()));
                beta = nearbyLocus.effectSize * correlation * (Math.random() > 0.5 ? 1 : -1);
            }
            maf = 0.05 + Math.random() * 0.4;
        } else {
            // Null distribution with genomic inflation
            const chi2 = this.generateInflatedChi2();
            pvalue = this.chi2ToPvalue(chi2);
            beta = (Math.random() - 0.5) * 0.1;
            maf = 0.01 + Math.random() * 0.49;
        }
        
        // Generate alleles
        const alleles = ['A', 'C', 'G', 'T'];
        const refIndex = Math.floor(Math.random() * 4);
        let altIndex;
        do {
            altIndex = Math.floor(Math.random() * 4);
        } while (altIndex === refIndex);
        
        return {
            id: `rs${10000000 + index}`,
            chr,
            pos,
            ref: alleles[refIndex],
            alt: alleles[altIndex],
            pvalue,
            logp: -Math.log10(pvalue),
            beta,
            se: Math.abs(beta) / (Math.sqrt(-2 * Math.log(pvalue)) + 0.1),
            maf,
            gene: this.assignNearestGene(chr, pos),
            isSignificant: pvalue < 5e-8
        };
    }

    /**
     * Generate chi-squared value with genomic inflation
     */
    generateInflatedChi2() {
        // Generate chi-squared(1) using Box-Muller transform
        const u1 = Math.random();
        const u2 = Math.random();
        const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        const chi2 = z * z;
        
        // Apply genomic inflation
        return chi2 * this.genomicInflation;
    }

    /**
     * Convert chi-squared to p-value
     */
    chi2ToPvalue(chi2) {
        // Approximation for chi-squared(1) CDF
        const z = Math.sqrt(chi2);
        const p = 2 * (1 - this.normalCDF(z));
        return Math.max(p, 1e-300); // Prevent underflow
    }

    /**
     * Standard normal CDF approximation
     */
    normalCDF(x) {
        const a1 = 0.254829592;
        const a2 = -0.284496736;
        const a3 = 1.421413741;
        const a4 = -1.453152027;
        const a5 = 1.061405429;
        const p = 0.3275911;
        
        const sign = x < 0 ? -1 : 1;
        x = Math.abs(x) / Math.sqrt(2);
        
        const t = 1.0 / (1.0 + p * x);
        const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
        
        return 0.5 * (1.0 + sign * y);
    }

    /**
     * Assign nearest gene (simulated)
     */
    assignNearestGene(chr, pos) {
        // Pseudo-random but deterministic gene assignment
        const seed = chr * 1000000 + Math.floor(pos / 1000000);
        const index = seed % this.sampleGenes.length;
        return this.sampleGenes[index];
    }

    /**
     * Add cumulative positions for Manhattan plot x-axis
     */
    addCumulativePositions(snps) {
        let cumOffset = 0;
        let currentChr = 0;
        const chrOffsets = {};
        
        // Calculate chromosome offsets
        for (let chr = 1; chr <= this.numChromosomes; chr++) {
            chrOffsets[chr] = cumOffset;
            cumOffset += this.chrLengths[chr - 1] * 1e6;
        }
        
        // Add cumulative position to each SNP
        snps.forEach(snp => {
            snp.cumPos = chrOffsets[snp.chr] + snp.pos;
        });
        
        // Store chromosome info for axis
        this.chromosomeInfo = [];
        cumOffset = 0;
        for (let chr = 1; chr <= this.numChromosomes; chr++) {
            const length = this.chrLengths[chr - 1] * 1e6;
            this.chromosomeInfo.push({
                chr,
                start: cumOffset,
                end: cumOffset + length,
                center: cumOffset + length / 2
            });
            cumOffset += length;
        }
    }

    /**
     * Calculate dataset metadata
     */
    calculateMetadata(snps) {
        const pvalues = snps.map(s => s.pvalue);
        const logp = snps.map(s => s.logp);
        
        // Calculate lambda (genomic inflation factor)
        const sortedPvals = [...pvalues].sort((a, b) => a - b);
        const medianIndex = Math.floor(sortedPvals.length / 2);
        const medianPval = sortedPvals[medianIndex];
        const expectedMedian = 0.5;
        
        // Lambda = observed_median_chi2 / expected_median_chi2
        const observedChi2 = this.pvalueToChiSquared(medianPval);
        const expectedChi2 = this.pvalueToChiSquared(expectedMedian);
        const lambda = observedChi2 / expectedChi2;
        
        // Count significant SNPs
        const significantCount = snps.filter(s => s.pvalue < 5e-8).length;
        const suggestiveCount = snps.filter(s => s.pvalue < 1e-5 && s.pvalue >= 5e-8).length;
        
        return {
            totalSNPs: snps.length,
            significantSNPs: significantCount,
            suggestiveSNPs: suggestiveCount,
            lambda: lambda,
            lambda1000: 1 + (lambda - 1) * (1/snps.length + 1/1000) / (1/snps.length),
            maxLogP: Math.max(...logp),
            minPvalue: Math.min(...pvalues),
            chromosomeInfo: this.chromosomeInfo
        };
    }

    /**
     * Convert p-value to chi-squared
     */
    pvalueToChiSquared(p) {
        // Inverse chi-squared CDF approximation
        const z = this.inverseNormalCDF(1 - p/2);
        return z * z;
    }

    /**
     * Inverse normal CDF approximation (Abramowitz and Stegun)
     */
    inverseNormalCDF(p) {
        if (p <= 0) return -Infinity;
        if (p >= 1) return Infinity;
        if (p === 0.5) return 0;
        
        const a = [
            -3.969683028665376e+01,
            2.209460984245205e+02,
            -2.759285104469687e+02,
            1.383577518672690e+02,
            -3.066479806614716e+01,
            2.506628277459239e+00
        ];
        const b = [
            -5.447609879822406e+01,
            1.615858368580409e+02,
            -1.556989798598866e+02,
            6.680131188771972e+01,
            -1.328068155288572e+01
        ];
        const c = [
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
            4.374664141464968e+00,
            2.938163982698783e+00
        ];
        const d = [
            7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996e+00,
            3.754408661907416e+00
        ];
        
        const pLow = 0.02425;
        const pHigh = 1 - pLow;
        let q, r;
        
        if (p < pLow) {
            q = Math.sqrt(-2 * Math.log(p));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        } else if (p <= pHigh) {
            q = p - 0.5;
            r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
        } else {
            q = Math.sqrt(-2 * Math.log(1 - p));
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        }
    }

    /**
     * Extract significant loci for quick navigation
     */
    extractSignificantLoci(snps) {
        const significant = snps.filter(s => s.pvalue < 5e-8);
        
        // Cluster nearby significant SNPs
        const loci = [];
        const clusterDistance = 1000000; // 1 Mb
        
        significant.sort((a, b) => a.pvalue - b.pvalue);
        
        significant.forEach(snp => {
            const existingLocus = loci.find(l => 
                l.chr === snp.chr && Math.abs(l.leadSNP.pos - snp.pos) < clusterDistance
            );
            
            if (existingLocus) {
                existingLocus.snps.push(snp);
            } else {
                loci.push({
                    chr: snp.chr,
                    start: snp.pos - clusterDistance/2,
                    end: snp.pos + clusterDistance/2,
                    leadSNP: snp,
                    snps: [snp]
                });
            }
        });
        
        return loci.sort((a, b) => a.leadSNP.pvalue - b.leadSNP.pvalue);
    }
}

// Export for use in other modules
window.GWASSimulator = GWASSimulator;
