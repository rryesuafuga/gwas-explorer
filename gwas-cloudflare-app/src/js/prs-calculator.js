/**
 * Polygenic Risk Score Calculator with TensorFlow.js
 * Browser-based ML inference for genetic risk prediction
 */

class PRSCalculator {
    constructor(options = {}) {
        this.options = {
            modelPath: null, // Path to pre-trained model
            numVariants: 50, // Number of risk variants
            ...options
        };
        
        this.model = null;
        this.isModelLoaded = false;
        this.riskVariants = [];
        this.populationStats = null;
        
        // Population-specific allele frequencies
        this.populationFrequencies = {
            AFR: { meanMAF: 0.25, sdMAF: 0.15 },
            EUR: { meanMAF: 0.20, sdMAF: 0.12 },
            EAS: { meanMAF: 0.18, sdMAF: 0.10 },
            SAS: { meanMAF: 0.22, sdMAF: 0.13 },
            AMR: { meanMAF: 0.21, sdMAF: 0.14 }
        };
        
        this.onModelLoaded = null;
        this.onError = null;
    }

    /**
     * Initialize the PRS calculator
     */
    async initialize() {
        try {
            // Generate risk variants with realistic weights
            this.generateRiskVariants();
            
            // Create or load the ML model
            await this.createModel();
            
            // Generate population statistics for percentile calculation
            this.generatePopulationStats();
            
            this.isModelLoaded = true;
            
            if (this.onModelLoaded) {
                this.onModelLoaded();
            }
            
            return true;
        } catch (error) {
            console.error('Failed to initialize PRS calculator:', error);
            if (this.onError) {
                this.onError(error);
            }
            return false;
        }
    }

    /**
     * Generate realistic risk variants
     */
    generateRiskVariants() {
        // Simulate SNPs associated with a complex trait
        // Using realistic effect sizes and allele frequencies
        const genes = [
            'APOE', 'PCSK9', 'LDLR', 'HMGCR', 'NPC1L1', 'ABCG5', 'ABCG8',
            'CETP', 'LIPC', 'LPL', 'APOB', 'ANGPTL3', 'ANGPTL4', 'SORT1',
            'TRIB1', 'GCKR', 'MLXIPL', 'FADS1', 'FADS2', 'ELOVL2',
            'TCF7L2', 'PPARG', 'KCNJ11', 'SLC30A8', 'CDKAL1', 'IGF2BP2',
            'MTNR1B', 'GLIS3', 'PROX1', 'ADCY5', 'GCK', 'HNF1A', 'HNF4A',
            'INS', 'IRS1', 'FTO', 'MC4R', 'BDNF', 'SEC16B', 'GNPDA2',
            'TMEM18', 'KCTD15', 'MTCH2', 'NEGR1', 'SH2B1', 'FAIM2',
            'LRRN6C', 'PRKD1', 'GPRC5B', 'MAP2K5'
        ];
        
        this.riskVariants = [];
        
        for (let i = 0; i < this.options.numVariants; i++) {
            const gene = genes[i % genes.length];
            const chr = Math.floor(Math.random() * 22) + 1;
            
            // Effect size follows realistic distribution
            // Most variants have small effects, few have larger effects
            const effectSize = this.generateEffectSize();
            
            // MAF typically 0.01 - 0.49
            const maf = 0.01 + Math.random() * 0.48;
            
            this.riskVariants.push({
                id: `rs${10000000 + i}`,
                gene,
                chr,
                weight: effectSize,
                maf,
                riskAllele: Math.random() > 0.5 ? 'A' : 'G',
                protectiveAllele: Math.random() > 0.5 ? 'T' : 'C'
            });
        }
        
        // Sort by absolute weight (largest effects first)
        this.riskVariants.sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight));
    }

    /**
     * Generate realistic effect size
     */
    generateEffectSize() {
        // Effect sizes follow exponential-like distribution
        // Most are small (0.01-0.05), few are large (0.1-0.3)
        const u = Math.random();
        
        if (u < 0.7) {
            // Small effect
            return (Math.random() * 0.05 + 0.01) * (Math.random() > 0.5 ? 1 : -1);
        } else if (u < 0.9) {
            // Medium effect
            return (Math.random() * 0.1 + 0.05) * (Math.random() > 0.5 ? 1 : -1);
        } else {
            // Large effect
            return (Math.random() * 0.2 + 0.1) * (Math.random() > 0.5 ? 1 : -1);
        }
    }

    /**
     * Create TensorFlow.js model
     */
    async createModel() {
        // Simple neural network for PRS prediction
        // In practice, this would be a pre-trained model loaded from file
        
        this.model = tf.sequential({
            layers: [
                // Input layer
                tf.layers.dense({
                    inputShape: [this.options.numVariants],
                    units: 32,
                    activation: 'relu',
                    kernelInitializer: 'glorotNormal'
                }),
                
                // Hidden layer with dropout
                tf.layers.dropout({ rate: 0.2 }),
                
                tf.layers.dense({
                    units: 16,
                    activation: 'relu'
                }),
                
                // Output layer (single PRS value)
                tf.layers.dense({
                    units: 1,
                    activation: 'linear'
                })
            ]
        });
        
        // Compile model
        this.model.compile({
            optimizer: tf.train.adam(0.001),
            loss: 'meanSquaredError'
        });
        
        // Initialize with realistic weights based on our risk variants
        await this.initializeModelWeights();
        
        console.log('TensorFlow.js model created');
        console.log('Model summary:');
        this.model.summary();
    }

    /**
     * Initialize model weights to produce realistic PRS
     */
    async initializeModelWeights() {
        // Get the weights of the first layer
        const weights = this.model.layers[0].getWeights();
        
        if (weights.length > 0) {
            // Create weight matrix based on our risk variant weights
            const numInputs = this.options.numVariants;
            const numUnits = 32;
            
            // Initialize with small random values + contribution from variant weights
            const weightData = new Float32Array(numInputs * numUnits);
            
            for (let i = 0; i < numInputs; i++) {
                const variantWeight = this.riskVariants[i]?.weight || 0;
                for (let j = 0; j < numUnits; j++) {
                    weightData[i * numUnits + j] = (Math.random() - 0.5) * 0.1 + variantWeight * 0.5;
                }
            }
            
            const newWeights = tf.tensor2d(weightData, [numInputs, numUnits]);
            const biases = weights[1];
            
            this.model.layers[0].setWeights([newWeights, biases]);
        }
    }

    /**
     * Generate population statistics for percentile calculation
     */
    generatePopulationStats() {
        // Simulate population PRS distribution
        const numSamples = 10000;
        const scores = [];
        
        for (let i = 0; i < numSamples; i++) {
            // Generate random genotypes
            const genotypes = this.riskVariants.map(v => {
                // Use Hardy-Weinberg to generate genotypes
                const p = v.maf;
                const r = Math.random();
                if (r < (1-p)**2) return 0;
                if (r < (1-p)**2 + 2*p*(1-p)) return 1;
                return 2;
            });
            
            // Calculate simple PRS (sum of weight * dosage)
            const prs = genotypes.reduce((sum, g, idx) => 
                sum + g * this.riskVariants[idx].weight, 0
            );
            
            scores.push(prs);
        }
        
        // Calculate statistics
        scores.sort((a, b) => a - b);
        
        this.populationStats = {
            scores,
            min: scores[0],
            max: scores[scores.length - 1],
            mean: scores.reduce((a, b) => a + b) / scores.length,
            median: scores[Math.floor(scores.length / 2)],
            percentiles: {}
        };
        
        // Pre-calculate percentile thresholds
        [1, 5, 10, 25, 50, 75, 90, 95, 99].forEach(p => {
            const idx = Math.floor(p / 100 * scores.length);
            this.populationStats.percentiles[p] = scores[idx];
        });
        
        // Calculate standard deviation
        const variance = scores.reduce((sum, s) => 
            sum + Math.pow(s - this.populationStats.mean, 2), 0
        ) / scores.length;
        this.populationStats.sd = Math.sqrt(variance);
    }

    /**
     * Calculate PRS from genotypes
     */
    async calculatePRS(genotypes) {
        if (!this.isModelLoaded) {
            throw new Error('Model not loaded');
        }
        
        if (genotypes.length !== this.options.numVariants) {
            throw new Error(`Expected ${this.options.numVariants} genotypes, got ${genotypes.length}`);
        }
        
        // Method 1: Simple weighted sum (traditional PRS)
        const simplePRS = genotypes.reduce((sum, g, idx) => 
            sum + g * this.riskVariants[idx].weight, 0
        );
        
        // Method 2: ML model prediction
        const inputTensor = tf.tensor2d([genotypes], [1, this.options.numVariants]);
        const prediction = this.model.predict(inputTensor);
        const mlPRS = (await prediction.data())[0];
        
        // Cleanup tensors
        inputTensor.dispose();
        prediction.dispose();
        
        // Combine both methods (weighted average)
        const combinedPRS = simplePRS * 0.7 + mlPRS * 0.3;
        
        // Calculate percentile
        const percentile = this.calculatePercentile(combinedPRS);
        
        // Calculate individual variant contributions
        const contributions = this.calculateContributions(genotypes);
        
        return {
            score: combinedPRS,
            standardizedScore: (combinedPRS - this.populationStats.mean) / this.populationStats.sd,
            percentile,
            riskCategory: this.getRiskCategory(percentile),
            contributions,
            simplePRS,
            mlPRS
        };
    }

    /**
     * Calculate percentile for a given PRS
     */
    calculatePercentile(score) {
        const scores = this.populationStats.scores;
        let count = 0;
        
        for (let i = 0; i < scores.length; i++) {
            if (scores[i] < score) count++;
            else break;
        }
        
        return (count / scores.length) * 100;
    }

    /**
     * Get risk category based on percentile
     */
    getRiskCategory(percentile) {
        if (percentile >= 95) return { label: 'Very High Risk', color: '#ef4444', level: 5 };
        if (percentile >= 80) return { label: 'High Risk', color: '#f97316', level: 4 };
        if (percentile >= 50) return { label: 'Average Risk', color: '#eab308', level: 3 };
        if (percentile >= 20) return { label: 'Low Risk', color: '#22c55e', level: 2 };
        return { label: 'Very Low Risk', color: '#10b981', level: 1 };
    }

    /**
     * Calculate individual variant contributions
     */
    calculateContributions(genotypes) {
        return this.riskVariants.map((variant, idx) => {
            const genotype = genotypes[idx];
            const contribution = genotype * variant.weight;
            
            return {
                id: variant.id,
                gene: variant.gene,
                genotype,
                weight: variant.weight,
                contribution,
                percentOfTotal: 0 // Will be calculated below
            };
        }).sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution))
        .map((c, idx, arr) => {
            const totalContribution = arr.reduce((sum, x) => sum + Math.abs(x.contribution), 0);
            c.percentOfTotal = totalContribution > 0 ? 
                (Math.abs(c.contribution) / totalContribution) * 100 : 0;
            return c;
        });
    }

    /**
     * Simulate random genotypes based on population
     */
    simulateGenotypes(population = 'EUR') {
        const popFreq = this.populationFrequencies[population] || this.populationFrequencies.EUR;
        
        return this.riskVariants.map(variant => {
            // Adjust MAF for population
            let maf = variant.maf;
            maf = Math.max(0.01, Math.min(0.49, 
                maf + (Math.random() - 0.5) * popFreq.sdMAF
            ));
            
            // Generate genotype using Hardy-Weinberg
            const p = maf;
            const r = Math.random();
            
            if (r < (1-p)**2) return 0;
            if (r < (1-p)**2 + 2*p*(1-p)) return 1;
            return 2;
        });
    }

    /**
     * Get risk variants for display
     */
    getRiskVariants() {
        return this.riskVariants;
    }

    /**
     * Get population statistics
     */
    getPopulationStats() {
        return this.populationStats;
    }

    /**
     * Check if model is ready
     */
    isReady() {
        return this.isModelLoaded;
    }

    /**
     * Dispose of TensorFlow resources
     */
    dispose() {
        if (this.model) {
            this.model.dispose();
        }
    }
}

// Export
window.PRSCalculator = PRSCalculator;
