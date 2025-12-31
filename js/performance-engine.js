/**
 * GWAS Performance Engine
 * Coordinates WebAssembly, WebGPU, and Web Workers for optimal performance
 *
 * Architecture:
 * ┌─────────────────────────────────────────────────────────────┐
 * │                    Main Thread (UI)                         │
 * │  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐ │
 * │  │ D3.js Axes  │  │ Tooltips    │  │ Event Handlers      │ │
 * │  └─────────────┘  └─────────────┘  └─────────────────────┘ │
 * └──────────────────────────┬──────────────────────────────────┘
 *                            │
 *         ┌──────────────────┼──────────────────┐
 *         │                  │                  │
 *         ▼                  ▼                  ▼
 * ┌───────────────┐  ┌───────────────┐  ┌───────────────┐
 * │  Web Worker   │  │   WebGPU      │  │  WebAssembly  │
 * │ (Background   │  │ (GPU Render)  │  │ (Fast Math)   │
 * │  Processing)  │  │               │  │               │
 * └───────────────┘  └───────────────┘  └───────────────┘
 */

class PerformanceEngine {
    constructor() {
        this.worker = null;
        this.wasmModule = null;
        this.isInitialized = false;
        this.pendingCallbacks = new Map();
        this.callId = 0;

        // Feature detection
        this.features = {
            webWorkers: typeof Worker !== 'undefined',
            webGPU: typeof navigator !== 'undefined' && 'gpu' in navigator,
            webAssembly: typeof WebAssembly !== 'undefined',
            sharedArrayBuffer: typeof SharedArrayBuffer !== 'undefined',
            transferableStreams: typeof ReadableStream !== 'undefined'
        };

        // Performance metrics
        this.metrics = {
            dataGeneration: 0,
            rendering: 0,
            statsCalculation: 0
        };
    }

    /**
     * Initialize all performance features
     */
    async initialize() {
        console.log('Initializing Performance Engine...');
        console.log('Available features:', this.features);

        const promises = [];

        // Initialize Web Worker
        if (this.features.webWorkers) {
            promises.push(this.initWorker());
        }

        // Initialize WebAssembly
        if (this.features.webAssembly) {
            promises.push(this.initWasm());
        }

        await Promise.allSettled(promises);

        this.isInitialized = true;
        console.log('Performance Engine initialized');

        return this.getStatus();
    }

    /**
     * Initialize Web Worker
     */
    async initWorker() {
        return new Promise((resolve, reject) => {
            try {
                this.worker = new Worker('js/workers/gwas-worker.js');

                this.worker.onmessage = (e) => {
                    const { id, success, result, error } = e.data;
                    const callback = this.pendingCallbacks.get(id);

                    if (callback) {
                        if (success) {
                            callback.resolve(result);
                        } else {
                            callback.reject(new Error(error));
                        }
                        this.pendingCallbacks.delete(id);
                    }
                };

                this.worker.onerror = (e) => {
                    console.error('Worker error:', e);
                };

                // Send init message
                this.sendToWorker('init').then(resolve).catch(reject);

            } catch (e) {
                console.warn('Web Worker initialization failed:', e);
                resolve(false);
            }
        });
    }

    /**
     * Initialize WebAssembly module
     */
    async initWasm() {
        try {
            // Check if WASM file exists
            const response = await fetch('wasm/pkg/gwas_wasm_bg.wasm');
            if (!response.ok) {
                console.log('WASM module not found, will use JS fallback');
                return false;
            }

            // Import the generated JS glue code
            const wasmModule = await import('./wasm/pkg/gwas_wasm.js');
            await wasmModule.default();
            this.wasmModule = wasmModule;
            console.log('WebAssembly module loaded');
            return true;

        } catch (e) {
            console.log('WASM not available:', e.message);
            return false;
        }
    }

    /**
     * Send message to worker and get response
     */
    sendToWorker(type, data = {}) {
        return new Promise((resolve, reject) => {
            if (!this.worker) {
                reject(new Error('Worker not initialized'));
                return;
            }

            const id = ++this.callId;
            this.pendingCallbacks.set(id, { resolve, reject });

            // Use transferable objects for typed arrays
            const transfer = [];
            if (data.buffer instanceof ArrayBuffer) {
                transfer.push(data.buffer);
            }

            this.worker.postMessage({ type, data, id }, transfer);
        });
    }

    /**
     * Generate GWAS data (uses Worker + optional WASM)
     */
    async generateGWAS(options = {}) {
        const startTime = performance.now();

        let result;
        if (this.worker) {
            result = await this.sendToWorker('generateGWAS', options);
        } else {
            // Fallback to main thread
            result = this.generateGWASSync(options);
        }

        this.metrics.dataGeneration = performance.now() - startTime;
        console.log(`GWAS generation: ${this.metrics.dataGeneration.toFixed(0)}ms`);

        return result;
    }

    /**
     * Synchronous GWAS generation (fallback)
     */
    generateGWASSync(options) {
        // Uses the existing GWASSimulator class
        if (typeof GWASSimulator !== 'undefined') {
            const simulator = new GWASSimulator(options);
            return simulator.generate();
        }
        throw new Error('GWASSimulator not available');
    }

    /**
     * Calculate QQ plot data in background
     */
    async calculateQQData(pvalues) {
        if (this.worker) {
            return this.sendToWorker('calculateQQ', { pvalues });
        }
        // Fallback
        return this.calculateQQSync(pvalues);
    }

    calculateQQSync(pvalues) {
        const sorted = [...pvalues].sort((a, b) => a - b);
        const n = sorted.length;

        const expected = [], observed = [];
        const step = n > 10000 ? Math.floor(n / 5000) : 1;

        for (let i = 0; i < n; i += step) {
            expected.push(-Math.log10((i + 0.5) / n));
            observed.push(-Math.log10(sorted[i]));
        }

        return { expected, observed };
    }

    /**
     * Filter SNPs by threshold
     */
    async filterSNPs(snps, threshold) {
        if (this.worker && snps.length > 10000) {
            return this.sendToWorker('filterSNPs', { snps, threshold });
        }
        return snps.filter(s => s.pvalue < threshold);
    }

    /**
     * Calculate PRS
     */
    async calculatePRS(genotypes, weights) {
        if (this.worker) {
            return this.sendToWorker('calculatePRS', { genotypes, weights });
        }

        const score = genotypes.reduce((sum, g, i) => sum + g * weights[i], 0);
        return { score, standardizedScore: 0, percentile: 50 };
    }

    /**
     * Get current performance status
     */
    getStatus() {
        return {
            initialized: this.isInitialized,
            features: { ...this.features },
            workerActive: !!this.worker,
            wasmActive: !!this.wasmModule,
            metrics: { ...this.metrics }
        };
    }

    /**
     * Benchmark current configuration
     */
    async benchmark() {
        const results = {};

        // Test data generation
        const start1 = performance.now();
        await this.generateGWAS({ numSNPs: 10000 });
        results.generate10k = performance.now() - start1;

        // Test with more SNPs
        const start2 = performance.now();
        await this.generateGWAS({ numSNPs: 50000 });
        results.generate50k = performance.now() - start2;

        // Test QQ calculation
        const pvals = Array.from({ length: 50000 }, () => Math.random());
        const start3 = performance.now();
        await this.calculateQQData(pvals);
        results.qqCalc50k = performance.now() - start3;

        console.log('Benchmark results:', results);
        return results;
    }

    /**
     * Cleanup resources
     */
    destroy() {
        if (this.worker) {
            this.worker.terminate();
            this.worker = null;
        }
        this.pendingCallbacks.clear();
    }
}

// Typed Array utilities for high-performance data handling
const TypedArrayUtils = {
    /**
     * Convert SNP objects to typed arrays for GPU upload
     */
    snpsToTypedArrays(snps) {
        const n = snps.length;
        const positions = new Float32Array(n * 2);
        const colors = new Float32Array(n * 4);

        for (let i = 0; i < n; i++) {
            const snp = snps[i];
            positions[i * 2] = snp.cumPos;
            positions[i * 2 + 1] = snp.logp;

            // Color encoding
            const isSignificant = snp.pvalue < 5e-8;
            const isSuggestive = snp.pvalue < 1e-5;
            const evenChr = snp.chr % 2 === 0;

            if (isSignificant) {
                colors[i*4] = 0.96; colors[i*4+1] = 0.62; colors[i*4+2] = 0.04; colors[i*4+3] = 1.0;
            } else if (isSuggestive) {
                colors[i*4] = 0.55; colors[i*4+1] = 0.36; colors[i*4+2] = 0.96; colors[i*4+3] = 0.8;
            } else if (evenChr) {
                colors[i*4] = 0.23; colors[i*4+1] = 0.51; colors[i*4+2] = 0.96; colors[i*4+3] = 0.7;
            } else {
                colors[i*4] = 0.02; colors[i*4+1] = 0.71; colors[i*4+2] = 0.83; colors[i*4+3] = 0.7;
            }
        }

        return { positions, colors };
    },

    /**
     * Downsample data for initial fast render
     */
    downsample(snps, targetCount = 10000) {
        if (snps.length <= targetCount) return snps;

        const step = Math.ceil(snps.length / targetCount);
        const sampled = [];

        // Always include significant SNPs
        const significant = snps.filter(s => s.pvalue < 5e-8);
        sampled.push(...significant);

        // Sample remaining
        for (let i = 0; i < snps.length; i += step) {
            if (!snps[i].isSignificant) {
                sampled.push(snps[i]);
            }
        }

        return sampled;
    }
};

// Performance monitoring utilities
const PerfMonitor = {
    marks: {},

    mark(name) {
        this.marks[name] = performance.now();
    },

    measure(name, startMark) {
        const duration = performance.now() - (this.marks[startMark] || 0);
        console.log(`[Perf] ${name}: ${duration.toFixed(2)}ms`);
        return duration;
    },

    async time(name, fn) {
        const start = performance.now();
        const result = await fn();
        console.log(`[Perf] ${name}: ${(performance.now() - start).toFixed(2)}ms`);
        return result;
    }
};

// Export
window.PerformanceEngine = PerformanceEngine;
window.TypedArrayUtils = TypedArrayUtils;
window.PerfMonitor = PerfMonitor;

// Auto-initialize if DOM is ready
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => {
        window.perfEngine = new PerformanceEngine();
        window.perfEngine.initialize();
    });
} else {
    window.perfEngine = new PerformanceEngine();
    window.perfEngine.initialize();
}
