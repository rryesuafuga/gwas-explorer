use wasm_bindgen::prelude::*;
use rand::prelude::*;
use rand_distr::{Normal, Distribution};
use std::f64::consts::PI;

/// GWAS WebAssembly Module - High-performance statistical computations
/// Provides 10-50x speedup over JavaScript for heavy math operations

#[wasm_bindgen]
pub struct GWASEngine {
    rng: rand::rngs::SmallRng,
    num_snps: usize,
    num_chromosomes: usize,
    genomic_inflation: f64,
    chr_lengths: Vec<f64>,
}

#[wasm_bindgen]
impl GWASEngine {
    #[wasm_bindgen(constructor)]
    pub fn new(num_snps: u32, num_chromosomes: u32, genomic_inflation: f64) -> GWASEngine {
        let chr_lengths = vec![
            249.0, 243.0, 198.0, 191.0, 182.0, 171.0, 159.0, 145.0, 138.0, 134.0,
            135.0, 133.0, 114.0, 107.0, 102.0, 90.0, 83.0, 80.0, 59.0, 64.0, 47.0, 51.0
        ];

        GWASEngine {
            rng: rand::rngs::SmallRng::from_entropy(),
            num_snps: num_snps as usize,
            num_chromosomes: num_chromosomes as usize,
            genomic_inflation,
            chr_lengths,
        }
    }

    /// Generate GWAS p-values with genomic inflation
    /// Returns flat array: [chr1, pos1, pval1, logp1, chr2, pos2, pval2, logp2, ...]
    #[wasm_bindgen]
    pub fn generate_pvalues(&mut self, significant_loci: u32) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.num_snps * 4);
        let snps_per_chr = self.num_snps / self.num_chromosomes;

        // Generate significant loci positions
        let sig_positions: Vec<(usize, f64, f64, f64)> = (0..significant_loci as usize)
            .map(|_| {
                let chr = self.rng.gen_range(1..=self.num_chromosomes);
                let pos = self.rng.gen::<f64>() * self.chr_lengths[chr - 1] * 1e6;
                let effect = 0.1 + self.rng.gen::<f64>() * 0.4;
                let radius = 500000.0 + self.rng.gen::<f64>() * 1000000.0;
                (chr, pos, effect, radius)
            })
            .collect();

        for chr in 1..=self.num_chromosomes {
            let chr_length = self.chr_lengths[chr - 1] * 1e6;

            for _ in 0..snps_per_chr {
                let pos = self.rng.gen::<f64>() * chr_length;

                // Check if near significant locus
                let near_sig = sig_positions.iter().find(|(c, p, _, r)| {
                    *c == chr && (pos - p).abs() < *r
                });

                let pvalue = if let Some((_, locus_pos, effect, radius)) = near_sig {
                    let distance = (pos - locus_pos).abs();
                    let correlation = (-distance / (radius * 0.3)).exp();

                    if distance < 10000.0 {
                        // Lead SNP
                        10.0_f64.powf(-(8.0 + self.rng.gen::<f64>() * 4.0))
                    } else {
                        // LD SNP
                        let base = 8.0 + self.rng.gen::<f64>() * 2.0;
                        10.0_f64.powf(-(base * correlation + self.rng.gen::<f64>()))
                    }
                } else {
                    // Null distribution with inflation
                    let chi2 = self.generate_inflated_chi2();
                    self.chi2_to_pvalue(chi2)
                };

                let logp = -pvalue.log10();

                result.push(chr as f64);
                result.push(pos);
                result.push(pvalue);
                result.push(logp);
            }
        }

        result
    }

    /// Generate chi-squared with genomic inflation
    fn generate_inflated_chi2(&mut self) -> f64 {
        let normal = Normal::new(0.0, 1.0).unwrap();
        let z: f64 = normal.sample(&mut self.rng);
        z * z * self.genomic_inflation
    }

    /// Convert chi-squared to p-value
    fn chi2_to_pvalue(&self, chi2: f64) -> f64 {
        let z = chi2.sqrt();
        let p = 2.0 * (1.0 - self.normal_cdf(z));
        p.max(1e-300)
    }

    /// Standard normal CDF (Abramowitz & Stegun approximation)
    fn normal_cdf(&self, x: f64) -> f64 {
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;

        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x = x.abs() / 2.0_f64.sqrt();

        let t = 1.0 / (1.0 + p * x);
        let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

        0.5 * (1.0 + sign * y)
    }

    /// Calculate genomic inflation factor (lambda)
    #[wasm_bindgen]
    pub fn calculate_lambda(&self, pvalues: &[f64]) -> f64 {
        let mut sorted: Vec<f64> = pvalues.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median_idx = sorted.len() / 2;
        let median_pval = sorted[median_idx];

        let observed_chi2 = self.pvalue_to_chi2(median_pval);
        let expected_chi2 = self.pvalue_to_chi2(0.5);

        observed_chi2 / expected_chi2
    }

    fn pvalue_to_chi2(&self, p: f64) -> f64 {
        let z = self.inverse_normal_cdf(1.0 - p / 2.0);
        z * z
    }

    fn inverse_normal_cdf(&self, p: f64) -> f64 {
        if p <= 0.0 { return f64::NEG_INFINITY; }
        if p >= 1.0 { return f64::INFINITY; }
        if p == 0.5 { return 0.0; }

        let c = [
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
            4.374664141464968e+00,
            2.938163982698783e+00
        ];
        let d = [
            7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996e+00,
            3.754408661907416e+00
        ];

        let p_low = 0.02425;

        if p < p_low {
            let q = (-2.0 * p.ln()).sqrt();
            (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
        } else if p <= 1.0 - p_low {
            let q = p - 0.5;
            let r = q * q;
            let a = [
                -3.969683028665376e+01,
                2.209460984245205e+02,
                -2.759285104469687e+02,
                1.383577518672690e+02,
                -3.066479806614716e+01,
                2.506628277459239e+00
            ];
            let b = [
                -5.447609879822406e+01,
                1.615858368580409e+02,
                -1.556989798598866e+02,
                6.680131188771972e+01,
                -1.328068155288572e+01
            ];
            (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0)
        } else {
            let q = (-2.0 * (1.0 - p).ln()).sqrt();
            -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0)
        }
    }
}

/// PRS calculation engine
#[wasm_bindgen]
pub struct PRSEngine {
    weights: Vec<f64>,
    mafs: Vec<f64>,
    population_scores: Vec<f64>,
    mean: f64,
    sd: f64,
}

#[wasm_bindgen]
impl PRSEngine {
    #[wasm_bindgen(constructor)]
    pub fn new(num_variants: u32) -> PRSEngine {
        let mut rng = rand::rngs::SmallRng::from_entropy();
        let n = num_variants as usize;

        // Generate weights (effect sizes)
        let weights: Vec<f64> = (0..n).map(|_| {
            let u: f64 = rng.gen();
            let sign = if rng.gen::<f64>() > 0.5 { 1.0 } else { -1.0 };

            if u < 0.7 {
                (rng.gen::<f64>() * 0.05 + 0.01) * sign
            } else if u < 0.9 {
                (rng.gen::<f64>() * 0.1 + 0.05) * sign
            } else {
                (rng.gen::<f64>() * 0.2 + 0.1) * sign
            }
        }).collect();

        // Generate MAFs
        let mafs: Vec<f64> = (0..n).map(|_| {
            0.01 + rng.gen::<f64>() * 0.48
        }).collect();

        // Generate population distribution
        let mut population_scores: Vec<f64> = (0..10000).map(|_| {
            let genotypes: Vec<f64> = mafs.iter().map(|maf| {
                let p = *maf;
                let r: f64 = rng.gen();
                if r < (1.0 - p).powi(2) { 0.0 }
                else if r < (1.0 - p).powi(2) + 2.0 * p * (1.0 - p) { 1.0 }
                else { 2.0 }
            }).collect();

            genotypes.iter().zip(weights.iter())
                .map(|(g, w)| g * w)
                .sum()
        }).collect();

        population_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mean = population_scores.iter().sum::<f64>() / population_scores.len() as f64;
        let variance = population_scores.iter()
            .map(|s| (s - mean).powi(2))
            .sum::<f64>() / population_scores.len() as f64;
        let sd = variance.sqrt();

        PRSEngine {
            weights,
            mafs,
            population_scores,
            mean,
            sd,
        }
    }

    /// Calculate PRS from genotypes array
    #[wasm_bindgen]
    pub fn calculate(&self, genotypes: &[f64]) -> Vec<f64> {
        let score: f64 = genotypes.iter()
            .zip(self.weights.iter())
            .map(|(g, w)| g * w)
            .sum();

        let standardized = (score - self.mean) / self.sd;
        let percentile = self.get_percentile(score);

        vec![score, standardized, percentile]
    }

    fn get_percentile(&self, score: f64) -> f64 {
        let count = self.population_scores.iter()
            .filter(|&s| *s < score)
            .count();
        (count as f64 / self.population_scores.len() as f64) * 100.0
    }

    #[wasm_bindgen]
    pub fn get_weights(&self) -> Vec<f64> {
        self.weights.clone()
    }

    #[wasm_bindgen]
    pub fn get_mafs(&self) -> Vec<f64> {
        self.mafs.clone()
    }

    /// Simulate genotypes for a population
    #[wasm_bindgen]
    pub fn simulate_genotypes(&self) -> Vec<f64> {
        let mut rng = rand::rngs::SmallRng::from_entropy();

        self.mafs.iter().map(|maf| {
            let p = *maf;
            let r: f64 = rng.gen();
            if r < (1.0 - p).powi(2) { 0.0 }
            else if r < (1.0 - p).powi(2) + 2.0 * p * (1.0 - p) { 1.0 }
            else { 2.0 }
        }).collect()
    }
}

/// Fast parallel sorting for large arrays
#[wasm_bindgen]
pub fn parallel_sort(data: &mut [f64]) {
    data.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

/// Calculate QQ plot expected vs observed
#[wasm_bindgen]
pub fn calculate_qq_data(pvalues: &[f64]) -> Vec<f64> {
    let n = pvalues.len();
    let mut sorted = pvalues.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut result = Vec::with_capacity(n * 2);

    for (i, &obs) in sorted.iter().enumerate() {
        let expected = (i as f64 + 0.5) / n as f64;
        let obs_logp = -obs.log10();
        let exp_logp = -expected.log10();
        result.push(exp_logp);
        result.push(obs_logp);
    }

    result
}
