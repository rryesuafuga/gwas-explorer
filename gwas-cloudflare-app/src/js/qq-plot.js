/**
 * QQ Plot - D3.js Implementation
 * Quantile-Quantile plot for assessing GWAS quality
 */

class QQPlot {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = d3.select(`#${containerId}`);
        
        this.options = {
            margin: { top: 30, right: 30, bottom: 60, left: 70 },
            pointRadius: 4,
            pointOpacity: 0.6,
            showCI: true,
            ciLevel: 0.95,
            animationDuration: 800,
            maxPoints: 10000, // Downsample for performance
            ...options
        };
        
        this.data = null;
        this.svg = null;
        
        this.init();
    }

    /**
     * Initialize the plot
     */
    init() {
        this.container.html('');
        
        const rect = this.container.node().getBoundingClientRect();
        this.width = rect.width || 500;
        this.height = rect.height || 450;
        
        this.innerWidth = this.width - this.options.margin.left - this.options.margin.right;
        this.innerHeight = this.height - this.options.margin.top - this.options.margin.bottom;
        
        this.svg = this.container
            .append('svg')
            .attr('width', this.width)
            .attr('height', this.height)
            .attr('class', 'qq-svg');
        
        this.mainGroup = this.svg
            .append('g')
            .attr('transform', `translate(${this.options.margin.left}, ${this.options.margin.top})`);
        
        // Create groups
        this.ciGroup = this.mainGroup.append('g').attr('class', 'ci');
        this.diagonalGroup = this.mainGroup.append('g').attr('class', 'diagonal');
        this.axisGroup = this.mainGroup.append('g').attr('class', 'axes');
        this.pointsGroup = this.mainGroup.append('g').attr('class', 'points');
        
        // Handle resize
        this.resizeObserver = new ResizeObserver(() => this.resize());
        this.resizeObserver.observe(this.container.node());
    }

    /**
     * Render the QQ plot with data
     */
    render(data) {
        if (!data || !data.snps || data.snps.length === 0) {
            this.showEmpty();
            return;
        }
        
        // Calculate QQ values
        this.data = this.calculateQQValues(data.snps);
        
        // Create scales
        this.createScales();
        
        // Render components
        this.renderAxes();
        this.renderDiagonal();
        if (this.options.showCI) {
            this.renderConfidenceInterval();
        }
        this.renderPoints();
        
        // Return statistics
        return this.calculateStatistics();
    }

    /**
     * Calculate expected vs observed -log10(p) values
     */
    calculateQQValues(snps) {
        // Sort p-values
        const pvalues = snps.map(s => s.pvalue).sort((a, b) => a - b);
        const n = pvalues.length;
        
        // Calculate expected values under null hypothesis
        const qqData = pvalues.map((pval, i) => {
            // Expected p-value under uniform distribution
            const expectedP = (i + 0.5) / n;
            return {
                observed: -Math.log10(pval),
                expected: -Math.log10(expectedP),
                pvalue: pval
            };
        });
        
        // Downsample if too many points
        if (qqData.length > this.options.maxPoints) {
            return this.downsample(qqData);
        }
        
        return qqData;
    }

    /**
     * Downsample data for performance while preserving extremes
     */
    downsample(data) {
        const n = data.length;
        const target = this.options.maxPoints;
        
        // Always keep the most significant points
        const topN = Math.min(500, Math.floor(target * 0.1));
        const top = data.slice(-topN);
        
        // Sample from the rest
        const rest = data.slice(0, -topN);
        const sampleRate = (target - topN) / rest.length;
        const sampled = rest.filter(() => Math.random() < sampleRate);
        
        return [...sampled, ...top].sort((a, b) => a.expected - b.expected);
    }

    /**
     * Create scales
     */
    createScales() {
        const maxExpected = d3.max(this.data, d => d.expected);
        const maxObserved = d3.max(this.data, d => d.observed);
        const maxVal = Math.max(maxExpected, maxObserved);
        
        this.xScale = d3.scaleLinear()
            .domain([0, maxVal * 1.05])
            .range([0, this.innerWidth]);
        
        this.yScale = d3.scaleLinear()
            .domain([0, maxVal * 1.05])
            .range([this.innerHeight, 0]);
    }

    /**
     * Render axes
     */
    renderAxes() {
        this.axisGroup.selectAll('*').remove();
        
        // X axis
        const xAxis = d3.axisBottom(this.xScale).ticks(6);
        this.axisGroup.append('g')
            .attr('class', 'x-axis qq-axis')
            .attr('transform', `translate(0, ${this.innerHeight})`)
            .call(xAxis);
        
        // X label
        this.axisGroup.append('text')
            .attr('class', 'axis-label')
            .attr('x', this.innerWidth / 2)
            .attr('y', this.innerHeight + 45)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text('Expected -log₁₀(p)');
        
        // Y axis
        const yAxis = d3.axisLeft(this.yScale).ticks(6);
        this.axisGroup.append('g')
            .attr('class', 'y-axis qq-axis')
            .call(yAxis);
        
        // Y label
        this.axisGroup.append('text')
            .attr('class', 'axis-label')
            .attr('transform', 'rotate(-90)')
            .attr('x', -this.innerHeight / 2)
            .attr('y', -50)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text('Observed -log₁₀(p)');
    }

    /**
     * Render diagonal reference line
     */
    renderDiagonal() {
        this.diagonalGroup.selectAll('*').remove();
        
        const maxVal = this.xScale.domain()[1];
        
        this.diagonalGroup.append('line')
            .attr('class', 'qq-diagonal')
            .attr('x1', this.xScale(0))
            .attr('y1', this.yScale(0))
            .attr('x2', this.xScale(maxVal))
            .attr('y2', this.yScale(maxVal))
            .attr('stroke', '#ef4444')
            .attr('stroke-width', 1.5)
            .attr('stroke-dasharray', '4,2');
    }

    /**
     * Render confidence interval band
     */
    renderConfidenceInterval() {
        this.ciGroup.selectAll('*').remove();
        
        const n = this.data.length;
        const alpha = 1 - this.options.ciLevel;
        
        // Calculate CI bounds using beta distribution
        const ciData = [];
        const numPoints = 100;
        
        for (let i = 0; i < numPoints; i++) {
            const p = (i + 0.5) / numPoints;
            const expected = -Math.log10(p);
            
            if (expected > this.xScale.domain()[1]) continue;
            
            // Beta distribution quantiles for CI
            const rank = Math.floor(p * n);
            const lower = this.betaQuantile(alpha / 2, rank + 1, n - rank);
            const upper = this.betaQuantile(1 - alpha / 2, rank + 1, n - rank);
            
            ciData.push({
                expected,
                lower: -Math.log10(upper),
                upper: -Math.log10(lower)
            });
        }
        
        // Create area generator
        const area = d3.area()
            .x(d => this.xScale(d.expected))
            .y0(d => this.yScale(Math.max(0, d.lower)))
            .y1(d => this.yScale(d.upper))
            .curve(d3.curveMonotoneX);
        
        this.ciGroup.append('path')
            .attr('class', 'qq-ci')
            .attr('d', area(ciData))
            .attr('fill', 'rgba(59, 130, 246, 0.15)')
            .attr('stroke', 'none');
    }

    /**
     * Beta distribution quantile approximation
     */
    betaQuantile(p, a, b) {
        // Simple approximation using normal approximation for large a, b
        if (a < 1 || b < 1) return p;
        
        const mean = a / (a + b);
        const variance = (a * b) / ((a + b) ** 2 * (a + b + 1));
        const std = Math.sqrt(variance);
        
        // Normal approximation
        const z = this.inverseNormalCDF(p);
        let result = mean + z * std;
        
        // Clamp to valid range
        return Math.max(0.0001, Math.min(0.9999, result));
    }

    /**
     * Inverse normal CDF
     */
    inverseNormalCDF(p) {
        if (p <= 0) return -Infinity;
        if (p >= 1) return Infinity;
        
        // Rational approximation
        const a = [
            -3.969683028665376e+01, 2.209460984245205e+02,
            -2.759285104469687e+02, 1.383577518672690e+02,
            -3.066479806614716e+01, 2.506628277459239e+00
        ];
        const b = [
            -5.447609879822406e+01, 1.615858368580409e+02,
            -1.556989798598866e+02, 6.680131188771972e+01,
            -1.328068155288572e+01
        ];
        
        const pLow = 0.02425;
        const pHigh = 1 - pLow;
        
        if (p < pLow) {
            const q = Math.sqrt(-2 * Math.log(p));
            return (((((a[0]*q+a[1])*q+a[2])*q+a[3])*q+a[4])*q+a[5]) /
                   ((((b[0]*q+b[1])*q+b[2])*q+b[3])*q+1);
        } else if (p <= pHigh) {
            const q = p - 0.5;
            const r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
        } else {
            const q = Math.sqrt(-2 * Math.log(1 - p));
            return -(((((a[0]*q+a[1])*q+a[2])*q+a[3])*q+a[4])*q+a[5]) /
                    ((((b[0]*q+b[1])*q+b[2])*q+b[3])*q+1);
        }
    }

    /**
     * Render data points
     */
    renderPoints() {
        this.pointsGroup.selectAll('*').remove();
        
        const points = this.pointsGroup.selectAll('.qq-point')
            .data(this.data)
            .enter()
            .append('circle')
            .attr('class', 'qq-point')
            .attr('cx', d => this.xScale(d.expected))
            .attr('cy', this.yScale(0))
            .attr('r', this.options.pointRadius)
            .attr('fill', d => this.getPointColor(d))
            .attr('opacity', 0);
        
        // Animate points
        points.transition()
            .duration(this.options.animationDuration)
            .delay((d, i) => Math.min(i * 0.05, 300))
            .attr('cy', d => this.yScale(d.observed))
            .attr('opacity', this.options.pointOpacity);
    }

    /**
     * Get point color based on deviation from expected
     */
    getPointColor(d) {
        const deviation = d.observed - d.expected;
        
        if (deviation > 2) {
            return '#f59e0b'; // Significant deviation (inflated)
        } else if (deviation > 1) {
            return '#8b5cf6'; // Moderate deviation
        } else if (deviation < -1) {
            return '#10b981'; // Deflated
        }
        return '#3b82f6'; // Normal
    }

    /**
     * Calculate QC statistics
     */
    calculateStatistics() {
        if (!this.data || this.data.length === 0) return null;
        
        const n = this.data.length;
        
        // Calculate lambda (genomic inflation factor)
        const medianIndex = Math.floor(n / 2);
        const sortedByExpected = [...this.data].sort((a, b) => a.expected - b.expected);
        const medianObserved = sortedByExpected[medianIndex].observed;
        const medianExpected = sortedByExpected[medianIndex].expected;
        
        // Convert -log10(p) back to chi-squared
        const observedChi2 = this.pvalueToChiSquared(Math.pow(10, -medianObserved));
        const expectedChi2 = this.pvalueToChiSquared(Math.pow(10, -medianExpected));
        const lambda = observedChi2 / expectedChi2;
        
        // Lambda scaled to 1000 cases and 1000 controls
        const lambda1000 = 1 + (lambda - 1) * ((1/1000 + 1/1000) / (1/n + 1/n));
        
        // Mean chi-squared
        const chiSquaredValues = this.data.map(d => 
            this.pvalueToChiSquared(Math.pow(10, -d.observed))
        );
        const meanChi2 = chiSquaredValues.reduce((a, b) => a + b, 0) / n;
        
        // Tail deviation (proportion of points above diagonal in tail)
        const tail = this.data.filter(d => d.expected > 2);
        const tailDeviation = tail.length > 0 ? 
            tail.filter(d => d.observed > d.expected).length / tail.length : 0;
        
        return {
            lambda,
            lambda1000,
            meanChi2,
            tailDeviation,
            nSNPs: n
        };
    }

    /**
     * Convert p-value to chi-squared
     */
    pvalueToChiSquared(p) {
        const z = this.inverseNormalCDF(1 - p/2);
        return z * z;
    }

    /**
     * Get interpretation text
     */
    getInterpretation(stats) {
        if (!stats) return 'No data available for interpretation.';
        
        const parts = [];
        
        // Lambda interpretation
        if (stats.lambda < 1.02) {
            parts.push('Genomic inflation factor (λ) indicates minimal population stratification.');
        } else if (stats.lambda < 1.1) {
            parts.push('Moderate genomic inflation detected (λ = ' + stats.lambda.toFixed(3) + '). Consider adjusting for population structure.');
        } else {
            parts.push('⚠️ High genomic inflation (λ = ' + stats.lambda.toFixed(3) + '). Strong population stratification or other systematic bias present.');
        }
        
        // Tail deviation
        if (stats.tailDeviation > 0.7) {
            parts.push('Strong signal enrichment in the tail suggests true associations.');
        } else if (stats.tailDeviation > 0.5) {
            parts.push('Moderate signal enrichment observed.');
        }
        
        // Mean chi-squared
        if (stats.meanChi2 > 1.5) {
            parts.push('Elevated mean χ² (' + stats.meanChi2.toFixed(2) + ') indicates polygenicity or inflation.');
        }
        
        return parts.join(' ');
    }

    /**
     * Show empty state
     */
    showEmpty() {
        this.pointsGroup.selectAll('*').remove();
        this.ciGroup.selectAll('*').remove();
        this.diagonalGroup.selectAll('*').remove();
        this.axisGroup.selectAll('*').remove();
        
        this.mainGroup.append('text')
            .attr('x', this.innerWidth / 2)
            .attr('y', this.innerHeight / 2)
            .attr('text-anchor', 'middle')
            .attr('fill', '#64748b')
            .attr('font-size', '14px')
            .text('No data loaded');
    }

    /**
     * Handle resize
     */
    resize() {
        const rect = this.container.node().getBoundingClientRect();
        if (rect.width === this.width && rect.height === this.height) return;
        
        this.width = rect.width;
        this.height = rect.height;
        this.innerWidth = this.width - this.options.margin.left - this.options.margin.right;
        this.innerHeight = this.height - this.options.margin.top - this.options.margin.bottom;
        
        this.svg
            .attr('width', this.width)
            .attr('height', this.height);
        
        if (this.data) {
            this.createScales();
            this.renderAxes();
            this.renderDiagonal();
            if (this.options.showCI) {
                this.renderConfidenceInterval();
            }
            this.renderPoints();
        }
    }

    /**
     * Cleanup
     */
    destroy() {
        if (this.resizeObserver) {
            this.resizeObserver.disconnect();
        }
        this.container.html('');
    }
}

// Export
window.QQPlot = QQPlot;
