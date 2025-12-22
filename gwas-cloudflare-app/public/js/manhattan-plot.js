/**
 * Manhattan Plot - D3.js Implementation
 * Interactive genome-wide visualization of association p-values
 */

class ManhattanPlot {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = d3.select(`#${containerId}`);
        
        // Default options
        this.options = {
            margin: { top: 30, right: 30, bottom: 60, left: 70 },
            pointRadius: 3,
            pointOpacity: 0.7,
            significanceThreshold: 5e-8,
            suggestiveThreshold: 1e-5,
            colorScheme: 'alternating',
            animationDuration: 500,
            ...options
        };
        
        // Color schemes
        this.colorSchemes = {
            alternating: (chr) => chr % 2 === 0 ? '#3b82f6' : '#06b6d4',
            gradient: d3.scaleSequential(d3.interpolateViridis),
            significance: (snp) => {
                if (snp.pvalue < 5e-8) return '#f59e0b';
                if (snp.pvalue < 1e-5) return '#8b5cf6';
                return snp.chr % 2 === 0 ? '#3b82f6' : '#06b6d4';
            }
        };
        
        this.data = null;
        this.svg = null;
        this.zoom = null;
        this.currentTransform = d3.zoomIdentity;
        
        // Callbacks
        this.onSNPClick = null;
        this.onSNPHover = null;
        
        this.init();
    }

    /**
     * Initialize the plot
     */
    init() {
        // Clear existing content
        this.container.html('');
        
        // Get container dimensions
        const rect = this.container.node().getBoundingClientRect();
        this.width = rect.width || 800;
        this.height = rect.height || 450;
        
        this.innerWidth = this.width - this.options.margin.left - this.options.margin.right;
        this.innerHeight = this.height - this.options.margin.top - this.options.margin.bottom;
        
        // Create SVG
        this.svg = this.container
            .append('svg')
            .attr('width', this.width)
            .attr('height', this.height)
            .attr('class', 'manhattan-svg');
        
        // Add gradient definitions
        this.addDefs();
        
        // Create main group
        this.mainGroup = this.svg
            .append('g')
            .attr('transform', `translate(${this.options.margin.left}, ${this.options.margin.top})`);
        
        // Create clip path for zooming
        this.svg.append('defs')
            .append('clipPath')
            .attr('id', 'manhattan-clip')
            .append('rect')
            .attr('width', this.innerWidth)
            .attr('height', this.innerHeight);
        
        // Create groups for different elements
        this.axisGroup = this.mainGroup.append('g').attr('class', 'axes');
        this.thresholdGroup = this.mainGroup.append('g').attr('class', 'thresholds');
        this.pointsGroup = this.mainGroup.append('g')
            .attr('class', 'points')
            .attr('clip-path', 'url(#manhattan-clip)');
        this.labelsGroup = this.mainGroup.append('g').attr('class', 'labels');
        
        // Initialize zoom behavior
        this.initZoom();
        
        // Handle resize
        this.resizeObserver = new ResizeObserver(() => this.resize());
        this.resizeObserver.observe(this.container.node());
    }

    /**
     * Add SVG definitions (gradients, etc.)
     */
    addDefs() {
        const defs = this.svg.append('defs');
        
        // Gradient for significant points
        const gradient = defs.append('radialGradient')
            .attr('id', 'significant-gradient');
        
        gradient.append('stop')
            .attr('offset', '0%')
            .attr('stop-color', '#fbbf24');
        
        gradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', '#f59e0b');
        
        // Glow filter
        const filter = defs.append('filter')
            .attr('id', 'glow')
            .attr('x', '-50%')
            .attr('y', '-50%')
            .attr('width', '200%')
            .attr('height', '200%');
        
        filter.append('feGaussianBlur')
            .attr('stdDeviation', '2')
            .attr('result', 'coloredBlur');
        
        const feMerge = filter.append('feMerge');
        feMerge.append('feMergeNode').attr('in', 'coloredBlur');
        feMerge.append('feMergeNode').attr('in', 'SourceGraphic');
    }

    /**
     * Initialize zoom behavior
     */
    initZoom() {
        this.zoom = d3.zoom()
            .scaleExtent([1, 50])
            .translateExtent([[0, 0], [this.innerWidth, this.innerHeight]])
            .extent([[0, 0], [this.innerWidth, this.innerHeight]])
            .on('zoom', (event) => this.handleZoom(event));
        
        this.svg.call(this.zoom);
    }

    /**
     * Handle zoom events
     */
    handleZoom(event) {
        this.currentTransform = event.transform;
        
        if (!this.data) return;
        
        // Update x scale
        const newXScale = this.currentTransform.rescaleX(this.xScale);
        
        // Update x axis
        this.axisGroup.select('.x-axis')
            .call(d3.axisBottom(newXScale)
                .tickValues(this.getChromosomeTicks(newXScale))
                .tickFormat((d, i) => this.chromosomeLabels[i] || ''));
        
        // Update points
        this.pointsGroup.selectAll('.snp-point')
            .attr('cx', d => newXScale(d.cumPos));
        
        // Update threshold lines
        this.thresholdGroup.selectAll('.threshold-line')
            .attr('x1', newXScale.range()[0])
            .attr('x2', newXScale.range()[1]);
    }

    /**
     * Get chromosome tick positions for current zoom level
     */
    getChromosomeTicks(scale) {
        if (!this.data) return [];
        
        const [xMin, xMax] = scale.domain();
        const visibleChrs = this.data.metadata.chromosomeInfo.filter(chr => 
            chr.center >= xMin && chr.center <= xMax
        );
        
        this.chromosomeLabels = visibleChrs.map(chr => chr.chr.toString());
        return visibleChrs.map(chr => chr.center);
    }

    /**
     * Render the plot with data
     */
    render(data) {
        this.data = data;
        
        if (!data || !data.snps || data.snps.length === 0) {
            this.showEmpty();
            return;
        }
        
        // Create scales
        this.createScales();
        
        // Render axes
        this.renderAxes();
        
        // Render threshold lines
        this.renderThresholds();
        
        // Render points
        this.renderPoints();
        
        // Render labels for top hits
        this.renderLabels();
    }

    /**
     * Create x and y scales
     */
    createScales() {
        const snps = this.data.snps;
        
        // X scale: cumulative genomic position
        const xExtent = d3.extent(snps, d => d.cumPos);
        this.xScale = d3.scaleLinear()
            .domain([xExtent[0] - 1e7, xExtent[1] + 1e7])
            .range([0, this.innerWidth]);
        
        // Y scale: -log10(p-value)
        const maxLogP = Math.max(d3.max(snps, d => d.logp), 10);
        this.yScale = d3.scaleLinear()
            .domain([0, maxLogP * 1.1])
            .range([this.innerHeight, 0]);
        
        // Store initial chromosome tick positions
        this.chromosomeLabels = this.data.metadata.chromosomeInfo.map(chr => chr.chr.toString());
    }

    /**
     * Render axes
     */
    renderAxes() {
        this.axisGroup.selectAll('*').remove();
        
        // X axis with chromosome labels
        const xAxis = d3.axisBottom(this.xScale)
            .tickValues(this.data.metadata.chromosomeInfo.map(chr => chr.center))
            .tickFormat((d, i) => this.data.metadata.chromosomeInfo[i]?.chr || '');
        
        this.axisGroup.append('g')
            .attr('class', 'x-axis manhattan-axis')
            .attr('transform', `translate(0, ${this.innerHeight})`)
            .call(xAxis);
        
        // X axis label
        this.axisGroup.append('text')
            .attr('class', 'axis-label')
            .attr('x', this.innerWidth / 2)
            .attr('y', this.innerHeight + 45)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text('Chromosome');
        
        // Y axis
        const yAxis = d3.axisLeft(this.yScale)
            .ticks(8)
            .tickFormat(d => d.toFixed(0));
        
        this.axisGroup.append('g')
            .attr('class', 'y-axis manhattan-axis')
            .call(yAxis);
        
        // Y axis label
        this.axisGroup.append('text')
            .attr('class', 'axis-label')
            .attr('transform', 'rotate(-90)')
            .attr('x', -this.innerHeight / 2)
            .attr('y', -50)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text('-log₁₀(p-value)');
        
        // Add chromosome alternating background
        this.renderChromosomeBackgrounds();
    }

    /**
     * Render alternating chromosome backgrounds
     */
    renderChromosomeBackgrounds() {
        const bgGroup = this.axisGroup.append('g').attr('class', 'chr-backgrounds');
        
        this.data.metadata.chromosomeInfo.forEach((chr, i) => {
            if (i % 2 === 0) {
                bgGroup.append('rect')
                    .attr('x', this.xScale(chr.start))
                    .attr('y', 0)
                    .attr('width', this.xScale(chr.end) - this.xScale(chr.start))
                    .attr('height', this.innerHeight)
                    .attr('fill', 'rgba(255, 255, 255, 0.02)');
            }
        });
    }

    /**
     * Render significance threshold lines
     */
    renderThresholds() {
        this.thresholdGroup.selectAll('*').remove();
        
        // Genome-wide significance (5e-8)
        const sigY = this.yScale(-Math.log10(this.options.significanceThreshold));
        if (sigY > 0 && sigY < this.innerHeight) {
            this.thresholdGroup.append('line')
                .attr('class', 'threshold-line')
                .attr('x1', 0)
                .attr('x2', this.innerWidth)
                .attr('y1', sigY)
                .attr('y2', sigY)
                .attr('stroke', '#f59e0b')
                .attr('stroke-width', 1.5)
                .attr('stroke-dasharray', '6,4');
            
            this.thresholdGroup.append('text')
                .attr('x', this.innerWidth - 5)
                .attr('y', sigY - 5)
                .attr('text-anchor', 'end')
                .attr('fill', '#f59e0b')
                .attr('font-size', '10px')
                .attr('font-family', 'JetBrains Mono, monospace')
                .text('p = 5×10⁻⁸');
        }
        
        // Suggestive threshold (1e-5)
        const sugY = this.yScale(-Math.log10(this.options.suggestiveThreshold));
        if (sugY > 0 && sugY < this.innerHeight) {
            this.thresholdGroup.append('line')
                .attr('class', 'threshold-line')
                .attr('x1', 0)
                .attr('x2', this.innerWidth)
                .attr('y1', sugY)
                .attr('y2', sugY)
                .attr('stroke', '#8b5cf6')
                .attr('stroke-width', 1)
                .attr('stroke-dasharray', '4,4');
            
            this.thresholdGroup.append('text')
                .attr('x', this.innerWidth - 5)
                .attr('y', sugY - 5)
                .attr('text-anchor', 'end')
                .attr('fill', '#8b5cf6')
                .attr('font-size', '10px')
                .attr('font-family', 'JetBrains Mono, monospace')
                .text('p = 1×10⁻⁵');
        }
    }

    /**
     * Render data points
     */
    renderPoints() {
        const self = this;
        const snps = this.data.snps;
        
        // Determine color function based on scheme
        const colorFn = this.getColorFunction();
        
        // Use canvas for large datasets, SVG for smaller ones
        if (snps.length > 100000) {
            this.renderPointsCanvas(snps, colorFn);
        } else {
            this.renderPointsSVG(snps, colorFn);
        }
    }

    /**
     * Render points using SVG (for interactivity)
     */
    renderPointsSVG(snps, colorFn) {
        const self = this;
        
        // Remove existing points
        this.pointsGroup.selectAll('.snp-point').remove();
        
        // Add points with animation
        const points = this.pointsGroup.selectAll('.snp-point')
            .data(snps)
            .enter()
            .append('circle')
            .attr('class', d => `snp-point ${d.isSignificant ? 'significant' : ''}`)
            .attr('cx', d => this.xScale(d.cumPos))
            .attr('cy', this.innerHeight)
            .attr('r', d => d.isSignificant ? this.options.pointRadius + 2 : this.options.pointRadius)
            .attr('fill', d => colorFn(d))
            .attr('opacity', 0)
            .style('filter', d => d.isSignificant ? 'url(#glow)' : 'none');
        
        // Animate points appearing
        points.transition()
            .duration(this.options.animationDuration)
            .delay((d, i) => Math.min(i * 0.01, 500))
            .attr('cy', d => this.yScale(d.logp))
            .attr('opacity', d => d.isSignificant ? 1 : this.options.pointOpacity);
        
        // Add interactivity
        points
            .on('mouseenter', function(event, d) {
                d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('r', self.options.pointRadius + 3)
                    .attr('opacity', 1);
                
                if (self.onSNPHover) {
                    self.onSNPHover(d, event);
                }
            })
            .on('mouseleave', function(event, d) {
                d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('r', d.isSignificant ? self.options.pointRadius + 2 : self.options.pointRadius)
                    .attr('opacity', d.isSignificant ? 1 : self.options.pointOpacity);
                
                if (self.onSNPHover) {
                    self.onSNPHover(null, event);
                }
            })
            .on('click', function(event, d) {
                if (self.onSNPClick) {
                    self.onSNPClick(d, event);
                }
            });
    }

    /**
     * Render points using Canvas (for large datasets)
     */
    renderPointsCanvas(snps, colorFn) {
        // For very large datasets, use a canvas overlay
        // This is more efficient than SVG for >100k points
        
        const canvas = this.pointsGroup.select('canvas');
        if (canvas.empty()) {
            this.pointsGroup.append('foreignObject')
                .attr('width', this.innerWidth)
                .attr('height', this.innerHeight)
                .append('xhtml:canvas')
                .attr('width', this.innerWidth)
                .attr('height', this.innerHeight);
        }
        
        const ctx = this.pointsGroup.select('canvas').node().getContext('2d');
        ctx.clearRect(0, 0, this.innerWidth, this.innerHeight);
        
        snps.forEach(snp => {
            const x = this.xScale(snp.cumPos);
            const y = this.yScale(snp.logp);
            const r = snp.isSignificant ? this.options.pointRadius + 2 : this.options.pointRadius;
            
            ctx.beginPath();
            ctx.arc(x, y, r, 0, 2 * Math.PI);
            ctx.fillStyle = colorFn(snp);
            ctx.globalAlpha = snp.isSignificant ? 1 : this.options.pointOpacity;
            ctx.fill();
        });
    }

    /**
     * Render labels for top hits
     */
    renderLabels() {
        this.labelsGroup.selectAll('*').remove();
        
        // Get top significant SNPs
        const topHits = this.data.snps
            .filter(s => s.isSignificant)
            .sort((a, b) => a.pvalue - b.pvalue)
            .slice(0, 10);
        
        // Add labels with collision avoidance
        const labelPositions = [];
        
        topHits.forEach(snp => {
            const x = this.xScale(snp.cumPos);
            const y = this.yScale(snp.logp);
            
            // Simple collision detection
            const collision = labelPositions.some(pos => 
                Math.abs(pos.x - x) < 60 && Math.abs(pos.y - y) < 20
            );
            
            if (!collision) {
                labelPositions.push({ x, y });
                
                this.labelsGroup.append('text')
                    .attr('x', x)
                    .attr('y', y - 10)
                    .attr('text-anchor', 'middle')
                    .attr('fill', '#f0f4f8')
                    .attr('font-size', '10px')
                    .attr('font-family', 'JetBrains Mono, monospace')
                    .attr('opacity', 0)
                    .text(snp.gene || snp.id)
                    .transition()
                    .delay(this.options.animationDuration)
                    .duration(300)
                    .attr('opacity', 1);
            }
        });
    }

    /**
     * Get color function based on current scheme
     */
    getColorFunction() {
        switch (this.options.colorScheme) {
            case 'gradient':
                const numChrs = this.data.metadata.chromosomeInfo.length;
                return (snp) => {
                    if (snp.isSignificant) return '#f59e0b';
                    return d3.interpolateViridis(snp.chr / numChrs);
                };
            case 'significance':
                return this.colorSchemes.significance;
            case 'alternating':
            default:
                return (snp) => {
                    if (snp.isSignificant) return '#f59e0b';
                    return this.colorSchemes.alternating(snp.chr);
                };
        }
    }

    /**
     * Update color scheme
     */
    setColorScheme(scheme) {
        this.options.colorScheme = scheme;
        if (this.data) {
            const colorFn = this.getColorFunction();
            this.pointsGroup.selectAll('.snp-point')
                .transition()
                .duration(300)
                .attr('fill', d => colorFn(d));
        }
    }

    /**
     * Update significance threshold
     */
    setThreshold(threshold) {
        this.options.significanceThreshold = threshold;
        this.renderThresholds();
    }

    /**
     * Reset zoom to initial state
     */
    resetZoom() {
        this.svg.transition()
            .duration(500)
            .call(this.zoom.transform, d3.zoomIdentity);
    }

    /**
     * Zoom to specific region
     */
    zoomToRegion(chr, start, end) {
        if (!this.data) return;
        
        const chrInfo = this.data.metadata.chromosomeInfo.find(c => c.chr === chr);
        if (!chrInfo) return;
        
        const x0 = chrInfo.start + start;
        const x1 = chrInfo.start + end;
        
        const scale = this.innerWidth / (this.xScale(x1) - this.xScale(x0));
        const translate = -this.xScale(x0) * scale;
        
        this.svg.transition()
            .duration(750)
            .call(this.zoom.transform, d3.zoomIdentity.translate(translate, 0).scale(scale));
    }

    /**
     * Highlight specific SNP
     */
    highlightSNP(snpId) {
        this.pointsGroup.selectAll('.snp-point')
            .classed('highlighted', d => d.id === snpId);
    }

    /**
     * Export plot as SVG
     */
    exportSVG() {
        const svgNode = this.svg.node().cloneNode(true);
        const serializer = new XMLSerializer();
        let svgString = serializer.serializeToString(svgNode);
        
        // Add XML declaration and styling
        svgString = '<?xml version="1.0" encoding="UTF-8"?>\n' + svgString;
        
        const blob = new Blob([svgString], { type: 'image/svg+xml' });
        const url = URL.createObjectURL(blob);
        
        const link = document.createElement('a');
        link.href = url;
        link.download = 'manhattan_plot.svg';
        link.click();
        
        URL.revokeObjectURL(url);
    }

    /**
     * Show empty state
     */
    showEmpty() {
        this.pointsGroup.selectAll('*').remove();
        this.axisGroup.selectAll('*').remove();
        this.thresholdGroup.selectAll('*').remove();
        this.labelsGroup.selectAll('*').remove();
        
        this.mainGroup.append('text')
            .attr('x', this.innerWidth / 2)
            .attr('y', this.innerHeight / 2)
            .attr('text-anchor', 'middle')
            .attr('fill', '#64748b')
            .attr('font-size', '16px')
            .text('No data loaded. Click "Simulate GWAS" or upload your data.');
    }

    /**
     * Handle container resize
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
        
        this.svg.select('#manhattan-clip rect')
            .attr('width', this.innerWidth)
            .attr('height', this.innerHeight);
        
        if (this.data) {
            this.createScales();
            this.renderAxes();
            this.renderThresholds();
            this.renderPoints();
            this.renderLabels();
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
window.ManhattanPlot = ManhattanPlot;
