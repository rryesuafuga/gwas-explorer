/**
 * LocusZoom Integration
 * Regional association plot with gene annotations
 */

class LocusZoomIntegration {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        
        this.options = {
            width: 800,
            height: 500,
            flankingRegion: 500000, // 500kb on each side
            ...options
        };
        
        this.data = null;
        this.plot = null;
        this.currentRegion = null;
        
        // Check if LocusZoom is available
        this.locusZoomAvailable = typeof LocusZoom !== 'undefined';
        
        if (!this.locusZoomAvailable) {
            console.warn('LocusZoom.js not loaded. Using fallback D3 implementation.');
        }
    }

    /**
     * Set GWAS data
     */
    setData(gwasData) {
        this.data = gwasData;
    }

    /**
     * View a specific genomic region
     */
    viewRegion(chr, start, end) {
        if (!this.data) {
            console.error('No GWAS data loaded');
            return;
        }
        
        this.currentRegion = { chr, start, end };
        
        // Filter SNPs in region
        const regionSNPs = this.data.snps.filter(snp => 
            snp.chr === chr && snp.pos >= start && snp.pos <= end
        );
        
        if (regionSNPs.length === 0) {
            this.showEmpty('No SNPs found in this region');
            return;
        }
        
        if (this.locusZoomAvailable) {
            this.renderWithLocusZoom(chr, start, end, regionSNPs);
        } else {
            this.renderFallback(chr, start, end, regionSNPs);
        }
    }

    /**
     * View region around a specific SNP
     */
    viewSNP(snp, flankingRegion = null) {
        const flanking = flankingRegion || this.options.flankingRegion;
        const start = Math.max(0, snp.pos - flanking);
        const end = snp.pos + flanking;
        
        this.viewRegion(snp.chr, start, end);
    }

    /**
     * Render using LocusZoom.js (if available)
     */
    renderWithLocusZoom(chr, start, end, snps) {
        // Clear previous plot
        this.container.innerHTML = '';
        
        // Create data source
        const assocData = snps.map(snp => ({
            variant: `${snp.chr}:${snp.pos}_${snp.ref}/${snp.alt}`,
            position: snp.pos,
            log_pvalue: snp.logp,
            ref_allele: snp.ref
        }));
        
        // Define data sources
        const dataSources = new LocusZoom.DataSources();
        
        // Add association data source
        dataSources.add('assoc', ['AssociationLZ', {
            data: assocData,
            id_field: 'variant',
            position_field: 'position',
            pvalue_field: 'log_pvalue'
        }]);
        
        // Define layout
        const layout = LocusZoom.Layouts.get('standard_association', {
            width: this.options.width,
            height: this.options.height,
            state: {
                chr: chr,
                start: start,
                end: end,
                genome_build: 'GRCh38'
            },
            panels: [
                LocusZoom.Layouts.get('association', {
                    title: { text: `Chr${chr}:${this.formatPosition(start)}-${this.formatPosition(end)}` },
                    margin: { top: 35, bottom: 40, left: 70, right: 50 }
                })
            ]
        });
        
        // Create plot
        try {
            this.plot = LocusZoom.populate(`#${this.containerId}`, dataSources, layout);
        } catch (e) {
            console.error('LocusZoom error:', e);
            this.renderFallback(chr, start, end, snps);
        }
    }

    /**
     * Fallback D3 implementation when LocusZoom is not available
     */
    renderFallback(chr, start, end, snps) {
        // Clear container
        this.container.innerHTML = '';
        
        const margin = { top: 40, right: 30, bottom: 60, left: 70 };
        const width = this.container.clientWidth || this.options.width;
        const height = this.options.height;
        const innerWidth = width - margin.left - margin.right;
        const innerHeight = height - margin.top - margin.bottom;
        
        // Create SVG
        const svg = d3.select(this.container)
            .append('svg')
            .attr('width', width)
            .attr('height', height);
        
        const g = svg.append('g')
            .attr('transform', `translate(${margin.left}, ${margin.top})`);
        
        // Title
        svg.append('text')
            .attr('x', width / 2)
            .attr('y', 25)
            .attr('text-anchor', 'middle')
            .attr('fill', '#f0f4f8')
            .attr('font-size', '14px')
            .attr('font-weight', '600')
            .text(`Chromosome ${chr}: ${this.formatPosition(start)} - ${this.formatPosition(end)}`);
        
        // Scales
        const xScale = d3.scaleLinear()
            .domain([start, end])
            .range([0, innerWidth]);
        
        const yMax = Math.max(d3.max(snps, d => d.logp), 8);
        const yScale = d3.scaleLinear()
            .domain([0, yMax * 1.1])
            .range([innerHeight, 0]);
        
        // Find lead SNP for LD coloring
        const leadSNP = snps.reduce((a, b) => a.logp > b.logp ? a : b);
        
        // Color scale for pseudo-LD (based on distance)
        const colorScale = d3.scaleSequential(d3.interpolateRdYlBu)
            .domain([1, 0]);
        
        // Axes
        const xAxis = d3.axisBottom(xScale)
            .ticks(6)
            .tickFormat(d => this.formatPosition(d));
        
        g.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0, ${innerHeight})`)
            .call(xAxis)
            .selectAll('text')
            .attr('fill', '#94a3b8')
            .attr('font-size', '10px');
        
        g.append('text')
            .attr('x', innerWidth / 2)
            .attr('y', innerHeight + 45)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text(`Position on Chromosome ${chr} (Mb)`);
        
        const yAxis = d3.axisLeft(yScale).ticks(6);
        
        g.append('g')
            .attr('class', 'y-axis')
            .call(yAxis)
            .selectAll('text')
            .attr('fill', '#94a3b8')
            .attr('font-size', '10px');
        
        g.append('text')
            .attr('transform', 'rotate(-90)')
            .attr('x', -innerHeight / 2)
            .attr('y', -50)
            .attr('text-anchor', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '12px')
            .text('-log₁₀(p-value)');
        
        // Significance threshold line
        const sigY = yScale(-Math.log10(5e-8));
        if (sigY > 0) {
            g.append('line')
                .attr('x1', 0)
                .attr('x2', innerWidth)
                .attr('y1', sigY)
                .attr('y2', sigY)
                .attr('stroke', '#f59e0b')
                .attr('stroke-dasharray', '4,4')
                .attr('stroke-width', 1);
        }
        
        // Points
        g.selectAll('.locus-point')
            .data(snps)
            .enter()
            .append('circle')
            .attr('class', 'locus-point')
            .attr('cx', d => xScale(d.pos))
            .attr('cy', d => yScale(d.logp))
            .attr('r', d => d.id === leadSNP.id ? 8 : 5)
            .attr('fill', d => {
                if (d.id === leadSNP.id) return '#8b5cf6';
                // Pseudo-LD based on distance
                const distance = Math.abs(d.pos - leadSNP.pos);
                const maxDist = (end - start) / 2;
                const r2 = Math.max(0, 1 - distance / maxDist);
                return colorScale(r2);
            })
            .attr('stroke', d => d.id === leadSNP.id ? '#fff' : 'none')
            .attr('stroke-width', 2)
            .attr('opacity', 0.8)
            .style('cursor', 'pointer')
            .on('mouseenter', (event, d) => this.showTooltip(event, d))
            .on('mouseleave', () => this.hideTooltip());
        
        // Lead SNP label
        g.append('text')
            .attr('x', xScale(leadSNP.pos))
            .attr('y', yScale(leadSNP.logp) - 15)
            .attr('text-anchor', 'middle')
            .attr('fill', '#f0f4f8')
            .attr('font-size', '11px')
            .attr('font-weight', '600')
            .text(leadSNP.id);
        
        // Color legend
        this.renderLegend(svg, width, height);
        
        // Gene track (simplified)
        this.renderGeneTrack(g, innerWidth, innerHeight, xScale, snps);
    }

    /**
     * Render LD color legend
     */
    renderLegend(svg, width, height) {
        const legendWidth = 150;
        const legendHeight = 12;
        const legendX = width - 180;
        const legendY = 45;
        
        const legendGroup = svg.append('g')
            .attr('transform', `translate(${legendX}, ${legendY})`);
        
        // Gradient definition
        const defs = svg.append('defs');
        const gradient = defs.append('linearGradient')
            .attr('id', 'ld-gradient');
        
        gradient.append('stop')
            .attr('offset', '0%')
            .attr('stop-color', d3.interpolateRdYlBu(0));
        
        gradient.append('stop')
            .attr('offset', '50%')
            .attr('stop-color', d3.interpolateRdYlBu(0.5));
        
        gradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', d3.interpolateRdYlBu(1));
        
        // Legend rectangle
        legendGroup.append('rect')
            .attr('width', legendWidth)
            .attr('height', legendHeight)
            .attr('fill', 'url(#ld-gradient)')
            .attr('rx', 2);
        
        // Labels
        legendGroup.append('text')
            .attr('x', 0)
            .attr('y', -5)
            .attr('fill', '#94a3b8')
            .attr('font-size', '10px')
            .text('LD (r²)');
        
        legendGroup.append('text')
            .attr('x', 0)
            .attr('y', legendHeight + 12)
            .attr('fill', '#94a3b8')
            .attr('font-size', '9px')
            .text('0');
        
        legendGroup.append('text')
            .attr('x', legendWidth)
            .attr('y', legendHeight + 12)
            .attr('text-anchor', 'end')
            .attr('fill', '#94a3b8')
            .attr('font-size', '9px')
            .text('1');
    }

    /**
     * Render simplified gene track
     */
    renderGeneTrack(g, width, height, xScale, snps) {
        const trackHeight = 30;
        const trackY = height + 20;
        
        // Get unique genes in region
        const genes = [...new Set(snps.map(s => s.gene).filter(Boolean))];
        
        if (genes.length === 0) return;
        
        // Simple gene representation
        const geneGroup = g.append('g')
            .attr('class', 'gene-track')
            .attr('transform', `translate(0, ${trackY})`);
        
        // Track background
        geneGroup.append('rect')
            .attr('x', 0)
            .attr('y', 0)
            .attr('width', width)
            .attr('height', trackHeight)
            .attr('fill', 'rgba(255, 255, 255, 0.02)');
        
        // Track label
        geneGroup.append('text')
            .attr('x', -10)
            .attr('y', trackHeight / 2 + 4)
            .attr('text-anchor', 'end')
            .attr('fill', '#64748b')
            .attr('font-size', '10px')
            .text('Genes');
        
        // Distribute genes along track
        const geneWidth = Math.min(width / genes.length - 10, 80);
        
        genes.slice(0, 10).forEach((gene, i) => {
            const snpsWithGene = snps.filter(s => s.gene === gene);
            const avgPos = snpsWithGene.reduce((a, s) => a + s.pos, 0) / snpsWithGene.length;
            const x = xScale(avgPos);
            
            // Gene body
            geneGroup.append('rect')
                .attr('x', x - geneWidth / 2)
                .attr('y', 8)
                .attr('width', geneWidth)
                .attr('height', 14)
                .attr('fill', '#3b82f6')
                .attr('rx', 2)
                .attr('opacity', 0.6);
            
            // Gene name
            geneGroup.append('text')
                .attr('x', x)
                .attr('y', 18)
                .attr('text-anchor', 'middle')
                .attr('fill', '#f0f4f8')
                .attr('font-size', '9px')
                .attr('font-style', 'italic')
                .text(gene);
        });
    }

    /**
     * Format genomic position
     */
    formatPosition(pos) {
        if (pos >= 1e6) {
            return (pos / 1e6).toFixed(2) + ' Mb';
        } else if (pos >= 1e3) {
            return (pos / 1e3).toFixed(1) + ' kb';
        }
        return pos.toString();
    }

    /**
     * Show tooltip
     */
    showTooltip(event, snp) {
        const tooltip = document.getElementById('tooltip');
        if (!tooltip) return;
        
        tooltip.innerHTML = `
            <div class="tooltip-title">${snp.id}</div>
            <div class="tooltip-content">
                <div>Position: ${snp.chr}:${snp.pos.toLocaleString()}</div>
                <div>P-value: ${snp.pvalue.toExponential(2)}</div>
                <div>-log₁₀(P): ${snp.logp.toFixed(2)}</div>
                <div>Gene: ${snp.gene || 'N/A'}</div>
                <div>MAF: ${snp.maf?.toFixed(3) || 'N/A'}</div>
            </div>
        `;
        
        tooltip.style.left = event.pageX + 15 + 'px';
        tooltip.style.top = event.pageY - 10 + 'px';
        tooltip.classList.add('visible');
    }

    /**
     * Hide tooltip
     */
    hideTooltip() {
        const tooltip = document.getElementById('tooltip');
        if (tooltip) {
            tooltip.classList.remove('visible');
        }
    }

    /**
     * Show empty state
     */
    showEmpty(message = 'Select a region to visualize') {
        this.container.innerHTML = `
            <div class="locus-placeholder">
                <div class="placeholder-icon">🔍</div>
                <p>${message}</p>
            </div>
        `;
    }

    /**
     * Get current region info
     */
    getCurrentRegion() {
        return this.currentRegion;
    }

    /**
     * Cleanup
     */
    destroy() {
        if (this.plot) {
            this.plot.destroy();
        }
        this.container.innerHTML = '';
    }
}

// Export
window.LocusZoomIntegration = LocusZoomIntegration;
