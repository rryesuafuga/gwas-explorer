/**
 * GWAS Explorer - Main Application
 * Orchestrates all visualization components and user interactions
 */

class GWASExplorerApp {
    constructor() {
        // Initialize components
        this.simulator = null;
        this.manhattanPlot = null;
        this.qqPlot = null;
        this.locusZoom = null;
        this.prsCalculator = null;

        // Application state
        this.gwasData = null;
        this.gwasMetadata = null;
        this.significantLoci = null;
        this.selectedSNP = null;
        this.significanceThreshold = 1e-5;
        this.colorScheme = 'alternating';

        // DOM Elements
        this.elements = {};

        // Initialize on DOM ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', () => this.init());
        } else {
            this.init();
        }
    }
    
    /**
     * Initialize the application
     */
    async init() {
        console.log('Initializing GWAS Explorer...');

        // Initialize simulator
        this.simulator = new GWASSimulator({
            numSNPs: 50000,
            significantLoci: 12,
            genomicInflation: 1.02
        });

        // Cache DOM elements
        this.cacheElements();

        // Set up event listeners
        this.setupEventListeners();

        // Initialize visualizations
        this.initializeVisualizations();

        // Generate initial demo data
        await this.loadDemoData();

        // Initialize PRS calculator with TensorFlow.js
        await this.initializePRS();

        // Populate chromosome selector
        this.populateChromosomeSelector();

        console.log('GWAS Explorer ready!');
    }
    
    /**
     * Cache DOM element references
     */
    cacheElements() {
        this.elements = {
            // Navigation
            navLinks: document.querySelectorAll('.nav-link'),
            
            // Buttons
            btnSimulate: document.getElementById('btn-simulate'),
            btnUpload: document.getElementById('btn-upload'),
            fileInput: document.getElementById('file-input'),
            btnResetZoom: document.getElementById('btn-reset-zoom'),
            btnDownloadManhattan: document.getElementById('btn-download-manhattan'),
            
            // Stats display
            statSnps: document.getElementById('stat-snps'),
            statSig: document.getElementById('stat-sig'),
            statLambda: document.getElementById('stat-lambda'),
            
            // Plot containers
            manhattanPlot: document.getElementById('manhattan-plot'),
            manhattanLoading: document.getElementById('manhattan-loading'),
            qqPlot: document.getElementById('qq-plot'),
            chromosomePreview: document.getElementById('chromosome-preview'),
            
            // Toolbar controls
            thresholdBtns: document.querySelectorAll('.threshold-btn'),
            colorSchemeSelect: document.getElementById('color-scheme'),
            
            // SNP info panel
            snpInfoPanel: document.getElementById('snp-info-panel'),
            closeSnpPanel: document.getElementById('close-snp-panel'),
            
            // LocusZoom controls
            locusChr: document.getElementById('locus-chr'),
            locusStart: document.getElementById('locus-start'),
            locusEnd: document.getElementById('locus-end'),
            btnViewRegion: document.getElementById('btn-view-region'),
            quickLoci: document.getElementById('quick-loci'),
            locusZoomPlot: document.getElementById('locus-zoom-plot'),
            
            // PRS controls
            modelStatus: document.getElementById('model-status'),
            prsTabs: document.querySelectorAll('.prs-tab'),
            genotypeGrid: document.getElementById('genotype-grid'),
            simPopulation: document.getElementById('sim-population'),
            btnSimulateGenotypes: document.getElementById('btn-simulate-genotypes'),
            vcfUploadZone: document.getElementById('vcf-upload-zone'),
            vcfInput: document.getElementById('vcf-input'),
            btnCalculatePRS: document.getElementById('btn-calculate-prs'),
            
            // QQ stats
            qqLambda: document.getElementById('qq-lambda'),
            qqLambdaStatus: document.getElementById('qq-lambda-status'),
            qqLambda1000: document.getElementById('qq-lambda-1000'),
            qqLambda1000Status: document.getElementById('qq-lambda-1000-status'),
            qqMeanChi2: document.getElementById('qq-mean-chi2'),
            qqDeviation: document.getElementById('qq-deviation'),
            qqDeviationStatus: document.getElementById('qq-deviation-status'),
            qqInterpretation: document.getElementById('qq-interpretation'),
            
            // Tooltip
            tooltip: document.getElementById('tooltip')
        };
    }
    
    /**
     * Set up event listeners
     */
    setupEventListeners() {
        // Navigation smooth scroll
        this.elements.navLinks.forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').slice(1);
                document.getElementById(targetId)?.scrollIntoView({ behavior: 'smooth' });
                
                // Update active state
                this.elements.navLinks.forEach(l => l.classList.remove('active'));
                link.classList.add('active');
            });
        });
        
        // Simulate GWAS button
        this.elements.btnSimulate?.addEventListener('click', () => this.simulateNewGWAS());
        
        // Upload button
        this.elements.btnUpload?.addEventListener('click', () => {
            this.elements.fileInput?.click();
        });
        
        // File input
        this.elements.fileInput?.addEventListener('change', (e) => {
            if (e.target.files.length > 0) {
                this.handleFileUpload(e.target.files[0]);
            }
        });
        
        // Reset zoom
        this.elements.btnResetZoom?.addEventListener('click', () => {
            this.manhattanPlot?.resetZoom();
        });
        
        // Download Manhattan plot
        this.elements.btnDownloadManhattan?.addEventListener('click', () => {
            this.manhattanPlot?.exportSVG();
        });
        
        // Significance threshold buttons
        this.elements.thresholdBtns.forEach(btn => {
            btn.addEventListener('click', () => {
                this.elements.thresholdBtns.forEach(b => b.classList.remove('active'));
                btn.classList.add('active');

                this.significanceThreshold = parseFloat(btn.dataset.threshold);
                this.manhattanPlot?.setThreshold(this.significanceThreshold);
            });
        });

        // Color scheme
        this.elements.colorSchemeSelect?.addEventListener('change', (e) => {
            this.colorScheme = e.target.value;
            this.manhattanPlot?.setColorScheme(this.colorScheme);
        });
        
        // Close SNP panel
        this.elements.closeSnpPanel?.addEventListener('click', () => {
            this.elements.snpInfoPanel?.classList.remove('visible');
        });
        
        // LocusZoom region viewer
        this.elements.btnViewRegion?.addEventListener('click', () => {
            const chr = this.elements.locusChr.value;
            const start = parseInt(this.elements.locusStart.value);
            const end = parseInt(this.elements.locusEnd.value);
            
            if (chr && start && end) {
                this.viewLocusRegion(parseInt(chr), start, end);
            }
        });
        
        // PRS tabs
        this.elements.prsTabs.forEach(tab => {
            tab.addEventListener('click', () => {
                this.elements.prsTabs.forEach(t => t.classList.remove('active'));
                tab.classList.add('active');
                
                document.querySelectorAll('.prs-tab-content').forEach(content => {
                    content.classList.remove('active');
                });
                document.getElementById(`tab-${tab.dataset.tab}`)?.classList.add('active');
            });
        });
        
        // Simulate genotypes button
        this.elements.btnSimulateGenotypes?.addEventListener('click', () => {
            this.simulateGenotypes();
        });
        
        // VCF upload zone
        const vcfZone = this.elements.vcfUploadZone;
        if (vcfZone) {
            vcfZone.addEventListener('click', () => this.elements.vcfInput?.click());
            vcfZone.addEventListener('dragover', (e) => {
                e.preventDefault();
                vcfZone.classList.add('dragover');
            });
            vcfZone.addEventListener('dragleave', () => {
                vcfZone.classList.remove('dragover');
            });
            vcfZone.addEventListener('drop', (e) => {
                e.preventDefault();
                vcfZone.classList.remove('dragover');
                if (e.dataTransfer.files.length > 0) {
                    this.handleVCFUpload(e.dataTransfer.files[0]);
                }
            });
        }
        
        this.elements.vcfInput?.addEventListener('change', (e) => {
            if (e.target.files.length > 0) {
                this.handleVCFUpload(e.target.files[0]);
            }
        });
        
        // Calculate PRS button
        this.elements.btnCalculatePRS?.addEventListener('click', () => {
            this.calculatePRS();
        });
        
        // Window resize
        window.addEventListener('resize', this.debounce(() => {
            this.manhattanPlot?.resize();
            this.qqPlot?.resize();
        }, 250));
        
        // Intersection observer for section highlighting
        this.setupScrollObserver();
    }
    
    /**
     * Initialize visualization components
     */
    initializeVisualizations() {
        // Manhattan plot
        if (this.elements.manhattanPlot) {
            this.manhattanPlot = new ManhattanPlot('manhattan-plot', {
                significanceThreshold: this.significanceThreshold,
                colorScheme: this.colorScheme
            });
            // Set up callbacks
            this.manhattanPlot.onSNPClick = (snp) => this.handleSNPClick(snp);
            this.manhattanPlot.onSNPHover = (snp, event) => this.handleSNPHover(snp, event);
        }

        // QQ plot
        if (this.elements.qqPlot) {
            this.qqPlot = new QQPlot('qq-plot', {});
        }

        // LocusZoom integration
        if (this.elements.locusZoomPlot) {
            this.locusZoom = new LocusZoomIntegration('locus-zoom-plot');
        }

        // Chromosome preview (mini visualization in hero)
        this.renderChromosomePreview();
    }
    
    /**
     * Load demo data
     */
    async loadDemoData() {
        this.showLoading(true);

        try {
            // Simulate GWAS data using the generator
            const result = this.simulator.generate();
            this.gwasData = result.snps;
            this.gwasMetadata = result.metadata;
            this.significantLoci = result.significantLoci;

            // Update visualizations
            await this.updateVisualizations();

            // Update stats
            this.updateStats();

            // Populate quick loci
            this.populateQuickLoci();

        } catch (error) {
            console.error('Error loading demo data:', error);
        } finally {
            this.showLoading(false);
        }
    }
    
    /**
     * Simulate new GWAS data
     */
    async simulateNewGWAS() {
        this.showLoading(true);

        // Add animation to button
        const btn = this.elements.btnSimulate;
        if (btn) {
            btn.disabled = true;
            const icon = btn.querySelector('.btn-icon');
            if (icon) icon.style.animation = 'spin 1s linear infinite';
        }

        try {
            // Generate new random parameters
            const numSNPs = 30000 + Math.floor(Math.random() * 40000);
            const significantLoci = 5 + Math.floor(Math.random() * 15);
            const genomicInflation = 1.0 + Math.random() * 0.08;

            // Create new simulator with random params
            this.simulator = new GWASSimulator({
                numSNPs,
                significantLoci,
                genomicInflation
            });

            const result = this.simulator.generate();
            this.gwasData = result.snps;
            this.gwasMetadata = result.metadata;
            this.significantLoci = result.significantLoci;

            await this.updateVisualizations();
            this.updateStats();
            this.populateQuickLoci();

        } catch (error) {
            console.error('Error simulating GWAS:', error);
        } finally {
            this.showLoading(false);
            if (btn) {
                btn.disabled = false;
                const icon = btn.querySelector('.btn-icon');
                if (icon) icon.style.animation = '';
            }
        }
    }
    
    /**
     * Update all visualizations with current data
     */
    async updateVisualizations() {
        if (!this.gwasData) return;

        // Prepare data in format expected by visualizations
        const plotData = {
            snps: this.gwasData,
            metadata: this.gwasMetadata
        };

        // Update Manhattan plot
        if (this.manhattanPlot) {
            await this.manhattanPlot.render(plotData);
        }

        // Update QQ plot
        if (this.qqPlot) {
            const stats = await this.qqPlot.render(plotData);
            if (stats) {
                this.updateQQStats(stats);
            }
        }

        // Update LocusZoom data
        if (this.locusZoom) {
            this.locusZoom.setData(plotData);
        }

        // Update chromosome preview
        this.renderChromosomePreview();
    }
    
    /**
     * Update statistics display
     */
    updateStats() {
        if (!this.gwasData || !this.gwasMetadata) return;

        const significant = this.gwasData.filter(d => d.pvalue < this.significanceThreshold);

        // Animate number updates
        this.animateNumber(this.elements.statSnps, this.gwasData.length);
        this.animateNumber(this.elements.statSig, significant.length);
        if (this.elements.statLambda) {
            this.elements.statLambda.textContent = this.gwasMetadata.lambda?.toFixed(3) || '-';
        }
    }
    
    /**
     * Update QQ plot statistics panel
     */
    updateQQStats(stats) {
        if (!stats) return;

        const { lambda, lambda1000, meanChi2, tailDeviation } = stats;

        // Update display
        if (this.elements.qqLambda) {
            this.elements.qqLambda.textContent = lambda?.toFixed(3) || '-';
        }
        if (this.elements.qqLambda1000) {
            this.elements.qqLambda1000.textContent = lambda1000?.toFixed(3) || '-';
        }
        if (this.elements.qqMeanChi2) {
            this.elements.qqMeanChi2.textContent = meanChi2?.toFixed(3) || '-';
        }
        if (this.elements.qqDeviation) {
            this.elements.qqDeviation.textContent = tailDeviation?.toFixed(2) || '-';
        }

        // Update status indicators
        if (this.elements.qqLambdaStatus) {
            if (lambda < 1.05) {
                this.elements.qqLambdaStatus.textContent = 'Excellent (no inflation)';
                this.elements.qqLambdaStatus.className = 'stat-card-status good';
            } else if (lambda < 1.1) {
                this.elements.qqLambdaStatus.textContent = 'Acceptable';
                this.elements.qqLambdaStatus.className = 'stat-card-status warning';
            } else {
                this.elements.qqLambdaStatus.textContent = 'Inflated - check stratification';
                this.elements.qqLambdaStatus.className = 'stat-card-status bad';
            }
        }

        // Interpretation
        if (this.elements.qqInterpretation) {
            let interpretation = '';
            if (lambda < 1.05 && tailDeviation > 0.5) {
                interpretation = 'The QQ plot shows minimal genomic inflation with clear deviation in the tail, suggesting true genetic associations. The study appears well-controlled with genuine signals.';
            } else if (lambda > 1.1) {
                interpretation = 'Elevated genomic inflation factor suggests possible population stratification or cryptic relatedness. Consider using mixed models or additional principal components.';
            } else {
                interpretation = 'The distribution shows expected behavior under the null hypothesis with some deviation indicating potential associations. Further fine-mapping is recommended for significant loci.';
            }
            this.elements.qqInterpretation.textContent = interpretation;
        }
    }
    
    /**
     * Initialize PRS calculator
     */
    async initializePRS() {
        try {
            this.prsCalculator = new PRSCalculator({
                numVariants: 20
            });

            this.prsCalculator.onModelLoaded = () => {
                this.updateModelStatus('ready', 'Model loaded successfully');
                if (this.elements.btnCalculatePRS) {
                    this.elements.btnCalculatePRS.disabled = false;
                }
            };

            this.prsCalculator.onError = (error) => {
                this.updateModelStatus('error', `Model error: ${error.message}`);
            };

            // Initialize the calculator
            await this.prsCalculator.initialize();

            // Populate genotype input grid
            this.populateGenotypeGrid();

        } catch (error) {
            console.error('Error initializing PRS:', error);
            this.updateModelStatus('error', 'Failed to initialize model');
        }
    }
    
    /**
     * Update model status display
     */
    updateModelStatus(status, message) {
        const indicator = this.elements.modelStatus.querySelector('.status-indicator');
        const text = this.elements.modelStatus.querySelector('span');
        
        indicator.className = `status-indicator ${status}`;
        text.textContent = message;
    }
    
    /**
     * Populate genotype input grid
     */
    populateGenotypeGrid() {
        if (!this.prsCalculator || !this.elements.genotypeGrid) return;

        const variants = this.prsCalculator.getRiskVariants();
        if (!variants) return;

        this.elements.genotypeGrid.innerHTML = variants.map((variant, idx) => `
            <div class="genotype-item">
                <span class="genotype-snp">${variant.id}</span>
                <input type="number"
                       class="genotype-input"
                       data-idx="${idx}"
                       data-snp="${variant.id}"
                       min="0" max="2" step="1"
                       value="0"
                       placeholder="0-2">
            </div>
        `).join('');
    }
    
    /**
     * Simulate genotypes based on population
     */
    simulateGenotypes() {
        if (!this.prsCalculator) return;

        const population = this.elements.simPopulation?.value || 'EUR';
        const genotypes = this.prsCalculator.simulateGenotypes(population);

        // Update input fields
        const inputs = this.elements.genotypeGrid?.querySelectorAll('.genotype-input');
        if (inputs) {
            inputs.forEach((input, idx) => {
                if (genotypes[idx] !== undefined) {
                    input.value = genotypes[idx];
                    // Animate the change
                    input.style.background = 'rgba(59, 130, 246, 0.3)';
                    setTimeout(() => {
                        input.style.background = '';
                    }, 500);
                }
            });
        }

        if (this.elements.btnCalculatePRS) {
            this.elements.btnCalculatePRS.disabled = false;
        }
    }
    
    /**
     * Calculate PRS
     */
    async calculatePRS() {
        if (!this.prsCalculator) return;

        // Collect genotypes from inputs as array
        const genotypes = [];
        const inputs = this.elements.genotypeGrid?.querySelectorAll('.genotype-input');
        if (inputs) {
            inputs.forEach(input => {
                genotypes.push(parseInt(input.value) || 0);
            });
        }

        try {
            const result = await this.prsCalculator.calculatePRS(genotypes);
            this.displayPRSResults(result);
        } catch (error) {
            console.error('Error calculating PRS:', error);
        }
    }
    
    /**
     * Display PRS results
     */
    displayPRSResults(result) {
        // Show results panel
        document.querySelector('.result-placeholder').style.display = 'none';
        document.getElementById('prs-results').classList.remove('hidden');
        
        // Update score display
        document.getElementById('prs-score-value').textContent = result.score.toFixed(3);
        
        // Animate ring progress
        const circumference = 2 * Math.PI * 90;
        const progress = Math.min(result.percentile / 100, 1);
        const offset = circumference * (1 - progress);
        
        const ring = document.getElementById('ring-progress');
        ring.style.strokeDashoffset = offset;
        
        // Update percentile bar
        document.getElementById('percentile-fill').style.width = `${result.percentile}%`;
        document.getElementById('percentile-marker').style.left = `${result.percentile}%`;
        document.getElementById('percentile-value').textContent = `${result.percentile.toFixed(1)}th`;
        
        // Render breakdown chart
        this.renderPRSBreakdown(result.contributions);
    }
    
    /**
     * Render PRS contribution breakdown chart
     */
    renderPRSBreakdown(contributions) {
        const container = document.getElementById('breakdown-chart');
        if (!container || !contributions) return;
        
        const width = container.clientWidth;
        const height = 150;
        const margin = { top: 10, right: 20, bottom: 30, left: 60 };
        
        // Clear previous
        container.innerHTML = '';
        
        const svg = d3.select(container)
            .append('svg')
            .attr('width', width)
            .attr('height', height);
        
        // Sort by absolute contribution
        const data = Object.entries(contributions)
            .map(([snp, value]) => ({ snp: snp.slice(0, 10), value }))
            .sort((a, b) => Math.abs(b.value) - Math.abs(a.value))
            .slice(0, 10);
        
        const x = d3.scaleLinear()
            .domain([d3.min(data, d => d.value) * 1.1, d3.max(data, d => d.value) * 1.1])
            .range([margin.left, width - margin.right]);
        
        const y = d3.scaleBand()
            .domain(data.map(d => d.snp))
            .range([margin.top, height - margin.bottom])
            .padding(0.2);
        
        // Zero line
        svg.append('line')
            .attr('x1', x(0))
            .attr('x2', x(0))
            .attr('y1', margin.top)
            .attr('y2', height - margin.bottom)
            .attr('stroke', 'rgba(255,255,255,0.2)');
        
        // Bars
        svg.selectAll('.bar')
            .data(data)
            .join('rect')
            .attr('class', 'bar')
            .attr('x', d => d.value >= 0 ? x(0) : x(d.value))
            .attr('y', d => y(d.snp))
            .attr('width', d => Math.abs(x(d.value) - x(0)))
            .attr('height', y.bandwidth())
            .attr('fill', d => d.value >= 0 ? '#3b82f6' : '#ef4444')
            .attr('rx', 2);
        
        // Labels
        svg.selectAll('.label')
            .data(data)
            .join('text')
            .attr('class', 'label')
            .attr('x', margin.left - 5)
            .attr('y', d => y(d.snp) + y.bandwidth() / 2)
            .attr('text-anchor', 'end')
            .attr('dominant-baseline', 'middle')
            .attr('fill', '#94a3b8')
            .attr('font-size', '10px')
            .attr('font-family', 'JetBrains Mono')
            .text(d => d.snp);
    }
    
    /**
     * Populate chromosome selector
     */
    populateChromosomeSelector() {
        if (!this.elements.locusChr) return;
        
        this.elements.locusChr.innerHTML = '<option value="">Chr</option>' +
            Array.from({ length: 22 }, (_, i) => i + 1)
                .map(chr => `<option value="${chr}">${chr}</option>`)
                .join('');
    }
    
    /**
     * Populate quick loci shortcuts
     */
    populateQuickLoci() {
        if (!this.gwasData || !this.elements.quickLoci) return;

        // Get significant SNPs - use isSignificant field and sort by pvalue
        const significant = this.gwasData
            .filter(d => d.isSignificant)
            .sort((a, b) => a.pvalue - b.pvalue)
            .slice(0, 8);

        if (significant.length === 0) {
            this.elements.quickLoci.innerHTML = '<span class="empty-state">No significant loci found</span>';
            return;
        }

        this.elements.quickLoci.innerHTML = significant.map(snp => `
            <button class="locus-chip"
                    data-chr="${snp.chr}"
                    data-pos="${snp.pos}"
                    title="${snp.id} (p=${snp.pvalue.toExponential(2)})">
                ${snp.chr}:${(snp.pos / 1e6).toFixed(1)}Mb
            </button>
        `).join('');

        // Add click handlers
        this.elements.quickLoci.querySelectorAll('.locus-chip').forEach(chip => {
            chip.addEventListener('click', () => {
                const chr = parseInt(chip.dataset.chr);
                const pos = parseInt(chip.dataset.pos);
                const windowSize = 500000; // 500kb window
                this.viewLocusRegion(chr, pos - windowSize, pos + windowSize);
            });
        });
    }
    
    /**
     * View a specific genomic region in LocusZoom
     */
    viewLocusRegion(chr, start, end) {
        if (!this.locusZoom || !this.gwasData) return;
        
        // Filter data for this region
        const regionData = this.gwasData.filter(d => 
            d.chr === chr && d.pos >= start && d.pos <= end
        );
        
        // Update LocusZoom
        this.locusZoom.render(regionData, { chr, start, end });
        
        // Update input fields
        this.elements.locusChr.value = chr;
        this.elements.locusStart.value = start;
        this.elements.locusEnd.value = end;
        
        // Scroll to section
        document.getElementById('locus')?.scrollIntoView({ behavior: 'smooth' });
    }
    
    /**
     * Handle SNP click in Manhattan plot
     */
    handleSNPClick(snp) {
        this.selectedSNP = snp;

        // Update info panel
        const setElement = (id, value) => {
            const el = document.getElementById(id);
            if (el) el.textContent = value;
        };

        setElement('info-snp-id', snp.id || '-');
        setElement('info-chr', snp.chr);
        setElement('info-pos', snp.pos?.toLocaleString() || '-');
        setElement('info-pval', snp.pvalue?.toExponential(2) || '-');
        setElement('info-logp', snp.logp?.toFixed(2) || '-');
        setElement('info-maf', snp.maf?.toFixed(3) || '-');
        setElement('info-beta', snp.beta?.toFixed(4) || '-');
        setElement('info-gene', snp.gene || '-');

        // Update external link
        const externalLink = document.getElementById('btn-external-link');
        if (externalLink) {
            const ensemblLink = `https://www.ensembl.org/Homo_sapiens/Location/View?r=${snp.chr}:${snp.pos - 50000}-${snp.pos + 50000}`;
            externalLink.href = ensemblLink;
        }

        // Show panel
        if (this.elements.snpInfoPanel) {
            this.elements.snpInfoPanel.classList.add('visible');
        }

        // Set up LocusZoom button
        const zoomBtn = document.getElementById('btn-zoom-locus');
        if (zoomBtn) {
            zoomBtn.onclick = () => {
                this.viewLocusRegion(snp.chr, snp.pos - 500000, snp.pos + 500000);
            };
        }
    }
    
    /**
     * Handle SNP hover
     */
    handleSNPHover(snp, event) {
        if (!snp) {
            this.hideTooltip();
            return;
        }
        const content = `
            <div class="tooltip-title">${snp.id || 'SNP'}</div>
            <div class="tooltip-content">
                Chr${snp.chr}:${snp.pos?.toLocaleString()}<br>
                P = ${snp.pvalue?.toExponential(2)}<br>
                -log10(P) = ${snp.logp?.toFixed(2)}
                ${snp.gene ? `<br>Gene: ${snp.gene}` : ''}
            </div>
        `;
        this.showTooltip(content, event);
    }
    
    /**
     * Show tooltip
     */
    showTooltip(content, event) {
        const tooltip = this.elements.tooltip;
        tooltip.innerHTML = content;
        tooltip.classList.add('visible');
        
        // Position tooltip
        const x = event.pageX + 15;
        const y = event.pageY - 10;
        
        tooltip.style.left = `${x}px`;
        tooltip.style.top = `${y}px`;
        
        // Adjust if off screen
        const rect = tooltip.getBoundingClientRect();
        if (rect.right > window.innerWidth) {
            tooltip.style.left = `${event.pageX - rect.width - 15}px`;
        }
    }
    
    /**
     * Hide tooltip
     */
    hideTooltip() {
        this.elements.tooltip.classList.remove('visible');
    }
    
    /**
     * Render chromosome preview in hero section
     */
    renderChromosomePreview() {
        if (!this.elements.chromosomePreview) return;

        const container = this.elements.chromosomePreview;
        const width = container.clientWidth;
        const height = container.clientHeight;

        // Clear previous
        container.innerHTML = '';

        if (!this.gwasData || this.gwasData.length === 0) {
            // Placeholder animation
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: var(--text-muted);">
                    <span>Loading preview...</span>
                </div>
            `;
            return;
        }

        const svg = d3.select(container)
            .append('svg')
            .attr('width', width)
            .attr('height', height);

        // Create a mini Manhattan plot
        const margin = { top: 20, right: 20, bottom: 20, left: 20 };
        const plotWidth = width - margin.left - margin.right;
        const plotHeight = height - margin.top - margin.bottom;

        const g = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Sample data for preview
        const sampleData = this.gwasData.filter((_, i) => i % 20 === 0);

        const x = d3.scaleLinear()
            .domain(d3.extent(sampleData, d => d.cumPos))
            .range([0, plotWidth]);

        const y = d3.scaleLinear()
            .domain([0, d3.max(sampleData, d => d.logp) * 1.1])
            .range([plotHeight, 0]);

        // Points
        g.selectAll('circle')
            .data(sampleData)
            .join('circle')
            .attr('cx', d => x(d.cumPos))
            .attr('cy', d => y(d.logp))
            .attr('r', d => d.isSignificant ? 3 : 1.5)
            .attr('fill', d => d.isSignificant ? '#f59e0b' : (d.chr % 2 === 0 ? '#3b82f6' : '#06b6d4'))
            .attr('opacity', d => d.isSignificant ? 1 : 0.5);

        // Threshold line
        const thresholdY = y(-Math.log10(this.significanceThreshold));
        if (thresholdY > 0 && thresholdY < plotHeight) {
            g.append('line')
                .attr('x1', 0)
                .attr('x2', plotWidth)
                .attr('y1', thresholdY)
                .attr('y2', thresholdY)
                .attr('stroke', '#f59e0b')
                .attr('stroke-dasharray', '4,2')
                .attr('opacity', 0.7);
        }
    }
    
    /**
     * Handle file upload
     */
    async handleFileUpload(file) {
        this.showLoading(true);
        
        try {
            const text = await file.text();
            let data;
            
            if (file.name.endsWith('.json')) {
                data = JSON.parse(text);
            } else {
                // Parse CSV/TSV
                data = this.parseGWASText(text);
            }
            
            if (data && data.length > 0) {
                this.gwasData = data;
                await this.updateVisualizations();
                this.updateStats();
                this.populateQuickLoci();
            }
        } catch (error) {
            console.error('Error loading file:', error);
            alert('Error loading file. Please check the format.');
        } finally {
            this.showLoading(false);
        }
    }
    
    /**
     * Parse GWAS text file (CSV/TSV)
     */
    parseGWASText(text) {
        const lines = text.trim().split('\n');
        const header = lines[0].toLowerCase().split(/[\t,]/);
        
        // Find column indices
        const chrIdx = header.findIndex(h => h.includes('chr'));
        const posIdx = header.findIndex(h => h.includes('pos') || h.includes('bp'));
        const pvalIdx = header.findIndex(h => h.includes('pval') || h.includes('p_value') || h === 'p');
        const snpIdx = header.findIndex(h => h.includes('snp') || h.includes('rsid') || h.includes('variant'));
        
        if (chrIdx === -1 || posIdx === -1 || pvalIdx === -1) {
            throw new Error('Could not find required columns (chr, pos, pval)');
        }
        
        const data = [];
        let cumPos = 0;
        const chrOffsets = {};
        
        for (let i = 1; i < lines.length; i++) {
            const cols = lines[i].split(/[\t,]/);
            const chr = parseInt(cols[chrIdx].replace('chr', ''));
            const pos = parseInt(cols[posIdx]);
            const pval = parseFloat(cols[pvalIdx]);
            
            if (isNaN(chr) || isNaN(pos) || isNaN(pval)) continue;
            
            if (!chrOffsets[chr]) {
                chrOffsets[chr] = cumPos;
                if (chr > 1) {
                    cumPos += this.simulator.chromosomeLengths[chr - 1] || 200000000;
                }
            }
            
            data.push({
                snp: snpIdx !== -1 ? cols[snpIdx] : `chr${chr}:${pos}`,
                chr,
                pos,
                pval,
                logP: -Math.log10(pval),
                cumPos: chrOffsets[chr] + pos,
                isSignificant: pval < this.significanceThreshold
            });
        }
        
        return data;
    }
    
    /**
     * Handle VCF file upload
     */
    async handleVCFUpload(file) {
        // Basic VCF parsing for demo purposes
        console.log('VCF upload:', file.name);
        alert('VCF parsing is a demonstration feature. Please use manual entry or simulation for now.');
    }
    
    /**
     * Show/hide loading indicator
     */
    showLoading(show) {
        if (this.elements.manhattanLoading) {
            this.elements.manhattanLoading.classList.toggle('visible', show);
        }
    }
    
    /**
     * Animate number change
     */
    animateNumber(element, target) {
        if (!element) return;
        
        const start = parseInt(element.textContent) || 0;
        const duration = 1000;
        const startTime = performance.now();
        
        const animate = (currentTime) => {
            const elapsed = currentTime - startTime;
            const progress = Math.min(elapsed / duration, 1);
            
            // Easing function
            const easeOut = 1 - Math.pow(1 - progress, 3);
            
            const current = Math.floor(start + (target - start) * easeOut);
            element.textContent = current.toLocaleString();
            
            if (progress < 1) {
                requestAnimationFrame(animate);
            }
        };
        
        requestAnimationFrame(animate);
    }
    
    /**
     * Set up intersection observer for section highlighting
     */
    setupScrollObserver() {
        const sections = document.querySelectorAll('.section');
        
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const id = entry.target.id;
                    this.elements.navLinks.forEach(link => {
                        link.classList.toggle('active', link.getAttribute('href') === `#${id}`);
                    });
                }
            });
        }, { threshold: 0.3 });
        
        sections.forEach(section => observer.observe(section));
    }
    
    /**
     * Debounce utility
     */
    debounce(func, wait) {
        let timeout;
        return function executedFunction(...args) {
            const later = () => {
                clearTimeout(timeout);
                func(...args);
            };
            clearTimeout(timeout);
            timeout = setTimeout(later, wait);
        };
    }
}

// Initialize application
const app = new GWASExplorerApp();
