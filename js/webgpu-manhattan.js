/**
 * WebGPU Manhattan Plot Renderer
 * GPU-accelerated rendering for millions of SNPs
 * Falls back to Canvas 2D if WebGPU unavailable
 */

class WebGPUManhattanPlot {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        this.options = {
            pointSize: 3.0,
            significanceThreshold: 5e-8,
            suggestiveThreshold: 1e-5,
            ...options
        };

        this.device = null;
        this.context = null;
        this.pipeline = null;
        this.vertexBuffer = null;
        this.colorBuffer = null;
        this.uniformBuffer = null;
        this.numPoints = 0;
        this.isWebGPUAvailable = false;
        this.canvas = null;

        // Transform state
        this.scale = 1.0;
        this.offsetX = 0;
        this.offsetY = 0;

        this.init();
    }

    async init() {
        // Create canvas
        this.canvas = document.createElement('canvas');
        this.canvas.style.width = '100%';
        this.canvas.style.height = '100%';
        this.container.appendChild(this.canvas);

        this.resize();

        // Try WebGPU first
        if (navigator.gpu) {
            try {
                const adapter = await navigator.gpu.requestAdapter();
                if (adapter) {
                    this.device = await adapter.requestDevice();
                    this.context = this.canvas.getContext('webgpu');

                    const format = navigator.gpu.getPreferredCanvasFormat();
                    this.context.configure({
                        device: this.device,
                        format: format,
                        alphaMode: 'premultiplied'
                    });

                    await this.createPipeline(format);
                    this.isWebGPUAvailable = true;
                    console.log('WebGPU initialized successfully');
                }
            } catch (e) {
                console.warn('WebGPU initialization failed:', e);
            }
        }

        if (!this.isWebGPUAvailable) {
            console.log('Falling back to Canvas 2D');
            this.ctx = this.canvas.getContext('2d');
        }

        // Setup resize observer
        new ResizeObserver(() => this.resize()).observe(this.container);

        // Setup mouse interactions
        this.setupInteractions();
    }

    async createPipeline(format) {
        // Vertex shader - transforms points and calculates size
        const vertexShader = `
            struct Uniforms {
                resolution: vec2f,
                scale: f32,
                offsetX: f32,
                offsetY: f32,
                pointSize: f32,
                maxLogP: f32,
                padding: f32,
            }

            @group(0) @binding(0) var<uniform> uniforms: Uniforms;

            struct VertexInput {
                @location(0) position: vec2f,  // cumPos, logP
                @location(1) color: vec4f,
            }

            struct VertexOutput {
                @builtin(position) position: vec4f,
                @location(0) color: vec4f,
                @location(1) @interpolate(flat) pointSize: f32,
            }

            @vertex
            fn vertexMain(input: VertexInput, @builtin(vertex_index) vertexIndex: u32) -> VertexOutput {
                var output: VertexOutput;

                // Transform to normalized device coordinates
                let x = (input.position.x * uniforms.scale + uniforms.offsetX) / uniforms.resolution.x * 2.0 - 1.0;
                let y = (input.position.y / uniforms.maxLogP) * 2.0 - 1.0;

                output.position = vec4f(x, y, 0.0, 1.0);
                output.color = input.color;
                output.pointSize = uniforms.pointSize;

                return output;
            }
        `;

        // Fragment shader - renders circular points with antialiasing
        const fragmentShader = `
            @fragment
            fn fragmentMain(
                @location(0) color: vec4f,
                @location(1) @interpolate(flat) pointSize: f32,
                @builtin(position) fragCoord: vec4f
            ) -> @location(0) vec4f {
                return color;
            }
        `;

        // For point rendering, we'll use a compute shader approach
        // Create vertex shader module
        const shaderModule = this.device.createShaderModule({
            code: `
                struct Uniforms {
                    resolution: vec2f,
                    scale: f32,
                    offsetX: f32,
                    offsetY: f32,
                    pointSize: f32,
                    maxLogP: f32,
                    minX: f32,
                    maxX: f32,
                    padding: f32,
                }

                @group(0) @binding(0) var<uniform> uniforms: Uniforms;

                struct VertexOutput {
                    @builtin(position) position: vec4f,
                    @location(0) color: vec4f,
                }

                @vertex
                fn vertexMain(
                    @location(0) pos: vec2f,
                    @location(1) color: vec4f
                ) -> VertexOutput {
                    var output: VertexOutput;

                    // Map cumPos to screen X (0 to width)
                    let normalizedX = (pos.x - uniforms.minX) / (uniforms.maxX - uniforms.minX);
                    let screenX = normalizedX * uniforms.scale + uniforms.offsetX;
                    let clipX = (screenX / uniforms.resolution.x) * 2.0 - 1.0;

                    // Map logP to screen Y (bottom to top)
                    let normalizedY = pos.y / uniforms.maxLogP;
                    let clipY = normalizedY * 2.0 - 1.0;

                    output.position = vec4f(clipX, clipY, 0.0, 1.0);
                    output.color = color;

                    return output;
                }

                @fragment
                fn fragmentMain(@location(0) color: vec4f) -> @location(0) vec4f {
                    return color;
                }
            `
        });

        // Create pipeline
        this.pipeline = this.device.createRenderPipeline({
            layout: 'auto',
            vertex: {
                module: shaderModule,
                entryPoint: 'vertexMain',
                buffers: [
                    {
                        arrayStride: 8, // 2 floats for position
                        attributes: [
                            { shaderLocation: 0, offset: 0, format: 'float32x2' }
                        ]
                    },
                    {
                        arrayStride: 16, // 4 floats for color
                        attributes: [
                            { shaderLocation: 1, offset: 0, format: 'float32x4' }
                        ]
                    }
                ]
            },
            fragment: {
                module: shaderModule,
                entryPoint: 'fragmentMain',
                targets: [{ format: format }]
            },
            primitive: {
                topology: 'point-list'
            }
        });

        // Create uniform buffer
        this.uniformBuffer = this.device.createBuffer({
            size: 48, // 10 floats + padding
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });
    }

    resize() {
        const rect = this.container.getBoundingClientRect();
        const dpr = window.devicePixelRatio || 1;
        this.canvas.width = rect.width * dpr;
        this.canvas.height = rect.height * dpr;
        this.width = rect.width;
        this.height = rect.height;
    }

    /**
     * Render SNP data
     * @param {Object} data - { snps: [{cumPos, logp, chr, isSignificant}], metadata: {maxLogP} }
     */
    async render(data) {
        if (!data || !data.snps || data.snps.length === 0) return;

        this.data = data;
        this.numPoints = data.snps.length;

        const positions = new Float32Array(this.numPoints * 2);
        const colors = new Float32Array(this.numPoints * 4);

        // Calculate bounds
        let minX = Infinity, maxX = -Infinity;
        data.snps.forEach(snp => {
            if (snp.cumPos < minX) minX = snp.cumPos;
            if (snp.cumPos > maxX) maxX = snp.cumPos;
        });
        this.minX = minX;
        this.maxX = maxX;
        this.maxLogP = data.metadata?.maxLogP || 15;

        // Fill buffers
        for (let i = 0; i < this.numPoints; i++) {
            const snp = data.snps[i];
            positions[i * 2] = snp.cumPos;
            positions[i * 2 + 1] = snp.logp;

            // Color based on significance and chromosome
            const color = this.getColor(snp);
            colors[i * 4] = color[0];
            colors[i * 4 + 1] = color[1];
            colors[i * 4 + 2] = color[2];
            colors[i * 4 + 3] = color[3];
        }

        if (this.isWebGPUAvailable) {
            await this.renderWebGPU(positions, colors);
        } else {
            this.renderCanvas2D(data.snps);
        }
    }

    getColor(snp) {
        if (snp.pvalue < this.options.significanceThreshold) {
            return [0.96, 0.62, 0.04, 1.0]; // Orange - significant
        }
        if (snp.pvalue < this.options.suggestiveThreshold) {
            return [0.55, 0.36, 0.96, 0.8]; // Purple - suggestive
        }
        // Alternating chromosome colors
        if (snp.chr % 2 === 0) {
            return [0.23, 0.51, 0.96, 0.7]; // Blue
        }
        return [0.02, 0.71, 0.83, 0.7]; // Cyan
    }

    async renderWebGPU(positions, colors) {
        // Create vertex buffer
        if (this.vertexBuffer) this.vertexBuffer.destroy();
        this.vertexBuffer = this.device.createBuffer({
            size: positions.byteLength,
            usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
        });
        this.device.queue.writeBuffer(this.vertexBuffer, 0, positions);

        // Create color buffer
        if (this.colorBuffer) this.colorBuffer.destroy();
        this.colorBuffer = this.device.createBuffer({
            size: colors.byteLength,
            usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
        });
        this.device.queue.writeBuffer(this.colorBuffer, 0, colors);

        this.draw();
    }

    draw() {
        if (!this.isWebGPUAvailable || !this.vertexBuffer) return;

        // Update uniforms
        const uniforms = new Float32Array([
            this.canvas.width,
            this.canvas.height,
            this.scale,
            this.offsetX,
            this.offsetY,
            this.options.pointSize,
            this.maxLogP,
            this.minX,
            this.maxX,
            0 // padding
        ]);
        this.device.queue.writeBuffer(this.uniformBuffer, 0, uniforms);

        // Create bind group
        const bindGroup = this.device.createBindGroup({
            layout: this.pipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: { buffer: this.uniformBuffer } }
            ]
        });

        // Create command encoder
        const commandEncoder = this.device.createCommandEncoder();

        const renderPass = commandEncoder.beginRenderPass({
            colorAttachments: [{
                view: this.context.getCurrentTexture().createView(),
                clearValue: { r: 0.06, g: 0.09, b: 0.16, a: 1.0 }, // Dark background
                loadOp: 'clear',
                storeOp: 'store'
            }]
        });

        renderPass.setPipeline(this.pipeline);
        renderPass.setBindGroup(0, bindGroup);
        renderPass.setVertexBuffer(0, this.vertexBuffer);
        renderPass.setVertexBuffer(1, this.colorBuffer);
        renderPass.draw(this.numPoints);
        renderPass.end();

        this.device.queue.submit([commandEncoder.finish()]);
    }

    renderCanvas2D(snps) {
        const ctx = this.ctx;
        const dpr = window.devicePixelRatio || 1;

        ctx.fillStyle = '#0f172a';
        ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

        const xScale = (this.canvas.width * this.scale) / (this.maxX - this.minX);
        const yScale = this.canvas.height / this.maxLogP;

        for (const snp of snps) {
            const x = (snp.cumPos - this.minX) * xScale + this.offsetX * dpr;
            const y = this.canvas.height - snp.logp * yScale;

            const color = this.getColor(snp);
            ctx.fillStyle = `rgba(${color[0]*255}, ${color[1]*255}, ${color[2]*255}, ${color[3]})`;
            ctx.beginPath();
            ctx.arc(x, y, this.options.pointSize * dpr, 0, Math.PI * 2);
            ctx.fill();
        }

        // Draw threshold lines
        this.drawThresholds(ctx, yScale, dpr);
    }

    drawThresholds(ctx, yScale, dpr) {
        // Significance threshold
        const sigY = this.canvas.height - (-Math.log10(this.options.significanceThreshold)) * yScale;
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2 * dpr;
        ctx.setLineDash([6 * dpr, 4 * dpr]);
        ctx.beginPath();
        ctx.moveTo(0, sigY);
        ctx.lineTo(this.canvas.width, sigY);
        ctx.stroke();

        // Suggestive threshold
        const sugY = this.canvas.height - (-Math.log10(this.options.suggestiveThreshold)) * yScale;
        ctx.strokeStyle = '#8b5cf6';
        ctx.lineWidth = 1 * dpr;
        ctx.beginPath();
        ctx.moveTo(0, sugY);
        ctx.lineTo(this.canvas.width, sugY);
        ctx.stroke();

        ctx.setLineDash([]);
    }

    setupInteractions() {
        let isDragging = false;
        let lastX = 0;

        this.canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            const mouseX = e.offsetX;

            // Zoom toward mouse position
            const oldScale = this.scale;
            this.scale *= delta;
            this.scale = Math.max(0.1, Math.min(50, this.scale));

            // Adjust offset to zoom toward cursor
            this.offsetX -= (mouseX * (this.scale - oldScale));

            this.draw();
        });

        this.canvas.addEventListener('mousedown', (e) => {
            isDragging = true;
            lastX = e.offsetX;
        });

        this.canvas.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const deltaX = e.offsetX - lastX;
            this.offsetX += deltaX;
            lastX = e.offsetX;
            this.draw();
        });

        this.canvas.addEventListener('mouseup', () => isDragging = false);
        this.canvas.addEventListener('mouseleave', () => isDragging = false);
    }

    resetZoom() {
        this.scale = 1.0;
        this.offsetX = 0;
        this.draw();
    }

    destroy() {
        if (this.vertexBuffer) this.vertexBuffer.destroy();
        if (this.colorBuffer) this.colorBuffer.destroy();
        if (this.uniformBuffer) this.uniformBuffer.destroy();
    }
}

window.WebGPUManhattanPlot = WebGPUManhattanPlot;
