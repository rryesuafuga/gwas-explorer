/**
 * GWAS Explorer API - Cloudflare Worker
 * Handles server-side processing for GWAS data
 */

export async function onRequest(context) {
    const { request, env } = context;
    const url = new URL(request.url);
    const path = url.pathname.replace('/api', '');
    
    // CORS headers
    const corsHeaders = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
        'Access-Control-Allow-Headers': 'Content-Type',
    };
    
    // Handle CORS preflight
    if (request.method === 'OPTIONS') {
        return new Response(null, { headers: corsHeaders });
    }
    
    try {
        // Route handling
        switch (path) {
            case '/health':
                return jsonResponse({ status: 'ok', timestamp: Date.now() }, corsHeaders);
            
            case '/simulate':
                return handleSimulate(request, corsHeaders);
            
            case '/calculate-prs':
                return handleCalculatePRS(request, corsHeaders);
            
            case '/lookup-snp':
                return handleLookupSNP(request, corsHeaders);
            
            case '/gene-annotation':
                return handleGeneAnnotation(request, corsHeaders);
            
            default:
                return jsonResponse({ 
                    error: 'Not found',
                    availableEndpoints: [
                        '/api/health',
                        '/api/simulate',
                        '/api/calculate-prs',
                        '/api/lookup-snp',
                        '/api/gene-annotation'
                    ]
                }, corsHeaders, 404);
        }
    } catch (error) {
        console.error('API Error:', error);
        return jsonResponse({ error: error.message }, corsHeaders, 500);
    }
}

/**
 * Generate simulated GWAS data
 */
async function handleSimulate(request, corsHeaders) {
    const params = await getRequestParams(request);
    
    const numSnps = Math.min(parseInt(params.numSnps) || 10000, 100000);
    const numSignificant = Math.min(parseInt(params.numSignificant) || 10, 50);
    const genomicInflation = parseFloat(params.genomicInflation) || 1.02;
    
    console.log(`Simulating GWAS: ${numSnps} SNPs, ${numSignificant} significant loci`);
    
    // Chromosome lengths (GRCh38)
    const chrLengths = {
        1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555,
        5: 181538259, 6: 170805979, 7: 159345973, 8: 145138636,
        9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
        13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345,
        17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
        21: 46709983, 22: 50818468
    };
    
    const snps = [];
    const significanceThreshold = 5e-8;
    
    // Generate significant loci first
    for (let i = 0; i < numSignificant; i++) {
        const chr = Math.floor(Math.random() * 22) + 1;
        const pos = Math.floor(Math.random() * chrLengths[chr]);
        const logP = -Math.log10(significanceThreshold) + Math.random() * 15;
        
        snps.push({
            snp: `rs${Math.floor(Math.random() * 900000000) + 100000000}`,
            chr,
            pos,
            pval: Math.pow(10, -logP),
            logP,
            beta: (Math.random() - 0.5) * 0.6,
            maf: 0.01 + Math.random() * 0.49,
            isSignificant: true,
            isLead: true
        });
    }
    
    // Generate null SNPs
    for (let i = 0; i < numSnps - numSignificant; i++) {
        const chr = Math.floor(Math.random() * 22) + 1;
        const pos = Math.floor(Math.random() * chrLengths[chr]);
        
        // Chi-squared distribution for null p-values
        const u1 = Math.random();
        const u2 = Math.random();
        const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        const chiSq = z * z * genomicInflation;
        
        // Convert to p-value (approximation)
        const pval = Math.exp(-chiSq / 2);
        const logP = -Math.log10(Math.max(pval, 1e-300));
        
        snps.push({
            snp: `rs${Math.floor(Math.random() * 900000000) + 100000000}`,
            chr,
            pos,
            pval: Math.max(pval, 1e-300),
            logP: Math.min(logP, 300),
            beta: (Math.random() - 0.5) * 0.1,
            maf: 0.01 + Math.random() * 0.49,
            isSignificant: false,
            isLead: false
        });
    }
    
    // Sort by chromosome and position
    snps.sort((a, b) => a.chr !== b.chr ? a.chr - b.chr : a.pos - b.pos);
    
    // Add cumulative positions
    let cumPos = 0;
    let currentChr = 0;
    for (const snp of snps) {
        if (snp.chr !== currentChr) {
            if (currentChr > 0) cumPos += chrLengths[currentChr];
            currentChr = snp.chr;
        }
        snp.cumPos = cumPos + snp.pos;
    }
    
    return jsonResponse({
        success: true,
        count: snps.length,
        significantCount: snps.filter(s => s.isSignificant).length,
        data: snps
    }, corsHeaders);
}

/**
 * Calculate Polygenic Risk Score
 */
async function handleCalculatePRS(request, corsHeaders) {
    if (request.method !== 'POST') {
        return jsonResponse({ error: 'POST method required' }, corsHeaders, 405);
    }
    
    const body = await request.json();
    const { genotypes, weights } = body;
    
    if (!genotypes || !weights) {
        return jsonResponse({ 
            error: 'Missing required fields: genotypes, weights' 
        }, corsHeaders, 400);
    }
    
    // Calculate weighted sum
    let score = 0;
    const contributions = {};
    
    for (const variant of weights) {
        const dosage = genotypes[variant.snp] || 0;
        const contribution = dosage * variant.weight;
        score += contribution;
        contributions[variant.snp] = contribution;
    }
    
    // Normalize and calculate percentile (simplified)
    // In practice, you'd use population-specific distributions
    const normalizedScore = score;
    
    // Approximate percentile using normal distribution
    // Assuming mean=0, sd=1 for demonstration
    const z = normalizedScore;
    const percentile = 100 * (0.5 * (1 + erf(z / Math.sqrt(2))));
    
    return jsonResponse({
        success: true,
        score: normalizedScore,
        percentile: Math.min(99.9, Math.max(0.1, percentile)),
        contributions,
        riskCategory: percentile > 80 ? 'elevated' : percentile > 50 ? 'average' : 'low'
    }, corsHeaders);
}

/**
 * Look up SNP information
 */
async function handleLookupSNP(request, corsHeaders) {
    const params = await getRequestParams(request);
    const snpId = params.snp || params.rsid;
    
    if (!snpId) {
        return jsonResponse({ error: 'Missing snp parameter' }, corsHeaders, 400);
    }
    
    // In production, this would query a database or external API
    // For demo, return mock data
    const mockData = {
        snp: snpId,
        chromosome: Math.floor(Math.random() * 22) + 1,
        position: Math.floor(Math.random() * 200000000),
        alleles: ['A', 'G'],
        maf: {
            global: 0.25,
            AFR: 0.35,
            EUR: 0.20,
            EAS: 0.28
        },
        genes: ['NEARBY_GENE'],
        consequence: 'intron_variant',
        publications: 5
    };
    
    return jsonResponse({
        success: true,
        data: mockData
    }, corsHeaders);
}

/**
 * Get gene annotations for a region
 */
async function handleGeneAnnotation(request, corsHeaders) {
    const params = await getRequestParams(request);
    const chr = parseInt(params.chr);
    const start = parseInt(params.start);
    const end = parseInt(params.end);
    
    if (!chr || !start || !end) {
        return jsonResponse({ 
            error: 'Missing required parameters: chr, start, end' 
        }, corsHeaders, 400);
    }
    
    // Mock gene data - in production, query Ensembl or similar
    const genes = [
        { name: 'GENE1', start: start + 10000, end: start + 50000, strand: '+' },
        { name: 'GENE2', start: start + 60000, end: start + 120000, strand: '-' },
    ].filter(g => g.start < end && g.end > start);
    
    return jsonResponse({
        success: true,
        region: { chr, start, end },
        genes
    }, corsHeaders);
}

/**
 * Helper: Get request parameters from GET or POST
 */
async function getRequestParams(request) {
    if (request.method === 'POST') {
        try {
            return await request.json();
        } catch {
            return {};
        }
    }
    
    const url = new URL(request.url);
    const params = {};
    url.searchParams.forEach((value, key) => {
        params[key] = value;
    });
    return params;
}

/**
 * Helper: Create JSON response
 */
function jsonResponse(data, corsHeaders, status = 200) {
    return new Response(JSON.stringify(data), {
        status,
        headers: {
            'Content-Type': 'application/json',
            ...corsHeaders
        }
    });
}

/**
 * Error function approximation for percentile calculation
 */
function erf(x) {
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;
    
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    
    return sign * y;
}
