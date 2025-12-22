# 🧬 GWAS Explorer

An interactive, production-ready Genome-Wide Association Study (GWAS) visualization platform built with modern web technologies. Features real-time Manhattan plots, QQ plots, regional association views, and ML-powered Polygenic Risk Score calculations — all running in the browser.

![GWAS Explorer Preview](https://via.placeholder.com/800x400/0a0f1a/3b82f6?text=GWAS+Explorer)

## ✨ Features

### 📊 Interactive Visualizations
- **Manhattan Plot** — Zoomable, pannable genome-wide visualization with D3.js
- **QQ Plot** — Quantile-quantile analysis with genomic inflation metrics
- **LocusZoom Integration** — Regional association plots with gene annotations
- **Chromosome Preview** — Mini visualization in the hero section

### 🤖 Machine Learning
- **TensorFlow.js Integration** — Browser-based ML inference (no server required!)
- **Polygenic Risk Score Calculator** — Calculate PRS with pre-trained models
- **Privacy-Preserving** — All data stays on your device

### 🚀 Modern Stack
- **Cloudflare Pages** — Global edge deployment with zero cold starts
- **Cloudflare Workers** — Serverless API functions
- **D3.js v7** — Custom, high-performance visualizations
- **TensorFlow.js** — Client-side machine learning

## 🛠️ Tech Stack

| Technology | Purpose |
|------------|---------|
| D3.js v7 | Manhattan plot, QQ plot, custom visualizations |
| TensorFlow.js | Browser-based ML inference for PRS |
| LocusZoom.js | Regional association plots |
| Cloudflare Pages | Static site hosting with global CDN |
| Cloudflare Workers | Serverless API endpoints |

## 📁 Project Structure

```
gwas-explorer/
├── public/                    # Static files (deployed to Cloudflare Pages)
│   └── index.html            # Main HTML file
├── src/
│   ├── css/
│   │   └── styles.css        # All styles (scientific theme)
│   └── js/
│       ├── app.js            # Main application orchestrator
│       ├── gwas-simulator.js # GWAS data simulation
│       ├── manhattan-plot.js # D3.js Manhattan plot
│       ├── qq-plot.js        # D3.js QQ plot
│       ├── locuszoom-integration.js # LocusZoom wrapper
│       └── prs-calculator.js # TensorFlow.js PRS calculator
├── functions/                 # Cloudflare Workers (serverless functions)
│   └── api/
│       └── [[path]].js       # Catch-all API handler
├── wrangler.toml             # Cloudflare configuration
├── package.json              # Node.js dependencies
└── README.md                 # This file
```

## 🚀 Quick Start

### Prerequisites
- Node.js 18+ 
- npm or yarn
- Cloudflare account (free tier works!)

### Local Development

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/gwas-explorer.git
   cd gwas-explorer
   ```

2. **Install dependencies**
   ```bash
   npm install
   ```

3. **Start local development server**
   ```bash
   npm run dev
   ```

4. **Open in browser**
   ```
   http://localhost:8788
   ```

### Deploy to Cloudflare Pages

#### Option 1: Git Integration (Recommended)

1. **Push to GitHub/GitLab**
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin https://github.com/yourusername/gwas-explorer.git
   git push -u origin main
   ```

2. **Connect to Cloudflare Pages**
   - Go to [Cloudflare Dashboard](https://dash.cloudflare.com)
   - Navigate to **Pages** → **Create a project**
   - Select **Connect to Git**
   - Choose your repository
   - Configure build settings:
     - **Build command:** (leave empty)
     - **Build output directory:** `public`
   - Click **Save and Deploy**

3. **Your site is live!**
   ```
   https://gwas-explorer.pages.dev
   ```

#### Option 2: Direct Upload (Wrangler CLI)

1. **Login to Cloudflare**
   ```bash
   npx wrangler login
   ```

2. **Deploy**
   ```bash
   npm run deploy
   ```

3. **Access your site**
   ```
   https://gwas-explorer.pages.dev
   ```

## 📊 Using the Application

### Simulate GWAS Data
Click the **"Simulate GWAS"** button to generate realistic test data with:
- 30,000-70,000 SNPs
- 5-20 significant loci
- Realistic genomic inflation

### Upload Your Data
Upload GWAS summary statistics in JSON or CSV format:

```json
[
  {"snp": "rs12345", "chr": 1, "pos": 1000000, "pval": 1e-8, "beta": 0.15},
  {"snp": "rs67890", "chr": 2, "pos": 2000000, "pval": 0.05, "beta": 0.02}
]
```

Or CSV:
```csv
snp,chr,pos,pval,beta
rs12345,1,1000000,1e-8,0.15
rs67890,2,2000000,0.05,0.02
```

### Calculate Polygenic Risk Score
1. Switch to the **PRS Calculator** section
2. Choose input method:
   - **Manual Entry** — Enter genotype dosages (0, 1, 2)
   - **Simulate** — Generate random genotypes by population
   - **Upload VCF** — Upload genetic data (demo only)
3. Click **Calculate Polygenic Risk Score**

## 🔧 API Endpoints

The Cloudflare Worker provides these API endpoints:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/health` | GET | Health check |
| `/api/simulate` | GET/POST | Generate simulated GWAS data |
| `/api/calculate-prs` | POST | Calculate PRS on server |
| `/api/lookup-snp` | GET | Look up SNP information |
| `/api/gene-annotation` | GET | Get gene annotations for region |

### Example API Usage

```javascript
// Simulate GWAS data
const response = await fetch('/api/simulate?numSnps=10000&numSignificant=5');
const data = await response.json();

// Calculate PRS
const prsResponse = await fetch('/api/calculate-prs', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    genotypes: { rs12345: 1, rs67890: 2 },
    weights: [
      { snp: 'rs12345', weight: 0.15 },
      { snp: 'rs67890', weight: -0.08 }
    ]
  })
});
```

## 🎨 Customization

### Color Schemes
The Manhattan plot supports three color schemes:
- **Alternating** — Classic chromosome alternation (blue/cyan)
- **Gradient** — Viridis color scale by position
- **Significance** — Color by p-value (gold = significant)

### Significance Thresholds
- **5×10⁻⁸** — Genome-wide significance
- **1×10⁻⁵** — Suggestive significance
- **1×10⁻³** — Nominal significance

## 🔬 Technical Details

### Browser Compatibility
- Chrome 90+
- Firefox 88+
- Safari 14+
- Edge 90+

### Performance
- Handles 100,000+ SNPs smoothly
- Uses WebGL acceleration where available
- Lazy loading for large datasets

### Security
- All ML inference runs client-side
- No data transmitted to external servers
- CORS-enabled API endpoints

## 🗺️ Roadmap

- [ ] Real LD data integration
- [ ] Ensembl API integration for gene annotations
- [ ] Multiple phenotype support
- [ ] Export to publication-quality figures
- [ ] VCF file parsing
- [ ] Pre-trained PRS models for common traits

## 📝 License

MIT License — see [LICENSE](LICENSE) for details.

## 🙏 Acknowledgments

- [D3.js](https://d3js.org/) — Data visualization
- [TensorFlow.js](https://www.tensorflow.org/js) — Machine learning
- [LocusZoom.js](https://statgen.github.io/locuszoom/) — Regional association plots
- [Cloudflare](https://www.cloudflare.com/) — Edge hosting

## 📧 Contact

For questions or feedback, open an issue on GitHub or contact the maintainers.

---

**Built with 🧬 for the genomics research community**
