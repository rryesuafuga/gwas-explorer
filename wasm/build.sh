#!/bin/bash
# Build script for GWAS WebAssembly module
# Requires: rustup, wasm-pack

set -e

echo "Building GWAS WebAssembly module..."

# Install wasm-pack if not present
if ! command -v wasm-pack &> /dev/null; then
    echo "Installing wasm-pack..."
    curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
fi

# Build the WASM module
wasm-pack build --target web --out-dir pkg

echo "Build complete! Output in wasm/pkg/"
echo ""
echo "To use in your project:"
echo "  import init, { GWASEngine, PRSEngine } from './wasm/pkg/gwas_wasm.js';"
echo "  await init();"
