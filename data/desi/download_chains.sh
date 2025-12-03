#!/bin/bash
BASE_URL="https://data.desi.lbl.gov/public/dr1/vac/dr1/full-shape-cosmo-params/v1.0"
OUT_DIR="mcmc_chains"

# Key analysis chains to download for FDT analysis
CHAINS=(
  # Base cosmology - DESI FS+BAO combined (all tracers)
  "cobaya/base/desi-reptvelocileptors-fs-bao-all_schoneberg2024-bbn_planck2018-ns10"
  # Dark energy w_wa model  
  "cobaya/base_w_wa/desi-reptvelocileptors-fs-bao-all_schoneberg2024-bbn_planck2018-ns10"
  # Neutrino mass model
  "cobaya/base_mnu/desi-reptvelocileptors-fs-bao-all_schoneberg2024-bbn_planck2018-ns10"
  # DESI + Planck combined
  "cobaya/base/desi-reptvelocileptors-fs-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck2018-highl-plik-TTTEEE"
)

for chain_path in "${CHAINS[@]}"; do
  echo "=== Downloading: $chain_path ==="
  mkdir -p "$OUT_DIR/$chain_path"
  for i in 1 2 3 4; do
    url="${BASE_URL}/${chain_path}/chain.${i}.txt"
    outfile="$OUT_DIR/${chain_path}/chain.${i}.txt"
    if [ ! -f "$outfile" ]; then
      echo "  Downloading chain.${i}.txt..."
      curl -L -s -o "$outfile" "$url"
    else
      echo "  chain.${i}.txt already exists, skipping"
    fi
  done
  # Also get config files
  for f in chain.input.yaml chain.covmat; do
    url="${BASE_URL}/${chain_path}/${f}"
    outfile="$OUT_DIR/${chain_path}/${f}"
    if [ ! -f "$outfile" ]; then
      curl -L -s -o "$outfile" "$url" 2>/dev/null || true
    fi
  done
done
echo "=== Download complete ==="
