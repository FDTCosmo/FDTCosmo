# DESI DR1 Full-Shape Analysis Data

Parameters and covariance matrices associated with DR1 full-shape analysis cosmology results.

## Download MCMC Chains

The MCMC chains are too large (~4GB) to include in this repository. Download them from the official DESI data release:

**Official Source**: https://data.desi.lbl.gov/doc/releases/dr1/vac/full-shape-cosmo-params/

### Quick Download

```bash
# Create directory structure
mkdir -p mcmc_chains/cobaya

# Download chains (requires DESI data access)
# Visit the link above for download instructions
```

### Expected Directory Structure

After downloading, your directory should look like:

```
data/desi/
├── README_full_shape.md
├── analyze_chains_fdt.py
├── download_chains.sh
├── fdt_analysis_results.json
└── mcmc_chains/
    └── cobaya/
        ├── base/
        ├── base_mnu/
        └── base_w_wa/
```

## References

- DESI Collaboration (2024): https://arxiv.org/abs/2404.03002
