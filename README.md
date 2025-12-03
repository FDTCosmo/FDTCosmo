# FDTCosmo

**Fundamental Density Theory (FDT) Analysis of Cosmological Datasets**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Julia 1.6+](https://img.shields.io/badge/julia-1.6+-purple.svg)](https://julialang.org/)

FDTCosmo provides tools for analyzing cosmological survey data through the framework of Fundamental Density Theory (FDT). The package includes analysis pipelines for BOSS DR12, eBOSS DR16, DESI DR1, and PT Challenge mock datasets.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [FDT Framework](#fdt-framework)
- [Datasets](#datasets)
- [Usage](#usage)
- [Results](#results)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)

## Overview

Fundamental Density Theory (FDT) provides a geometric framework for understanding physics based on the fundamental quantity *density* (ρ = m/V). This repository contains:

- **Python analysis scripts** for processing MCMC chains and power spectrum data
- **Julia package (FDTEffort.jl)** for FDT physics calculations
- **Cosmological datasets** from major galaxy surveys
- **LaTeX reports** summarizing FDT analysis results

### Key Features

- Transform standard cosmological parameters (H₀, Ωₘ, σ₈) to FDT-native quantities (α, R_universe)
- Verify FDT predictions against observational data
- Analyze MCMC chains from DESI full-shape analysis
- Process power spectrum data from BOSS and eBOSS surveys

## Installation

### Prerequisites

- Python 3.8 or higher
- Julia 1.6 or higher (optional, for FDTEffort.jl)

### Python Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/FDTCosmo.git
cd FDTCosmo

# Install Python dependencies
pip install -r requirements.txt
```

### Julia Setup (Optional)

```julia
# In Julia REPL
using Pkg
Pkg.develop(path="FDTEffort.jl")
```

## Quick Start

### Run Full Analysis

```bash
# Analyze all datasets
python analyze_all_datasets_fdt.py
```

This will:
1. Process BOSS DR12, eBOSS DR16, PT Challenge, and DESI DR1 data
2. Convert parameters to FDT framework
3. Generate JSON results (`fdt_all_datasets_results.json`)
4. Create LaTeX summary (`FDT_All_Datasets_Summary.tex`)

### Analyze DESI Chains Only

```bash
cd data/desi
python analyze_chains_fdt.py
```

## FDT Framework

### Central Invariant

The cornerstone of FDT is the central invariant:

```
E/m = c² = 1/(ε₀μ₀) = (d/t)²
```

This identity unifies mechanics, electromagnetism, and geometry.

### Density Parameter α

The fundamental FDT variable is the dimensionless density parameter:

```
α = Rₛ/r ∈ [0, 1)
```

where Rₛ = 2GM/c² is the Schwarzschild radius.

| α Range | Regime |
|---------|--------|
| < 0.001 | Gravitational |
| 0.001 - 0.3 | Electromagnetic |
| 0.3 - 0.7 | Weak |
| > 0.7 | Strong |

### Key FDT Predictions

| Prediction | FDT Value | Observed |
|------------|-----------|----------|
| c² = 1/(ε₀μ₀) | Exact | Δ < 10⁻¹³ |
| T_CMB/T_universe | 0.10 | 0.101 |
| DM/Baryon ratio | 5.4 | 5.26-5.48 |
| Maximum force | c⁴/4G | 3.03×10⁴³ N |

### Parameter Transformations

| Standard | FDT | Physical Meaning |
|----------|-----|------------------|
| H₀ | c/R_universe | Cosmic Schwarzschild scale |
| Ωₘ | α_total | Total density parameter |
| Ωᵦ | α_baryon | Visible density |
| Ωc | α_cdm | Shadow density |
| σ₈ | 10×α_initial | Initial fluctuation amplitude |
| nₛ | n_α | Density spectral index |

## Datasets

### BOSS DR12

- **Samples**: CMASS NGC/SGC (z~0.57), LOWZ NGC/SGC (z~0.32)
- **Data**: Power spectrum multipoles (ℓ=0,2)
- **k range**: 0.01-0.23 h/Mpc

### eBOSS DR16

- **Samples**: LRG (z~0.7), ELG (z~0.85), QSO (z~1.48)
- **Data**: Power spectrum multipoles
- **k range**: 0.01-0.24 h/Mpc

### DESI DR1

- **Models**: Base ΛCDM, Base+mν, w₀wₐCDM
- **Likelihoods**: Full-shape + BAO + Planck 2018
- **Chains**: ~500k-1.1M samples per model
- **Note**: MCMC chains (~4GB) must be downloaded separately from [DESI Data Release](https://data.desi.lbl.gov/doc/releases/dr1/vac/full-shape-cosmo-params/). See `data/desi/README_full_shape.md` for instructions.

### PT Challenge

- **Mocks**: LOWZ, CMASS1, CMASS2
- **Realizations**: 10 per sample
- **Purpose**: Pipeline validation (blinded cosmology)

## Usage

### Python API

```python
from analyze_all_datasets_fdt import convert_to_fdt, central_invariant_check

# Verify central invariant
ci = central_invariant_check()
print(f"c² = {ci['c_squared']:.6e}")
print(f"Verified: {ci['verified']}")

# Convert cosmology to FDT
fdt = convert_to_fdt(
    name="My_Analysis",
    H0=67.4,
    omega_m=0.315,
    omega_b=0.049,
    sigma8=0.811,
    ns=0.965
)

print(f"R_universe = {fdt.R_universe_Gly:.2f} Gly")
print(f"α_total = {fdt.alpha_total:.4f}")
print(f"Regime: {fdt.regime}")
```

### Julia API

```julia
using FDTEffort

# Central invariant
c², em, diff = central_invariant_check()
println("Fractional difference: $diff")

# Density parameter
α = α_outside(r, R_s)
regime = triple_identity_regime(α)

# FDT force
F = fdt_force(α₁, α₂)
```

### Command Line

```bash
# Full analysis with JSON output
python analyze_all_datasets_fdt.py

# DESI-only analysis
python data/desi/analyze_chains_fdt.py

# View results
cat fdt_all_datasets_results.json | python -m json.tool
```

## Results

### Summary Statistics (All Datasets)

- **α_total range**: [0.236, 0.297]
- **α_total mean**: 0.291
- **DM/Baryon ratio**: 5.42 (FDT prediction: 5.4)
- **All datasets in electromagnetic regime**: ✓

### DESI DR1 Results

| Model | H₀ | Ω_m | α_total | DM/Baryon |
|-------|-----|-----|---------|-----------|
| Base+Planck | 68.14±0.40 | 0.305±0.005 | 0.291 | 5.26 |
| Base+BBN | 68.57±0.75 | 0.296±0.009 | 0.283 | 5.31 |
| Base+mν | 68.36±0.42 | 0.302±0.005 | 0.290 | 5.27 |
| w₀wₐCDM | 76.53±3.24 | 0.245±0.021 | 0.236 | 5.48 |

### Dark Energy (w₀wₐCDM)

- **w₀** = -1.52 ± 0.22 (phantom)
- **wₐ** = 0.99 ± 0.61 (evolving)

In FDT: interpreted as release of stored potential energy in the density field.

## Repository Structure

```
FDTCosmo/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── .gitignore                         # Git ignore patterns
│
├── analyze_all_datasets_fdt.py        # Main analysis script
├── fdt_all_datasets_results.json      # Analysis results (JSON)
├── FDT_All_Datasets_Summary.tex       # LaTeX summary
├── FDT_All_Datasets_Summary.pdf       # Compiled PDF
├── FDT_Cosmological_Analysis_Report.tex  # Detailed report
├── FDT_Cosmological_Analysis_Report.pdf  # Compiled PDF
│
├── FDTEffort.jl/                      # Julia package
│   ├── Project.toml
│   ├── README.md
│   ├── src/
│   │   ├── FDTEffort.jl               # Main module
│   │   ├── fdt_physics.jl             # Core physics
│   │   ├── shadow_matter.jl           # Shadow matter model
│   │   ├── parent_universe.jl         # Parent universe model
│   │   └── fdt_emulator.jl            # Emulator
│   ├── test/
│   │   └── runtests.jl
│   └── examples/
│       └── fdt_effort_example.jl
│
└── data/
    ├── boss_dr12/                     # BOSS DR12 data
    │   ├── boss_dr12_2pt.h5           # Power spectrum data
    │   ├── config/                    # Analysis configs
    │   ├── likelihood_config/         # Likelihood configs
    │   └── run_config/                # Run configs
    │
    ├── eboss_dr16/                    # eBOSS DR16 data
    │   ├── eboss_dr16_2pt.h5          # Power spectrum data
    │   └── config/
    │
    ├── desi/                          # DESI DR1 data
    │   ├── README_full_shape.md
    │   ├── analyze_chains_fdt.py      # DESI analysis script
    │   ├── fdt_analysis_results.json  # DESI results
    │   ├── download_chains.sh         # Chain download script
    │   └── mcmc_chains/               # MCMC chains
    │       └── cobaya/
    │           ├── base/
    │           ├── base_mnu/
    │           └── base_w_wa/
    │
    └── pt_challenge/                  # PT Challenge mocks
        ├── camb_params_challenge.ini
        └── R*_*.dat                   # Mock realizations
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use FDTCosmo in your research, please cite:

```bibtex
@software{fdtcosmo2024,
  author = {FDTCosmo Contributors},
  title = {FDTCosmo: Fundamental Density Theory Analysis of Cosmological Datasets},
  year = {2024},
  url = {https://github.com/yourusername/FDTCosmo}
}
```

### Effort.jl Citation

The FDTEffort.jl module builds upon Effort.jl. If you use the Julia components, please also cite:

```bibtex
@article{effortjl2024,
  author = {Bonici, Marco and D'Amico, Guido and Bel, Julien and Carbone, Carmelita},
  title = {Effort.jl: a fast and differentiable emulator for the Effective Field Theory of the Large Scale Structure of the Universe},
  year = {2024},
  url = {https://github.com/CosmologicalEmulators/Effort.jl}
}
```

### Data Citations

This work uses public data from:

- **BOSS DR12**: [Alam et al. (2017)](https://arxiv.org/abs/1607.03155)
- **eBOSS DR16**: [eBOSS Collaboration (2020)](https://arxiv.org/abs/2007.08991)
- **DESI DR1**: [DESI Collaboration (2024)](https://arxiv.org/abs/2404.03002)
- **Planck 2018**: [Planck Collaboration (2018)](https://arxiv.org/abs/1807.06209)

## Acknowledgments

- **[Effort.jl](https://github.com/CosmologicalEmulators/Effort.jl)** by Marco Bonici, Guido D'Amico, Julien Bel, and Carmelita Carbone - the fast and differentiable EFT emulator that FDTEffort.jl builds upon
- SDSS/BOSS/eBOSS collaborations for public galaxy survey data
- DESI collaboration for full-shape analysis chains
- Planck collaboration for CMB data products
