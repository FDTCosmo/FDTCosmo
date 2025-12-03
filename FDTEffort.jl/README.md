# FDTEffort.jl - Fundamental Density Theory Integration with Effort.jl

A Julia package implementing Fundamental Density Theory (FDT) physics for cosmological emulators.

## Acknowledgment

This package builds upon and integrates with [**Effort.jl**](https://github.com/CosmologicalEmulators/Effort.jl) by Marco Bonici et al. - a fast and differentiable emulator for the Effective Field Theory of the Large Scale Structure of the Universe.

If you use this package, please also cite the original Effort.jl paper:

> M. Bonici, G. D'Amico, J. Bel, C. Carbone, "Effort.jl: a fast and differentiable emulator for the Effective Field Theory of the Large Scale Structure of the Universe"

## Overview

FDT proposes that reality is pure geometry, with the central invariant:

```
E/m = c² = 1/(ε₀μ₀) = (d/t)²
```

This package implements:
- **Core FDT Physics**: Maximum force principle, density parameter α, triple identity
- **Shadow Matter**: Dark matter as inseparable geometric shadow of regular matter
- **Stored Potential Energy**: Dark energy as stored compression energy
- **CMB as Density Waves**: The CMB as interference pattern of all universal density waves
- **Cosmological Integration**: Parameter transformations for Effort.jl emulators

## Key Physical Insights

### CMB Temperature (from Triple Identity Paper)

The CMB is the **interference pattern of ALL universal density waves, past and present**.

| Parameter | Value | Description |
|-----------|-------|-------------|
| T_universe | 27 K | Actual self-pressure temperature |
| T_CMB | 2.7 K | Observed (10% of actual) |
| Baryon fraction | 5% | Visible matter |
| Observable fraction | 10% | 2 × baryon (photon=graviton) |

**Why we see only 10%**: The photon = graviton identity means we see baryonic matter twice (once as photon, once as graviton channel). 5% × 2 = 10%.

Formula: `T_universe = [3c⁸/(64πaG³M²)]^(1/4) ≈ 27 K`

### Dark Matter (Shadow Matter)

**Dark matter is simply heavier, denser, compressed spacetime.**

Key properties:
- **INSEPARABLE from regular matter** - like its shadow
- Heavier and denser than visible matter
- Compressed spacetime geometry (NOT particles)
- Appears at α > α_threshold ≈ 0.156
- DM/baryon ratio ≈ 5.4 from α-space partition

**No DM particles will ever be found** - dark matter is the geometric shadow cast by regular matter.

### Dark Energy (Stored Potential Energy)

**Dark energy is the stored potential energy of compressed spacetime and matter.**

Key properties:
- Stored in the compression of spacetime geometry
- **Released when bodies lose mass** (radiation, decay, etc.)
- Drives cosmic acceleration as compressed regions relax
- NOT an external force - intrinsic stored energy

### Force Hierarchy

- `F_max = c⁴/4G` - Maximum force bound (~3.03×10⁴³ N)
- `F_Planck = c⁴/G` - Planck force (4 × F_max)
- This creates the visible/shadow partition at α_threshold

## Installation

```julia
using Pkg
Pkg.add(url="path/to/FDTEffort.jl")
```

Or for development:
```julia
include("src/fdt_physics.jl")
using .FDTPhysics

include("src/shadow_matter.jl")
using .ShadowMatter

include("src/parent_universe.jl")
using .CMBDensityWaves
```

## Quick Start

```julia
using FDTPhysics, ShadowMatter, CMBDensityWaves

# CMB Temperature relationship
println("T_universe = $(T_UNIVERSE) K")
println("T_CMB = $(T_CMB) K")
println("Observable = $(OBSERVABLE_FRACTION * 100)%")

# Shadow matter (dark matter)
ratio = shadow_mass_ratio()  # ≈ 5.4
regime = α_regime(0.3)  # :shadow

# Stored potential energy (dark energy)
E_stored = stored_potential_energy(1e30, 0.5)

# Cosmic budget summary
cosmic_budget_summary()

# CMB density wave model
cmb_summary()
```

## Module Reference

### FDTPhysics

Core constants and functions:

| Constant | Value | Description |
|----------|-------|-------------|
| `T_UNIVERSE` | 27.0 K | Actual universe temperature |
| `T_CMB` | 2.7255 K | Observed CMB temperature |
| `BARYON_FRACTION` | 0.05 | Baryonic matter fraction |
| `OBSERVABLE_FRACTION` | 0.10 | Observable (photon=graviton) |
| `F_MAX` | 3.03×10⁴³ N | Maximum force |

| Function | Description |
|----------|-------------|
| `schwarzschild_radius(M)` | R_s = 2GM/c² |
| `α_outside(r, R_s)` | α = R_s/r for r > R_s |
| `cmb_from_universe_temperature(T)` | T_CMB = T × 0.10 |
| `triple_identity_regime(α)` | Maps α to force regime |

### ShadowMatter

Dark matter as geometric shadow:

| Function | Description |
|----------|-------------|
| `α_regime(α)` | :visible, :shadow, or :boundary |
| `shadow_mass_ratio()` | ≈ 5.4 |
| `compression_factor(α)` | Spacetime compression at α |
| `stored_potential_energy(M, α)` | Dark energy stored in compression |
| `dark_energy_release_rate(M, dM_dt, α)` | Energy released on mass loss |
| `cosmic_budget_summary()` | Print cosmic energy budget |

### CMBDensityWaves

CMB as interference pattern:

| Function | Description |
|----------|-------------|
| `universe_self_pressure_temperature(M)` | T from self-pressure |
| `cmb_observable_temperature(T)` | Observable T (10% of actual) |
| `density_wave_contribution(α, phase)` | Single wave contribution |
| `cmb_interference_pattern(αs, phases)` | Total interference |
| `cmb_summary()` | Print CMB model summary |

## Physical Predictions

### Dark Matter (Shadow Matter)
1. **No particles will be found** - dark matter is compressed spacetime
2. **Dark matter is INSEPARABLE from regular matter** (like its shadow)
3. **DM/baryon ≈ 5.4 universal** across all scales
4. **Halos are extended α-fields**, not particle clouds

### Dark Energy (Stored Potential)
1. **Released when bodies lose mass** (radiation, decay)
2. **Drives cosmic acceleration** as compression relaxes
3. **Stored in spacetime geometry** itself

### CMB
1. **NOT fossil radiation** from early universe
2. **IS interference pattern** of ALL density waves (past and present)
3. **T_universe = 27 K** (actual), **T_CMB = 2.7 K** (10% observable)
4. **Pattern reflects current structure**, not ancient history

## Testing

```julia
include("test/runtests.jl")
```

All core physics tests verify the updated definitions.

## References

- Alfaro, M. (2023). "Foundations of Fundamental Density Theory (FDT)"
- Alfaro, M. (2024). "The Photon Is the Graviton Is the Gluon: Complete Force Unification"
- Alfaro, M. (2024). "Shadow Matter: Dark Matter as Inseparable Geometry"
- Alfaro, M. (2024). "CMB as Density Wave Interference Pattern"

## License

MIT License - See LICENSE file

## Citation

If you use FDTEffort.jl in your research, please cite both this package and the original Effort.jl:

```bibtex
@software{FDTEffort,
  author = {Alfaro, Manuel},
  title = {FDTEffort.jl: Fundamental Density Theory for Cosmological Emulators},
  year = {2024},
  url = {https://github.com/user/FDTEffort.jl}
}

@article{Effort.jl,
  author = {Bonici, Marco and D'Amico, Guido and Bel, Julien and Carbone, Carmelita},
  title = {Effort.jl: a fast and differentiable emulator for the Effective Field Theory of the Large Scale Structure of the Universe},
  year = {2024},
  url = {https://github.com/CosmologicalEmulators/Effort.jl}
}
```
