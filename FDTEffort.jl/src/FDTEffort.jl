"""
    FDTEffort.jl

    Fundamental Density Theory Extension for Cosmological Emulation

    This package extends Effort.jl with FDT physics, providing:

    1. PARAMETER TRANSFORMATION
       - Standard cosmology ↔ FDT density parameters
       - Bias parameters → α ∈ (0,1) naturally bounded

    2. FDT PHYSICS
       - Central invariant: E/m = c² = 1/(ε₀μ₀) = (d/t)²
       - Maximum force: F_max = c⁴/4G (natural regularization)
       - Triple identity: photon = graviton = gluon
       - Universal pressure: P(r) = c⁴/(16πGr²)

    3. COSMOLOGICAL INTERPRETATION
       - CMB as interference pattern of ALL density waves (past & present)
       - T_universe = 27 K (actual), T_CMB = 2.7 K (10% observable)
       - Dark matter = shadow of matter (heavier/denser compressed spacetime)
       - Dark energy = stored potential energy (released when bodies lose mass)
       - No Big Bang singularity (α never reaches 1)

    4. EMULATOR CORRECTIONS
       - Power spectrum corrections from FDT
       - Regularization without renormalization
       - Density-dependent growth factor

    Usage:
    ```julia
    using FDTEffort

    # Check e-Kelvin convergence
    ratio, percent = e_kelvin_convergence()
    println("Universe is \$(percent)% above perfect fractal point")

    # Convert cosmology to FDT
    cosmo_std = [0.5, 3.0, 0.96, 67.0, 0.022, 0.12, 0.06, -1.0, 0.0]
    fdt_cosmo = convert_cosmology_to_fdt(cosmo_std)

    # Check central invariant
    c², em, diff = central_invariant_check()
    println("c² = \$c², 1/(ε₀μ₀) = \$em, diff = \$diff")
    ```

    Author: Manuel Alfaro
    Based on: Effort.jl by Bonici et al. (arXiv:2501.04639)
"""
module FDTEffort

# Include submodules
include("fdt_physics.jl")
include("fdt_emulator.jl")

# Re-export from submodules
using .FDTPhysics
using .FDTEmulator

# Re-export physics constants and functions
export F_MAX, C_LIGHT, G_NEWTON, PLANCK_HBAR, T_CMB, T_UNIVERSE, T_EULER, M_UNIVERSE_E,
       BARYON_FRACTION, OBSERVABLE_FRACTION,
       α_outside, α_inside, schwarzschild_radius, fdt_force, universal_pressure,
       central_invariant_check, triple_identity_regime, e_kelvin_convergence,
       bias_to_alpha, alpha_to_bias, fdt_power_spectrum_correction,
       regularize_with_fmax, fdt_growth_factor, cmb_from_universe_temperature

# Re-export emulator types and functions
export FDTCosmology, FDTBiasParameters,
       convert_cosmology_to_fdt, convert_bias_to_fdt,
       fdt_to_standard_cosmology, fdt_to_standard_bias,
       get_Pℓ_fdt, fdt_stochastic_model, fdt_bias_combination

#=============================================================================
                        PACKAGE INITIALIZATION
=============================================================================#

function __init__()
    # Verify central invariant on load
    c², em, diff = central_invariant_check()
    if diff > 1e-10
        @warn "Central invariant check failed: c²=$c², 1/(ε₀μ₀)=$em"
    end
    
    # Report e-Kelvin status
    ratio, percent = e_kelvin_convergence()
    @info "FDTEffort loaded. Universe at $(round(percent, digits=2))% above T=e equilibrium"
end

#=============================================================================
                        CONVENIENCE FUNCTIONS
=============================================================================#

"""
    fdt_summary()

Print summary of FDT framework and current universe state.
"""
function fdt_summary()
    println("=" ^ 60)
    println("FUNDAMENTAL DENSITY THEORY (FDT) SUMMARY")
    println("=" ^ 60)

    println("\n📐 CENTRAL INVARIANT:")
    c², em, diff = central_invariant_check()
    println("   E/m = c² = 1/(ε₀μ₀) = (d/t)²")
    println("   c² = $(c²) m²/s²")
    println("   1/(ε₀μ₀) = $(em) m²/s²")
    println("   Fractional difference: $(diff)")

    println("\n⚡ MAXIMUM FORCE:")
    println("   F_max = c⁴/4G = $(F_MAX) N")
    println("   This bounds all forces naturally")

    println("\n🌡️ CMB TEMPERATURE:")
    println("   T_universe (actual) = $(T_UNIVERSE) K")
    println("   T_CMB (observed) = $(T_CMB) K")
    println("   Observable fraction = $(OBSERVABLE_FRACTION * 100)%")
    println("   CMB = interference pattern of ALL density waves")

    println("\n🌑 DARK MATTER (Shadow):")
    println("   = Heavier, denser, compressed spacetime")
    println("   = INSEPARABLE from matter (like its shadow)")
    println("   NOT particles - geometry!")

    println("\n⚡ DARK ENERGY (Stored Potential):")
    println("   = Stored potential energy of compressed spacetime")
    println("   = Released when bodies lose mass")
    println("   NOT external force - intrinsic stored energy!")

    println("\n🔺 TRIPLE IDENTITY:")
    println("   photon = graviton = gluon")
    println("   All spin-1 at different α regimes")
    for α in [0.0001, 0.1, 0.5, 0.9]
        regime = triple_identity_regime(α)
        println("   α = $α → $regime")
    end

    println("\n" * "=" ^ 60)
end

"""
    analyze_cosmology(cosmo_standard)

Analyze cosmological parameters through FDT lens.
"""
function analyze_cosmology(cosmo::AbstractVector)
    fdt = convert_cosmology_to_fdt(cosmo)
    
    println("FDT Cosmological Analysis")
    println("-" ^ 40)
    println("Redshift: $(fdt.z)")
    println("α_initial: $(round(fdt.α_initial, digits=4))")
    println("n_α: $(fdt.n_α)")
    println("R_universe: $(round(fdt.R_universe/1e26, digits=2)) × 10²⁶ m")
    println("α_baryon: $(round(fdt.α_baryon, digits=4))")
    println("α_cdm: $(round(fdt.α_cdm, digits=4))")
    println("δα_ν: $(round(fdt.δα_ν, digits=4))")
    
    # Determine dominant regime
    α_total = sqrt(fdt.α_baryon^2 + fdt.α_cdm^2)
    regime = triple_identity_regime(α_total)
    println("\nDominant regime: $regime (α_total = $(round(α_total, digits=4)))")
    
    return fdt
end

end # module FDTEffort
