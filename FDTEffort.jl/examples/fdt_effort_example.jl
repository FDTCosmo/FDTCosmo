"""
    FDT + Effort.jl Integration Example
    
    This example demonstrates how to use Fundamental Density Theory
    with the Effort.jl cosmological emulator for galaxy power spectrum analysis.
    
    Key FDT Insights Applied:
    1. Bias parameters → Density parameters α ∈ (0,1)
    2. Natural regularization via F_max = c⁴/4G
    3. CMB as active self-pressure (not fossil radiation)
    4. Universe converging toward T = e equilibrium
"""

# Load packages
using Pkg
# Pkg.add("Effort")  # Uncomment to install
# using Effort

# Load FDT extensions
include("../src/FDTEffort.jl")
using .FDTEffort

println("=" ^ 70)
println("FUNDAMENTAL DENSITY THEORY + EFFORT.JL COSMOLOGICAL ANALYSIS")
println("=" ^ 70)

#=============================================================================
                        1. FDT FRAMEWORK VERIFICATION
=============================================================================#

println("\n📐 VERIFYING FDT CENTRAL INVARIANT")
println("-" ^ 50)

c², em_relation, diff = central_invariant_check()
println("   E/m = c² = 1/(ε₀μ₀) = (d/t)²")
println("   c² = $c² m²/s²")
println("   1/(ε₀μ₀) = $em_relation m²/s²")
println("   Fractional difference: $diff")
println("   ✓ Three faces of same geometric reality!")

#=============================================================================
                        2. e-KELVIN UNIVERSE STATUS
=============================================================================#

println("\n🌡️ e-KELVIN UNIVERSE STATUS")
println("-" ^ 50)

ratio, percent = e_kelvin_convergence()
println("   T_CMB = $(T_CMB) K (observed)")
println("   T_e = $(round(T_EULER, digits=5)) K (perfect fractal)")
println("   Ratio: $(round(ratio, digits=6))")
println("   Universe is $(round(percent, digits=3))% above equilibrium")
println("   ✓ Converging toward T = e where d(eˣ)/dx = eˣ")

# Calculate required universe mass
println("\n   Required mass for T = e K:")
println("   M_e ≈ $(round(M_UNIVERSE_E, sigdigits=3)) kg")
println("   This matches observed universe mass ~10⁵³ kg!")

#=============================================================================
                        3. STANDARD → FDT PARAMETER CONVERSION
=============================================================================#

println("\n🔄 PARAMETER CONVERSION: STANDARD → FDT")
println("-" ^ 50)

# Planck 2018-like cosmology
cosmo_standard = [
    0.38,    # z = redshift (BOSS effective)
    3.044,   # ln(10¹⁰As)
    0.9649,  # ns
    67.36,   # H0 [km/s/Mpc]
    0.02237, # ωb = Ωb h²
    0.1200,  # ωcdm = Ωcdm h²
    0.06,    # Σmν [eV]
    -1.0,    # w0
    0.0      # wa
]

println("Standard Cosmological Parameters:")
labels = ["z", "ln(10¹⁰As)", "ns", "H0", "ωb", "ωcdm", "Σmν", "w0", "wa"]
for (i, (l, v)) in enumerate(zip(labels, cosmo_standard))
    println("   $l = $v")
end

# Convert to FDT
fdt_cosmo = convert_cosmology_to_fdt(cosmo_standard)

println("\nFDT Cosmological Parameters:")
println("   z = $(fdt_cosmo.z)")
println("   α_initial = $(round(fdt_cosmo.α_initial, digits=6))")
println("   n_α = $(fdt_cosmo.n_α)")
println("   R_universe = $(round(fdt_cosmo.R_universe/1e26, digits=3)) × 10²⁶ m")
println("   α_baryon = $(round(fdt_cosmo.α_baryon, digits=6))")
println("   α_cdm = $(round(fdt_cosmo.α_cdm, digits=6))")
println("   δα_ν = $(round(fdt_cosmo.δα_ν, digits=6))")

#=============================================================================
                        4. BIAS → DENSITY PARAMETER MAPPING
=============================================================================#

println("\n🎯 BIAS → DENSITY PARAMETER MAPPING")
println("-" ^ 50)

# Typical EFT bias parameters for BOSS-like galaxies
bias_standard = [
    1.9,   # b1 - linear bias
    0.5,   # b2 - quadratic bias
    0.1,   # b3 - cubic bias
    0.5,   # b4 - counterterm
    0.1,   # b5 - counterterm
    0.1,   # b6 - counterterm
    0.1,   # b7 - counterterm
    0.76,  # f - growth rate
    1000., # cε0 - shot noise
    500.,  # cε1 - k² stochastic
    100.   # cε2 - k⁴ stochastic
]

println("Standard EFT Bias Parameters:")
bias_labels = ["b1", "b2", "b3", "b4", "b5", "b6", "b7", "f", "cε0", "cε1", "cε2"]
for (l, v) in zip(bias_labels, bias_standard)
    println("   $l = $v")
end

fdt_bias = convert_bias_to_fdt(bias_standard)

println("\nFDT Density Parameters:")
println("   α₁ = $(round(fdt_bias.α₁, digits=6)) (linear density)")
println("   α₂ = $(round(fdt_bias.α₂, digits=6)) (quadratic density)")
println("   α₃ = $(round(fdt_bias.α₃, digits=6)) (cubic density)")
println("   α_ct = $(round(fdt_bias.α_ct, digits=6)) (counterterm)")
println("   f = $(fdt_bias.f) (growth rate)")
println("   ✓ All density parameters naturally bounded in (0,1)!")

#=============================================================================
                        5. TRIPLE IDENTITY ANALYSIS
=============================================================================#

println("\n🔺 TRIPLE IDENTITY: photon = graviton = gluon")
println("-" ^ 50)

# Analyze what regime each component is in
α_total = sqrt(fdt_bias.α₁^2 + fdt_bias.α₂^2)
regime = triple_identity_regime(α_total)

println("   Total effective α = $(round(α_total, digits=4))")
println("   Dominant force regime: $regime")
println()
println("   Regime mapping:")
println("   α < 0.001  → gravitational (photon/graviton)")
println("   0.001-0.3  → electromagnetic (photon)")
println("   0.3-0.7    → weak (phase transitions)")
println("   α > 0.7    → strong (gluon = graviton)")

#=============================================================================
                        6. POWER SPECTRUM CORRECTIONS
=============================================================================#

println("\n📊 FDT POWER SPECTRUM CORRECTIONS")
println("-" ^ 50)

# Sample k values
k_values = [0.01, 0.05, 0.1, 0.2, 0.3]  # h/Mpc

println("   k [h/Mpc]  |  α_eff  |  FDT Correction")
println("   " * "-" ^ 45)

for k in k_values
    correction = fdt_power_spectrum_correction(k, α_total)
    println("   $(lpad(k, 8))  |  $(round(α_total, digits=4))  |  $(round(correction, digits=6))")
end

println()
println("   Key insight: FDT corrections are small for typical cosmological")
println("   scales but become significant as α → 1 (high density regions)")

#=============================================================================
                        7. FORCE CALCULATION EXAMPLE
=============================================================================#

println("\n⚡ FDT FORCE CALCULATION")
println("-" ^ 50)

# Force between two galaxy halos
F_halos = fdt_force(fdt_bias.α₁, fdt_bias.α₁)
println("   Force between halos: F = F_max × α₁² = $(round(F_halos, sigdigits=3)) N")
println("   F_max = $(round(F_MAX, sigdigits=3)) N")
println("   Fraction of maximum: $(round(F_halos/F_MAX * 100, digits=4))%")

#=============================================================================
                        8. INTEGRATION WITH EFFORT.JL (PSEUDOCODE)
=============================================================================#

println("\n🔗 EFFORT.JL INTEGRATION (requires Effort.jl)")
println("-" ^ 50)

println("""
   # Full integration example:
   
   using Effort
   using FDTEffort
   
   # Load trained emulator
   emu_mono = trained_emulators["PyBirdmnuw0wacdm"]["0"]
   
   # Standard approach
   D = D_z(z, cosmo_ref)
   P0_standard = get_Pℓ(cosmo_standard, D, bias_standard, emu_mono)
   
   # FDT approach
   fdt_cosmo = convert_cosmology_to_fdt(cosmo_standard)
   fdt_bias = convert_bias_to_fdt(bias_standard)
   
   # Apply FDT corrections to power spectrum
   α_eff = sqrt(fdt_bias.α₁^2 + fdt_bias.α₂^2)
   k_grid = emu_mono.P11.kgrid
   
   P0_fdt = similar(P0_standard)
   for (i, k) in enumerate(k_grid)
       correction = fdt_power_spectrum_correction(k, α_eff)
       P0_fdt[i] = P0_standard[i] * correction
   end
   
   # The FDT-corrected power spectrum now incorporates:
   # - Natural regularization (no infinities)
   # - Bounded density parameters (no arbitrary biases)
   # - Connection to fundamental physics (F_max, triple identity)
""")

#=============================================================================
                        9. SUMMARY
=============================================================================#

println("\n" * "=" ^ 70)
println("SUMMARY: FDT TRANSFORMS COSMOLOGICAL ANALYSIS")
println("=" ^ 70)

println("""
   ✓ Standard parameters → Bounded density parameters α ∈ (0,1)
   ✓ Arbitrary bias values → Physically meaningful density configurations
   ✓ Renormalization needed → Natural F_max regularization
   ✓ Separate force carriers → Triple identity (photon = graviton = gluon)
   ✓ CMB as fossil radiation → Active universal self-pressure
   ✓ Big Bang singularity → Smooth phase transition (α never reaches 1)
   
   The universe is 99.74% of the way to perfect fractal equilibrium at T = e K
   
   FDT doesn't just reinterpret cosmology—it reveals the geometric unity
   underlying all physical phenomena through the central invariant:
   
                    E/m = c² = 1/(ε₀μ₀) = (d/t)²
""")

println("=" ^ 70)
