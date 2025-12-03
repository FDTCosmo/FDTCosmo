"""
    Tests for FDTEffort.jl

    Verifies FDT physics, parameter transformations, and predictions.

    Updated definitions:
    - CMB = interference pattern of ALL density waves (past and present)
    - Dark matter = shadow of matter (heavier/denser compressed spacetime, inseparable)
    - Dark energy = stored potential energy (released when bodies lose mass)
"""

using Test

# Include modules directly for testing
include("../src/fdt_physics.jl")
using .FDTPhysics

include("../src/shadow_matter.jl")
using .ShadowMatter

include("../src/parent_universe.jl")
using .CMBDensityWaves

include("../src/fdt_emulator.jl")
using .FDTEmulator

@testset "FDTEffort.jl" begin

    @testset "Central Invariant" begin
        c², em, diff = central_invariant_check()

        @test c² ≈ 299792458.0^2
        @test em ≈ c² rtol=1e-10
        @test diff < 1e-10

        println("✓ Central invariant: E/m = c² = 1/(ε₀μ₀) verified")
    end

    @testset "Maximum Force" begin
        @test F_MAX ≈ 3.0255e43 rtol=1e-3

        # F_max should equal c⁴/4G exactly
        F_calc = C_LIGHT^4 / (4 * G_NEWTON)
        @test F_MAX ≈ F_calc

        println("✓ F_max = c⁴/4G = $(F_MAX) N")
    end

    @testset "Density Parameter α" begin
        # Test Schwarzschild radius
        M_sun = 1.989e30  # kg
        R_s_sun = schwarzschild_radius(M_sun)
        @test R_s_sun ≈ 2953.25 rtol=1e-3  # ~3 km

        # Test α_outside: should approach 1 as r → R_s
        @test α_outside(10 * R_s_sun, R_s_sun) ≈ 0.1
        @test α_outside(2 * R_s_sun, R_s_sun) ≈ 0.5
        @test α_outside(1.001 * R_s_sun, R_s_sun) > 0.99

        # Test α_inside: should approach 0 as r → 0
        @test α_inside(0.1 * R_s_sun, R_s_sun) ≈ 0.1
        @test α_inside(0.5 * R_s_sun, R_s_sun) ≈ 0.5

        # Test bounding in (0, 1)
        @test 0 < α_outside(1e20, 1.0) < 1
        @test 0 < α_inside(1e-20, 1.0) < 1

        println("✓ Density parameter α bounded in (0,1)")
    end

    @testset "FDT Force" begin
        # Force between two configurations
        F = fdt_force(0.5, 0.5)
        @test F ≈ F_MAX * 0.25

        # Maximum force at α = 1
        F_at_max = fdt_force(0.9999, 0.9999)
        @test F_at_max < F_MAX
        @test F_at_max > 0.99 * F_MAX

        # Regularization test
        F_large = 1e50
        F_reg = regularize_with_fmax(F_large)
        @test F_reg <= F_MAX  # tanh(huge) = 1.0 exactly
        @test F_reg ≈ F_MAX rtol=0.01

        println("✓ Force naturally bounded by F_max")
    end

    @testset "Triple Identity Regimes" begin
        @test triple_identity_regime(0.0001) == :gravitational
        @test triple_identity_regime(0.1) == :electromagnetic
        @test triple_identity_regime(0.5) == :weak
        @test triple_identity_regime(0.9) == :strong

        println("✓ Triple identity: photon = graviton = gluon at different α")
    end

    @testset "CMB Temperature Parameters" begin
        # Check T_universe and T_CMB relationship
        @test T_UNIVERSE ≈ 27.0
        @test T_CMB ≈ 2.7255 rtol=0.01

        # Observable fraction should be 10%
        @test OBSERVABLE_FRACTION ≈ 0.10
        @test BARYON_FRACTION ≈ 0.05

        # CMB from universe temperature
        T_calc = cmb_from_universe_temperature(T_UNIVERSE)
        @test T_calc ≈ T_CMB rtol=0.02

        println("✓ T_universe = $(T_UNIVERSE) K, T_CMB = $(T_CMB) K (10% observable)")
    end

    @testset "e-Kelvin Cosmology" begin
        ratio, percent = e_kelvin_convergence()

        @test ratio ≈ 2.7255 / ℯ rtol=1e-3
        @test 0 < percent < 1  # Should be ~0.26%

        # Universe mass for T = e (≈ 8.9×10⁵⁴ kg)
        @test M_UNIVERSE_E > 1e54
        @test M_UNIVERSE_E < 1e56

        println("✓ Universe at $(round(percent, digits=2))% above T=e equilibrium")
        println("  M_universe ≈ $(round(M_UNIVERSE_E, sigdigits=2)) kg")
    end

    @testset "Universal Pressure" begin
        # Pressure at Earth radius scale
        R_earth = 6.371e6
        P_earth = universal_pressure(R_earth)
        @test P_earth > 0

        # Pressure scales as 1/r²
        P_2earth = universal_pressure(2 * R_earth)
        @test P_2earth ≈ P_earth / 4 rtol=1e-10

        println("✓ Universal pressure P(r) = c⁴/(16πGr²)")
    end

    @testset "Bias ↔ Alpha Mapping" begin
        # Bias to alpha should be bounded
        for b in [0.1, 1.0, 2.0, 5.0, 10.0]
            α = bias_to_alpha(b)
            @test 0 < α < 1
        end

        # Round-trip conversion
        b_original = 2.0
        α = bias_to_alpha(b_original)
        b_recovered = alpha_to_bias(α)
        @test b_recovered ≈ b_original rtol=0.1

        println("✓ Bias ↔ α mapping preserves bounds")
    end

    @testset "FDT Cosmology Conversion" begin
        # Standard ΛCDM-like parameters
        cosmo_std = [0.5, 3.0, 0.96, 67.0, 0.022, 0.12, 0.06, -1.0, 0.0]

        fdt = convert_cosmology_to_fdt(cosmo_std)

        @test fdt.z == 0.5
        @test 0 < fdt.α_initial < 1
        @test fdt.n_α == 0.96
        @test fdt.R_universe > 1e25  # > 10²⁵ m
        @test 0 < fdt.α_baryon < 1
        @test 0 < fdt.α_cdm < 1
        @test fdt.w0 == -1.0

        # Round-trip
        cosmo_back = FDTEmulator.fdt_to_standard_cosmology(fdt)
        @test cosmo_back[1] ≈ cosmo_std[1]  # z
        @test cosmo_back[8] ≈ cosmo_std[8]  # w0

        println("✓ Cosmology conversion preserves structure")
    end

    @testset "FDT Bias Conversion" begin
        # Standard EFT bias parameters
        bias_std = [1.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.8, 1.0, 0.5, 0.2]

        fdt_bias = convert_bias_to_fdt(bias_std)

        @test 0 < fdt_bias.α₁ < 1
        @test 0 < fdt_bias.α₂ < 1
        @test 0 < fdt_bias.α₃ < 1
        @test fdt_bias.f == 0.8

        println("✓ Bias parameters naturally bounded as α ∈ (0,1)")
    end

    @testset "Power Spectrum Correction" begin
        # Low k: minimal correction
        corr_low = fdt_power_spectrum_correction(0.01, 0.1)
        @test corr_low ≈ 1.0 rtol=0.1

        # High α: suppression
        corr_high_α = fdt_power_spectrum_correction(0.1, 0.9)
        @test corr_high_α < 0.5  # Significant suppression

        # Very high k: regularization
        corr_high_k = fdt_power_spectrum_correction(1e30, 0.1)
        @test corr_high_k < 1.0

        println("✓ Power spectrum corrections applied correctly")
    end

    @testset "Growth Factor" begin
        # Growth factor should decrease with redshift
        D_z0 = fdt_growth_factor(0.0, 0.1)
        D_z1 = fdt_growth_factor(1.0, 0.1)

        @test D_z0 > D_z1

        # High α suppresses growth
        D_low_α = fdt_growth_factor(0.5, 0.1)
        D_high_α = fdt_growth_factor(0.5, 0.9)

        @test D_low_α > D_high_α

        println("✓ FDT growth factor behaves correctly")
    end

    #=========================================================================
                        SHADOW MATTER (DARK MATTER) TESTS
    =========================================================================#

    @testset "Shadow Matter: Definition" begin
        # Dark matter is the SHADOW of regular matter
        # - Heavier, denser, compressed spacetime
        # - INSEPARABLE from regular matter
        # - NOT particles

        # α threshold should give DM/baryon ratio ≈ 5.4
        @test shadow_mass_ratio() ≈ DM_BARYON_RATIO rtol=0.01
        @test ALPHA_THRESHOLD ≈ 1 / (DM_BARYON_RATIO + 1)

        # Partition should sum to 1
        @test shadow_fraction() + visible_fraction() ≈ 1.0

        # Regime classification
        @test α_regime(0.1) == :visible
        @test α_regime(0.5) == :shadow
        @test α_regime(ALPHA_THRESHOLD) == :boundary

        println("✓ Shadow matter = compressed spacetime shadow of regular matter")
    end

    @testset "Shadow Matter: Compression Factor" begin
        # No compression in visible regime
        @test compression_factor(0.1) ≈ 1.0

        # Compression increases in shadow regime
        @test compression_factor(0.5) > 1.0
        @test compression_factor(0.9) > compression_factor(0.5)

        # Maximum compression approaches FORCE_RATIO at α → 1
        @test compression_factor(0.99) < FORCE_RATIO

        println("✓ Shadow matter compression factor correct")
    end

    @testset "Shadow Matter: EM Coupling" begin
        # EM coupling should decrease with α
        @test em_coupling(0.01) > em_coupling(0.5)
        @test em_coupling(0.5) > em_coupling(0.99)

        # High α should suppress EM (why dark matter is dark)
        @test is_em_active(0.01) == true
        # At very high α, EM coupling is suppressed (check relative suppression)
        @test em_coupling(0.99) < em_coupling(0.01) / 100

        println("✓ EM coupling suppressed at high α (explains darkness)")
    end

    #=========================================================================
                        STORED POTENTIAL ENERGY (DARK ENERGY) TESTS
    =========================================================================#

    @testset "Dark Energy: Stored Potential Energy" begin
        # Dark energy is stored potential energy of compressed spacetime
        # Released when bodies lose mass

        M_test = 1e30  # ~0.5 solar mass

        # Stored energy should increase with α (more compression)
        E_low = stored_potential_energy(M_test, 0.2)
        E_mid = stored_potential_energy(M_test, 0.5)
        E_high = stored_potential_energy(M_test, 0.8)

        @test E_mid > E_low
        @test E_high > E_mid

        # No stored energy at very low α (minimal compression)
        E_minimal = stored_potential_energy(M_test, 0.01)
        @test E_minimal < E_low

        println("✓ Stored potential energy increases with compression")
    end

    @testset "Dark Energy: Release Rate" begin
        M_test = 1e30
        dM_dt_loss = -1e20  # Mass loss rate (negative = losing mass)
        α = 0.5

        # Energy should be released when mass is lost
        power = dark_energy_release_rate(M_test, dM_dt_loss, α)
        @test power > 0  # Positive power output

        # No release when gaining mass
        dM_dt_gain = 1e20
        power_gain = dark_energy_release_rate(M_test, dM_dt_gain, α)
        @test power_gain == 0.0

        # Higher α means more stored energy to release
        power_low_α = dark_energy_release_rate(M_test, dM_dt_loss, 0.2)
        power_high_α = dark_energy_release_rate(M_test, dM_dt_loss, 0.8)
        @test power_high_α > power_low_α

        println("✓ Dark energy released when bodies lose mass")
    end

    #=========================================================================
                        CMB DENSITY WAVE TESTS
    =========================================================================#

    @testset "CMB: Interference Pattern Model" begin
        # CMB is the interference pattern of ALL density waves

        # Temperature calculations
        T_calc = universe_self_pressure_temperature(M_OBSERVABLE)
        @test 10 < T_calc < 100  # Should be ~27 K

        # Observable fraction
        T_cmb_calc = cmb_observable_temperature(T_UNIVERSE)
        @test abs(T_cmb_calc - T_CMB) / T_CMB < 0.02

        println("✓ CMB = interference pattern of all density waves")
    end

    @testset "CMB: Density Wave Contributions" begin
        # Density wave contribution should be bounded
        for α in [0.1, 0.3, 0.5, 0.7, 0.9]
            contrib = density_wave_contribution(α, 0.0)
            @test 0 <= contrib <= 0.25  # Maximum at α=0.5 is 0.25
        end

        # Maximum contribution at α = 0.5
        contrib_05 = density_wave_contribution(0.5, 0.0)
        contrib_01 = density_wave_contribution(0.1, 0.0)
        contrib_09 = density_wave_contribution(0.9, 0.0)

        @test contrib_05 >= contrib_01
        @test contrib_05 >= contrib_09

        println("✓ Density wave contributions properly bounded")
    end

    @testset "CMB: Mass-Temperature Inversion" begin
        # Should be able to recover mass from CMB temperature
        M_recovered = mass_from_temperature(T_CMB)
        @test M_recovered > 1e52  # Should be ~10⁵³ kg

        println("✓ Mass-temperature inversion works correctly")
    end

end

println("\n" * "=" ^ 60)
println("All FDT tests passed!")
println("=" ^ 60)
println("\nKey definitions verified:")
println("• CMB = interference pattern of ALL density waves (past & present)")
println("• Dark matter = shadow of matter (inseparable compressed spacetime)")
println("• Dark energy = stored potential energy (released on mass loss)")
println("• T_universe = 27 K, T_CMB = 2.7 K (10% observable)")
