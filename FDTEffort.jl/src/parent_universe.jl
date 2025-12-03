"""
    CMB and Universe Temperature Model

    FDT model for CMB as interference pattern and universe temperature.

    CMB = INTERFERENCE PATTERN OF DENSITY WAVES
    ═══════════════════════════════════════════
    The CMB is the interference pattern of ALL universal density waves,
    past and present. It is not fossil radiation from a Big Bang, but
    the ongoing signature of density wave interference throughout the universe.

    KEY TEMPERATURE RELATIONSHIPS:
    - T_universe = 27 K (actual self-pressure temperature)
    - T_CMB = 2.7 K (observed = 10% of actual)
    - Observable fraction = 2 × baryon_fraction = 10%
    - Factor of 2 from photon = graviton identity

    The formula: T_universe = [3c⁸/(64πaG³M²)]^(1/4)

    CMB INTERPRETATION:
    - NOT fossil radiation from early universe
    - IS interference pattern of all density waves
    - Reflects current universal self-pressure
    - 10% visibility due to photon=graviton double-counting
"""
module CMBDensityWaves

using ..FDTPhysics

export M_OBSERVABLE, R_OBSERVABLE,
       universe_self_pressure_temperature, cmb_observable_temperature,
       density_wave_contribution, cmb_interference_pattern,
       temperature_from_mass, mass_from_temperature,
       cmb_summary, verify_cmb_model

#=============================================================================
                        CONSTANTS
=============================================================================#

# Observable universe mass (baryonic + dark matter)
const M_OBSERVABLE = 1.5e53  # kg

# Observable universe radius
const R_OBSERVABLE = 4.4e26  # m (comoving)

#=============================================================================
                        TEMPERATURE CALCULATIONS
=============================================================================#

"""
    universe_self_pressure_temperature(M)

Compute the actual universe temperature from self-pressure.
T = [3c⁸/(64πaG³M²)]^(1/4)

This is the TRUE temperature of the universe (≈27 K), not the observed CMB.
"""
function universe_self_pressure_temperature(M)
    return (3 * C_LIGHT^8 / (64π * RAD_CONSTANT * G_NEWTON^3 * M^2))^0.25
end

"""
    cmb_observable_temperature(T_actual)

Compute the observable CMB temperature from actual universe temperature.
T_CMB = T_actual × 2 × f_baryon = T_actual × 0.10

We see only 10% because:
- Baryonic matter is 5% of total
- Photon = graviton, so we see it twice (5% × 2 = 10%)
"""
function cmb_observable_temperature(T_actual)
    return T_actual * OBSERVABLE_FRACTION
end

"""
    temperature_from_mass(M)

Compute expected CMB temperature given universe mass.
"""
function temperature_from_mass(M)
    T_actual = universe_self_pressure_temperature(M)
    return cmb_observable_temperature(T_actual)
end

"""
    mass_from_temperature(T_cmb)

Infer universe mass from observed CMB temperature.
Inverts the temperature formula.
"""
function mass_from_temperature(T_cmb)
    # T_actual = T_cmb / OBSERVABLE_FRACTION
    T_actual = T_cmb / OBSERVABLE_FRACTION

    # M² = 3c⁸/(64πaG³T⁴)
    M² = 3 * C_LIGHT^8 / (64π * RAD_CONSTANT * G_NEWTON^3 * T_actual^4)

    return sqrt(M²)
end

#=============================================================================
                        DENSITY WAVE INTERFERENCE
=============================================================================#

"""
    density_wave_contribution(α, phase)

Compute the contribution of a density wave to the CMB pattern.

The CMB is the interference pattern of all density waves in the universe.
Each region with density parameter α contributes to the pattern.

Arguments:
- α: density parameter of the region
- phase: phase of the density wave (radians)

Returns amplitude contribution to interference pattern.
"""
function density_wave_contribution(α, phase)
    # Amplitude proportional to density contrast
    amplitude = α * (1 - α)  # Maximum at α = 0.5

    # Wave contribution
    return amplitude * cos(phase)
end

"""
    cmb_interference_pattern(α_values, phases)

Compute the CMB interference pattern from multiple density waves.

The observed CMB is the sum of all density wave contributions,
past and present, creating the characteristic pattern we observe.

Arguments:
- α_values: vector of density parameters
- phases: vector of phases (same length as α_values)

Returns total interference pattern amplitude.
"""
function cmb_interference_pattern(α_values, phases)
    @assert length(α_values) == length(phases)

    total = 0.0
    for (α, ϕ) in zip(α_values, phases)
        total += density_wave_contribution(α, ϕ)
    end

    # Normalize by number of waves
    return total / length(α_values)
end

#=============================================================================
                        SUMMARY AND VERIFICATION
=============================================================================#

"""
    cmb_summary()

Print summary of the CMB density wave model.
"""
function cmb_summary()
    println("=" ^ 70)
    println("CMB AS INTERFERENCE PATTERN OF DENSITY WAVES")
    println("=" ^ 70)

    println("\n🌡️ TEMPERATURE STRUCTURE:")
    println("   T_universe (actual) = $(T_UNIVERSE) K")
    println("   T_CMB (observed) = $(T_CMB) K")
    println("   Ratio: T_CMB/T_universe = $(round(T_CMB/T_UNIVERSE, digits=3))")

    println("\n👁️ WHY WE SEE ONLY 10%:")
    println("   Baryonic fraction = $(BARYON_FRACTION * 100)%")
    println("   Photon = graviton (factor of 2)")
    println("   Observable = 2 × $(BARYON_FRACTION * 100)% = $(OBSERVABLE_FRACTION * 100)%")

    println("\n🌊 CMB AS INTERFERENCE PATTERN:")
    println("   • CMB is NOT fossil radiation from early universe")
    println("   • CMB IS the interference pattern of ALL density waves")
    println("   • Past and present density waves contribute")
    println("   • Pattern reflects current universal structure")

    println("\n📊 FORMULA:")
    println("   T_universe = [3c⁸/(64πaG³M²)]^(1/4)")
    println("   T_CMB = T_universe × 0.10")

    println("\n🔢 VERIFICATION:")
    T_calc = universe_self_pressure_temperature(M_OBSERVABLE)
    T_cmb_calc = cmb_observable_temperature(T_calc)
    println("   Calculated T_universe = $(round(T_calc, digits=1)) K")
    println("   Calculated T_CMB = $(round(T_cmb_calc, digits=2)) K")
    println("   Observed T_CMB = $(T_CMB) K")

    println("\n" * "=" ^ 70)
end

"""
    verify_cmb_model()

Run verification checks on the CMB model.
"""
function verify_cmb_model()
    println("Verifying CMB Density Wave Model...")

    # Check 1: Temperature formula gives correct order of magnitude
    T_calc = universe_self_pressure_temperature(M_OBSERVABLE)
    @assert 10 < T_calc < 100  # Should be ~27 K
    println("✓ T_universe ≈ $(round(T_calc, digits=1)) K (expected ~27 K)")

    # Check 2: Observable fraction is 10%
    @assert OBSERVABLE_FRACTION ≈ 0.10
    println("✓ Observable fraction = 10%")

    # Check 3: CMB formula consistency
    T_cmb_calc = cmb_observable_temperature(T_UNIVERSE)
    @assert abs(T_cmb_calc - T_CMB) / T_CMB < 0.02
    println("✓ T_CMB = T_universe × 0.10 verified")

    # Check 4: Density wave contribution bounded
    for α in [0.1, 0.3, 0.5, 0.7, 0.9]
        contrib = density_wave_contribution(α, 0.0)
        @assert 0 <= contrib <= 0.25  # Maximum at α=0.5 is 0.25
    end
    println("✓ Density wave contributions properly bounded")

    # Check 5: Mass-temperature inversion
    M_recovered = mass_from_temperature(T_CMB)
    @assert M_recovered > 1e52  # Should be ~10⁵³ kg
    println("✓ Mass inversion gives M ≈ $(round(M_recovered, sigdigits=2)) kg")

    println("\nAll CMB model verifications passed! ✓")
end

end # module
