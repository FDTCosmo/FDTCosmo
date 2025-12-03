"""
    FDT Physics Module

    Core constants, equations, and physical relationships from Fundamental Density Theory.
    Central Invariant: E/m = c² = 1/(ε₀μ₀) = (d/t)²

    CMB INTERPRETATION:
    The CMB is the interference pattern of all universal density waves, past and present.
    - T_universe = 27 K (actual self-pressure temperature)
    - T_CMB = 2.7 K (10% observable via photon=graviton identity)
    - We see 5% baryonic matter twice (as photon and graviton channels)
"""
module FDTPhysics

export F_MAX, C_LIGHT, G_NEWTON, PLANCK_HBAR, T_CMB, T_UNIVERSE, T_EULER, M_UNIVERSE_E,
       RAD_CONSTANT, STEFAN_BOLTZMANN, BARYON_FRACTION, OBSERVABLE_FRACTION,
       α_outside, α_inside, schwarzschild_radius, fdt_force, universal_pressure,
       central_invariant_check, triple_identity_regime, e_kelvin_convergence,
       bias_to_alpha, alpha_to_bias, fdt_power_spectrum_correction,
       regularize_with_fmax, fdt_growth_factor, cmb_from_universe_temperature

# Fundamental Constants
const C_LIGHT = 299792458.0
const G_NEWTON = 6.67430e-11
const PLANCK_HBAR = 1.054571817e-34
const EPSILON_0 = 8.8541878128e-12
const MU_0 = 1.25663706212e-6
const K_BOLTZMANN = 1.380649e-23
const STEFAN_BOLTZMANN = 5.670374419e-8
const RAD_CONSTANT = 4 * STEFAN_BOLTZMANN / C_LIGHT

# FDT Quantities
const F_MAX = C_LIGHT^4 / (4 * G_NEWTON)  # ~3.03×10⁴³ N
const PLANCK_LENGTH = sqrt(PLANCK_HBAR * G_NEWTON / C_LIGHT^3)

# CMB Temperature Parameters (from triple identity paper)
# The universe's actual self-pressure temperature is 27 K
# We observe only 10% (2.7 K) because:
#   - Baryonic matter = 5% of total
#   - Photon = graviton, so we see it twice (5% × 2 = 10%)
# CMB is the interference pattern of ALL universal density waves, past and present
const T_UNIVERSE = 27.0  # Actual universe temperature (K) from self-pressure
const T_CMB = 2.7255     # Observed CMB temperature (K) = 10% of T_universe
const BARYON_FRACTION = 0.05  # 5% baryonic matter
const OBSERVABLE_FRACTION = 2 * BARYON_FRACTION  # 10% observable (photon=graviton)

const T_EULER = ℯ
const M_UNIVERSE_E = sqrt(3 * C_LIGHT^8 / (64π * RAD_CONSTANT * G_NEWTON^3 * ℯ^4))

# Density Parameter
schwarzschild_radius(M) = 2 * G_NEWTON * M / C_LIGHT^2
α_outside(r, R_s) = clamp(R_s / r, 0.0, 0.9999)
α_inside(r, R_s) = clamp(r / R_s, 0.0, 0.9999)

# Force and Pressure
fdt_force(α₁, α₂) = F_MAX * α₁ * α₂
universal_pressure(r) = C_LIGHT^4 / (16π * G_NEWTON * r^2)

# Central Invariant Check
function central_invariant_check()
    c² = C_LIGHT^2
    em = 1 / (EPSILON_0 * MU_0)
    return (c², em, abs(c² - em) / c²)
end

# Triple Identity Regime
function triple_identity_regime(α)
    α < 0.001 ? :gravitational : α < 0.3 ? :electromagnetic : α < 0.7 ? :weak : :strong
end

# e-Kelvin Cosmology
e_kelvin_convergence() = (T_CMB / T_EULER, (T_CMB / T_EULER - 1) * 100)

"""
    cmb_from_universe_temperature(T_univ)

Compute observed CMB temperature from actual universe temperature.
T_CMB = T_universe × 2 × f_baryon = T_universe × 0.10

The factor of 2 comes from the photon=graviton identity:
we see baryonic matter twice (once as photon, once as graviton).
"""
cmb_from_universe_temperature(T_univ) = T_univ * OBSERVABLE_FRACTION

# Regularization
regularize_with_fmax(F) = F_MAX * tanh(F / F_MAX)

# Growth Factor
function fdt_growth_factor(z, α₀; Ωm=0.3)
    a = 1 / (1 + z)
    D = a * (Ωm * a^(-3) + (1 - Ωm))^(-0.5)
    return D * (1 - α₀^2)
end

# Bias ↔ Alpha Mapping
bias_to_alpha(b; b_max=5.0) = tanh(abs(b) / b_max)
alpha_to_bias(α; b_max=5.0) = b_max * atanh(clamp(α, 0.001, 0.999))

# Power Spectrum Correction
function fdt_power_spectrum_correction(k, α; k_max=C_LIGHT^4/(4*G_NEWTON*PLANCK_HBAR))
    α_corr = 1 - α^2
    k_reg = exp(-(k / k_max)^2)
    return α_corr * (1 + k_reg * (1 - α_corr))
end

end # module
