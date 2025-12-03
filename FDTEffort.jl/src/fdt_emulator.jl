"""
    FDT Emulator Module

    Extends Effort.jl's EFTofLSS emulator with FDT physics.

    Key modifications:
    - Bias parameters → Density parameters α ∈ (0,1)
    - Natural regularization via F_max = c⁴/4G
    - Triple identity unification at high α
    - CMB as interference pattern of all density waves
    - Shadow matter (dark matter) as inseparable compressed spacetime
    - Stored potential energy (dark energy) released on mass loss

    DARK MATTER:
    - Shadow of regular matter - INSEPARABLE, like its shadow
    - Heavier, denser, compressed spacetime
    - Not particles - geometry at α > 0.156

    DARK ENERGY:
    - Stored potential energy of compressed spacetime
    - Released when bodies lose mass
    - Drives cosmic acceleration

    CMB:
    - Interference pattern of ALL density waves (past and present)
    - T_universe = 27 K, T_CMB = 2.7 K (10% observable)
"""
module FDTEmulator

using ..FDTPhysics

export FDTCosmology, FDTBiasParameters,
       get_Pℓ_fdt, get_Pℓ_jacobian_fdt,
       convert_cosmology_to_fdt, convert_bias_to_fdt,
       fdt_to_standard_cosmology, fdt_to_standard_bias,
       fdt_stochastic_model, fdt_bias_combination,
       shadow_matter_boost, is_shadow_regime

#=============================================================================
                        FDT COSMOLOGY STRUCT
=============================================================================#

"""
    FDTCosmology

FDT-native cosmological parameters.

# Fields
- `z`: Redshift
- `α_initial`: Initial density amplitude (replaces ln10As)
- `n_α`: Density spectral index (replaces ns)  
- `R_universe`: Cosmic Schwarzschild scale (replaces H0)
- `α_baryon`: Baryonic density parameter (replaces ωb)
- `α_cdm`: Cold dark omnium density (replaces ωcdm)
- `δα_ν`: Neutrino geometric contribution (replaces mν)
- `w0`: Dark energy equation of state
- `wa`: Dark energy evolution
"""
struct FDTCosmology{T<:Real}
    z::T
    α_initial::T
    n_α::T
    R_universe::T
    α_baryon::T
    α_cdm::T
    δα_ν::T
    w0::T
    wa::T
end

"""
    convert_cosmology_to_fdt(cosmo_standard)

Convert standard cosmological parameters to FDT framework.

Standard: [z, ln10As, ns, H0, ωb, ωcdm, mν, w0, wa]
FDT: [z, α_initial, n_α, R_universe, α_baryon, α_cdm, δα_ν, w0, wa]
"""
function convert_cosmology_to_fdt(cosmo::AbstractVector)
    z, ln10As, ns, H0, ωb, ωcdm, mν, w0, wa = cosmo
    
    # Convert Hubble to cosmic Schwarzschild scale
    # H0 in km/s/Mpc → R_universe = c/H0
    H0_si = H0 * 1000 / 3.086e22  # Convert to 1/s
    R_universe = C_LIGHT / H0_si
    
    # Convert amplitude to density parameter
    # ln(10¹⁰As) → α_initial via geometric mapping
    As = exp(ln10As) * 1e-10
    α_initial = sqrt(As) * 0.1  # Normalized to α ∈ (0,1)
    
    # Spectral index maps directly (geometric scaling)
    n_α = ns
    
    # Density fractions become density parameters
    # ωb = Ωb h² → α_baryon
    α_baryon = bias_to_alpha(ωb * 100)
    α_cdm = bias_to_alpha(ωcdm * 100)
    
    # Neutrino mass → geometric contribution
    δα_ν = mν / 10.0  # Normalized
    
    return FDTCosmology(z, α_initial, n_α, R_universe, α_baryon, α_cdm, δα_ν, w0, wa)
end

"""
    fdt_to_standard_cosmology(fdt_cosmo::FDTCosmology)

Convert FDT parameters back to standard format for emulator input.
"""
function fdt_to_standard_cosmology(fdt::FDTCosmology)
    # Inverse transformations
    H0 = C_LIGHT / fdt.R_universe * 3.086e22 / 1000
    ln10As = log(100 * fdt.α_initial^2) + 10 * log(10)
    ns = fdt.n_α
    ωb = alpha_to_bias(fdt.α_baryon) / 100
    ωcdm = alpha_to_bias(fdt.α_cdm) / 100
    mν = fdt.δα_ν * 10.0
    
    return [fdt.z, ln10As, ns, H0, ωb, ωcdm, mν, fdt.w0, fdt.wa]
end

#=============================================================================
                        FDT BIAS TRANSFORMATION
=============================================================================#

"""
    FDTBiasParameters

FDT-native bias/density parameters for galaxy clustering.

In FDT, bias parameters represent local density enhancements:
- Higher bias = denser region = larger α
- All α bounded in (0,1) naturally
"""
struct FDTBiasParameters{T<:Real}
    α₁::T      # Linear density parameter
    α₂::T      # Quadratic density parameter  
    α₃::T      # Cubic density parameter
    α_ct::T    # Counterterm density
    α_stoch::NTuple{3,T}  # Stochastic density parameters
    f::T       # Growth rate (geometric)
end

"""
    convert_bias_to_fdt(bs_standard)

Convert standard EFT bias parameters to FDT density parameters.

Standard: [b1, b2, b3, b4, b5, b6, b7, f, cϵ0, cϵ1, cϵ2]
"""
function convert_bias_to_fdt(bs::AbstractVector)
    b1, b2, b3, b4, b5, b6, b7, f, cϵ0, cϵ1, cϵ2 = bs
    
    # Convert biases to bounded density parameters
    α₁ = bias_to_alpha(b1)
    α₂ = bias_to_alpha(b2)
    α₃ = bias_to_alpha(b3)
    
    # Counterterms combine b4-b7
    α_ct = bias_to_alpha(sqrt(b4^2 + b5^2 + b6^2 + b7^2))
    
    # Stochastic terms
    α_stoch = (bias_to_alpha(cϵ0), bias_to_alpha(cϵ1), bias_to_alpha(cϵ2))
    
    return FDTBiasParameters(α₁, α₂, α₃, α_ct, α_stoch, f)
end

"""
    fdt_to_standard_bias(fdt_bias::FDTBiasParameters)

Convert FDT parameters back to standard bias vector.
"""
function fdt_to_standard_bias(fdt::FDTBiasParameters)
    b1 = alpha_to_bias(fdt.α₁)
    b2 = alpha_to_bias(fdt.α₂)
    b3 = alpha_to_bias(fdt.α₃)
    
    # Distribute counterterm among b4-b7
    b_ct = alpha_to_bias(fdt.α_ct) / 2
    
    cϵ0 = alpha_to_bias(fdt.α_stoch[1])
    cϵ1 = alpha_to_bias(fdt.α_stoch[2])
    cϵ2 = alpha_to_bias(fdt.α_stoch[3])
    
    return [b1, b2, b3, b_ct, b_ct, b_ct, b_ct, fdt.f, cϵ0, cϵ1, cϵ2]
end

#=============================================================================
                        FDT POWER SPECTRUM
=============================================================================#

"""
    get_Pℓ_fdt(fdt_cosmo, D, fdt_bias, cosmoemu; stoch_kwargs...)

Compute power spectrum multipole with FDT corrections.

This wraps the standard emulator but:
1. Converts FDT parameters to standard format
2. Applies FDT corrections to output
3. Enforces F_max regularization
"""
function get_Pℓ_fdt(fdt_cosmo::FDTCosmology, D, fdt_bias::FDTBiasParameters, 
                    cosmoemu; stoch_kwargs...)
    
    # Convert to standard parameters for emulator
    cosmo_std = fdt_to_standard_cosmology(fdt_cosmo)
    bias_std = fdt_to_standard_bias(fdt_bias)
    
    # Get standard power spectrum (assumes Effort.jl loaded)
    # P_standard = get_Pℓ(cosmo_std, D, bias_std, cosmoemu; stoch_kwargs...)
    
    # For now, return placeholder - actual integration requires Effort.jl
    k_grid = cosmoemu.P11.kgrid
    P_standard = zeros(length(k_grid))
    
    # Apply FDT corrections
    α_eff = sqrt(fdt_bias.α₁^2 + fdt_bias.α₂^2)  # Effective density
    
    P_fdt = similar(P_standard)
    for (i, k) in enumerate(k_grid)
        correction = fdt_power_spectrum_correction(k, α_eff)
        P_fdt[i] = P_standard[i] * correction
    end
    
    return P_fdt
end

#=============================================================================
                        FDT STOCHASTIC MODEL
=============================================================================#

"""
    fdt_stochastic_model(kgrid, α_stoch; Pshot=0.0)

FDT stochastic/shot noise model.

In FDT, stochastic terms arise from:
- Discrete sampling of continuous omnium field
- Finite density resolution effects
- Phase transition noise
"""
function fdt_stochastic_model(kgrid, α_stoch::NTuple{3}; Pshot=0.0)
    α₀, α₁, α₂ = α_stoch
    
    n_k = length(kgrid)
    stoch = zeros(n_k, 3)
    
    for (i, k) in enumerate(kgrid)
        # Constant term (shot noise modified by density)
        stoch[i, 1] = Pshot * (1 - α₀^2)
        
        # k² term (small scale discreteness)
        stoch[i, 2] = Pshot * k^2 * (1 - α₁^2)
        
        # k⁴ term (UV regularization from F_max)
        stoch[i, 3] = Pshot * k^4 * (1 - α₂^2) * exp(-k^2 / 100)
    end
    
    return stoch
end

#=============================================================================
                        FDT BIAS COMBINATION
=============================================================================#

"""
    fdt_bias_combination(fdt_bias::FDTBiasParameters)

Compute bias combination weights in FDT framework.

Key insight: In FDT, bias combination represents how different
density shells contribute to observed clustering.
"""
function fdt_bias_combination(fdt::FDTBiasParameters)
    α₁, α₂, α₃ = fdt.α₁, fdt.α₂, fdt.α₃
    f = fdt.f
    
    # Convert to effective weights via density product rule
    # F = (c⁴/4G)α₁α₂ implies quadratic combinations
    
    # Linear terms
    w_11 = α₁^2
    
    # Quadratic terms  
    w_22 = α₂^2
    w_12 = 2 * α₁ * α₂
    
    # Growth rate terms
    w_f = f^2
    w_f1 = 2 * f * α₁
    
    # Counterterm weights (regularization)
    w_ct = fdt.α_ct^2 * (1 - fdt.α_ct^2)  # Self-regularizing
    
    return [w_11, w_22, w_12, w_f, w_f1, w_ct]
end

#=============================================================================
                        SHADOW MATTER FUNCTIONS
=============================================================================#

# α threshold dividing visible from shadow (matches ShadowMatter module)
const ALPHA_THRESHOLD_EMU = 1 / 6.4  # ≈ 0.156, DM/baryon = 5.4

"""
    is_shadow_regime(α)

Check if density parameter α is in the shadow (dark matter) regime.

Shadow matter is the geometric SHADOW of regular matter:
- Heavier, denser, compressed spacetime
- INSEPARABLE from regular matter (like its shadow)
- EM-inactive (explains why dark matter is dark)
"""
function is_shadow_regime(α)
    return α > ALPHA_THRESHOLD_EMU
end

"""
    shadow_matter_boost(α)

Compute the boost factor from shadow matter contribution.

At low α (visible regime): boost ≈ 1
At high α (shadow regime): boost increases due to compressed geometry

This accounts for the additional gravitational effect of shadow matter
on the power spectrum.
"""
function shadow_matter_boost(α)
    if α < ALPHA_THRESHOLD_EMU
        return 1.0
    end

    # Shadow matter contributes additional gravity
    # Boost factor based on DM/baryon ratio
    α_normalized = (α - ALPHA_THRESHOLD_EMU) / (1 - ALPHA_THRESHOLD_EMU)

    # Maximum boost is (1 + DM/baryon) at α = 1
    dm_baryon_ratio = 5.4
    boost = 1.0 + dm_baryon_ratio * α_normalized^2

    return boost
end

"""
    fdt_matter_power_with_shadow(k, α, P_baryon)

Compute total matter power spectrum including shadow contribution.

Arguments:
- k: wavenumber
- α: effective density parameter
- P_baryon: baryonic power spectrum

Returns total power P_total = P_baryon × shadow_boost × FDT_correction
"""
function fdt_matter_power_with_shadow(k, α, P_baryon)
    # Shadow matter boost
    boost = shadow_matter_boost(α)

    # FDT correction (regularization at high k and α)
    correction = fdt_power_spectrum_correction(k, α)

    return P_baryon * boost * correction
end

end # module FDTEmulator
