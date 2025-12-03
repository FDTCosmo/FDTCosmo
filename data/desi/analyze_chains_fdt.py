#!/usr/bin/env python3
"""
FDT Analysis of DESI MCMC Chains

Analyzes DESI DR1 cosmological parameter posteriors through the lens of
Fundamental Density Theory (FDT), extracting density parameters and
checking consistency with FDT predictions.

Key FDT Transformations:
- H0 -> R_universe = c/H0 (cosmic Schwarzschild scale)
- Omega_m h^2 -> alpha parameters (density parameters in [0,1])
- w, wa -> dark energy as stored potential energy

FDT Predictions to Test:
1. Central invariant: E/m = c^2 = 1/(eps0*mu0)
2. Dark matter as shadow of baryonic matter (inseparable)
3. CMB temperature: T_universe = 27K, T_CMB = 10% = 2.7K
4. Maximum force regularization: F_max = c^4/4G

Author: FDT Universe Simulator
"""

import numpy as np
import pandas as pd
from pathlib import Path
import warnings
from dataclasses import dataclass
from typing import Tuple, Dict, List, Optional
import json

# Physical constants
C_LIGHT = 299792458.0  # m/s
G_NEWTON = 6.67430e-11  # m^3 kg^-1 s^-2
EPSILON_0 = 8.8541878128e-12  # F/m
MU_0 = 1.25663706212e-6  # H/m
PLANCK_HBAR = 1.054571817e-34  # J s
K_BOLTZMANN = 1.380649e-23  # J/K

# FDT quantities
F_MAX = C_LIGHT**4 / (4 * G_NEWTON)  # Maximum force ~3.03e43 N
T_UNIVERSE_FDT = 27.0  # Actual universe temperature (K)
T_CMB_FDT = 2.7255  # Observed CMB (10% of T_universe)
BARYON_FRACTION = 0.05  # 5% baryonic matter
DM_BARYON_RATIO = 5.4  # Dark matter to baryon ratio
ALPHA_THRESHOLD = 1 / 6.4  # ~0.156, threshold between visible/shadow regimes


@dataclass
class FDTCosmology:
    """FDT-native cosmological parameters."""
    z: float = 0.0
    alpha_initial: float = 0.0  # Density amplitude
    n_alpha: float = 0.0  # Density spectral index
    R_universe: float = 0.0  # Cosmic Schwarzschild scale (m)
    alpha_baryon: float = 0.0  # Baryonic density parameter
    alpha_cdm: float = 0.0  # CDM density parameter
    delta_alpha_nu: float = 0.0  # Neutrino contribution
    w0: float = -1.0  # Dark energy equation of state
    wa: float = 0.0  # Dark energy evolution


def bias_to_alpha(b: float, b_max: float = 5.0) -> float:
    """Convert bias to bounded density parameter alpha in [0,1).

    In FDT, alpha = R_s/r represents the density parameter.
    alpha = 0: empty space (r >> R_s)
    alpha -> 1: approaching event horizon (r -> R_s)
    alpha must NEVER reach or exceed 1.
    """
    return np.clip(np.tanh(np.abs(b) / b_max), 0.0, 0.9999)


def alpha_to_bias(alpha: float, b_max: float = 5.0) -> float:
    """Convert density parameter alpha to bias."""
    alpha = np.clip(alpha, 0.001, 0.999)
    return b_max * np.arctanh(alpha)


def convert_cosmology_to_fdt(
    H0: float,
    omega_m: float = None,
    omega_b: float = None,
    omega_cdm: float = None,
    sigma8: float = None,
    ns: float = None,
    mnu: float = 0.0,
    w0: float = -1.0,
    wa: float = 0.0,
    z: float = 0.0
) -> FDTCosmology:
    """
    Convert standard cosmological parameters to FDT framework.

    Parameters from DESI chains:
    - H0: Hubble constant (km/s/Mpc)
    - omega_m: Total matter density Omega_m
    - omega_b: Baryon density Omega_b h^2 (if available)
    - omega_cdm: CDM density Omega_cdm h^2 (if available)
    - sigma8: Matter fluctuation amplitude
    - ns: Spectral index
    - mnu: Neutrino mass sum (eV)
    - w0, wa: Dark energy parameters
    """
    # Convert H0 to cosmic Schwarzschild scale
    H0_si = H0 * 1000 / 3.086e22  # Convert km/s/Mpc to 1/s
    R_universe = C_LIGHT / H0_si

    # Convert amplitude to density parameter
    if sigma8 is not None:
        alpha_initial = sigma8 / 10.0  # Normalized to alpha in (0,1)
    else:
        alpha_initial = 0.08  # Default ~sigma8=0.8

    # Spectral index maps directly
    n_alpha = ns if ns is not None else 0.96

    # Density fractions become density parameters
    # In FDT, alpha = R_s/r is the density parameter
    # For cosmic average: alpha relates to Omega through the cosmic density
    # alpha_i = Omega_i * (critical density effect)
    # We use a scaling that maps Omega_i h^2 to alpha in (0,1)
    #
    # Physical interpretation:
    # - alpha_baryon ~ 0.05 (baryons are 5% of critical density)
    # - alpha_cdm ~ 0.27 (CDM is ~27% of critical density)
    # - alpha_total < 1 always (never at event horizon on cosmic scales)

    h = H0 / 100.0

    # Direct mapping: alpha = Omega_component (bounded to 0.9999)
    # This preserves the physical meaning of density fractions
    if omega_b is not None:
        # omega_b is Omega_b * h^2, typically ~0.022
        # alpha_baryon should reflect baryon fraction ~0.05
        alpha_baryon = np.clip(omega_b / h**2, 0.001, 0.9999) if h > 0 else 0.05
    elif omega_m is not None:
        alpha_baryon = np.clip(omega_m * BARYON_FRACTION / (1 + DM_BARYON_RATIO * BARYON_FRACTION), 0.001, 0.9999)
    else:
        alpha_baryon = 0.05

    if omega_cdm is not None:
        # omega_cdm is Omega_cdm * h^2, typically ~0.12
        # alpha_cdm should reflect CDM fraction ~0.27
        alpha_cdm = np.clip(omega_cdm / h**2, 0.001, 0.9999) if h > 0 else 0.27
    elif omega_m is not None:
        alpha_cdm = np.clip(omega_m * DM_BARYON_RATIO * BARYON_FRACTION / (1 + DM_BARYON_RATIO * BARYON_FRACTION), 0.001, 0.9999)
    else:
        alpha_cdm = 0.27

    # Neutrino contribution
    delta_alpha_nu = mnu / 10.0 if mnu is not None else 0.0

    return FDTCosmology(
        z=z,
        alpha_initial=alpha_initial,
        n_alpha=n_alpha,
        R_universe=R_universe,
        alpha_baryon=alpha_baryon,
        alpha_cdm=alpha_cdm,
        delta_alpha_nu=delta_alpha_nu,
        w0=w0,
        wa=wa
    )


def central_invariant_check() -> Tuple[float, float, float]:
    """Verify the central FDT invariant: c^2 = 1/(eps0*mu0)."""
    c_squared = C_LIGHT**2
    em_product = 1 / (EPSILON_0 * MU_0)
    fractional_diff = abs(c_squared - em_product) / c_squared
    return c_squared, em_product, fractional_diff


def triple_identity_regime(alpha: float) -> str:
    """Determine which regime the density parameter is in."""
    if alpha < 0.001:
        return "gravitational"
    elif alpha < 0.3:
        return "electromagnetic"
    elif alpha < 0.7:
        return "weak"
    else:
        return "strong"


def shadow_matter_boost(alpha: float) -> float:
    """Compute boost factor from shadow matter contribution."""
    if alpha < ALPHA_THRESHOLD:
        return 1.0
    alpha_normalized = (alpha - ALPHA_THRESHOLD) / (1 - ALPHA_THRESHOLD)
    return 1.0 + DM_BARYON_RATIO * alpha_normalized**2


def fdt_power_spectrum_correction(k: float, alpha: float, k_max: float = 1e50) -> float:
    """FDT correction to power spectrum."""
    alpha_corr = 1 - alpha**2
    k_reg = np.exp(-(k / k_max)**2)
    return alpha_corr * (1 + k_reg * (1 - alpha_corr))


def load_chain(chain_path: Path) -> pd.DataFrame:
    """Load a single MCMC chain file."""
    try:
        # Read header to get column names
        with open(chain_path, 'r') as f:
            header = f.readline().strip()

        # Parse column names (space-separated, starting with #)
        if header.startswith('#'):
            columns = header[1:].split()
        else:
            columns = None

        # Load data
        df = pd.read_csv(chain_path, sep=r'\s+', comment='#', header=None, names=columns)
        return df
    except Exception as e:
        print(f"Error loading {chain_path}: {e}")
        return None


def extract_cosmology_params(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Extract relevant cosmological parameters from chain DataFrame."""
    params = {}

    # Map of common parameter names to look for
    param_map = {
        'H0': ['H0', 'h0', 'hubble'],
        'omegam': ['omegam', 'omega_m', 'Om', 'Omega_m'],
        'omegamh2': ['omegamh2', 'omega_m_h2', 'Omh2'],
        'ombh2': ['ombh2', 'omega_b_h2', 'Obh2'],
        'omch2': ['omch2', 'omega_cdm_h2', 'Ocdmh2'],
        'sigma8': ['sigma8', 's8', 'S8'],
        'ns': ['ns', 'n_s', 'spectral_index'],
        'logA': ['logA', 'ln10As', 'As', 'ln10^{10}As'],
        'w': ['w', 'w0', 'w_0'],
        'wa': ['wa', 'w_a'],
        'mnu': ['mnu', 'm_nu', 'sum_mnu', 'omnuh2'],
        'weight': ['weight'],
        'minuslogpost': ['minuslogpost', '-logpost', 'chi2']
    }

    for param_name, aliases in param_map.items():
        for alias in aliases:
            if alias in df.columns:
                params[param_name] = df[alias].values
                break

    return params


def analyze_chain_fdt(chain_dir: Path) -> Dict:
    """
    Analyze all chains in a directory through FDT lens.

    Returns dictionary with:
    - Parameter statistics (mean, std, median, 68% CI)
    - FDT-converted parameters
    - Regime analysis
    - Consistency checks
    """
    results = {
        'chain_dir': str(chain_dir),
        'n_samples': 0,
        'params': {},
        'fdt_params': {},
        'regime_analysis': {},
        'consistency': {}
    }

    # Load all chains
    all_data = []
    for chain_file in sorted(chain_dir.glob('chain.*.txt')):
        if 'post' in chain_file.name:
            continue  # Skip post-processed chains
        df = load_chain(chain_file)
        if df is not None:
            all_data.append(df)

    if not all_data:
        print(f"No chains found in {chain_dir}")
        return results

    # Combine chains
    combined = pd.concat(all_data, ignore_index=True)
    results['n_samples'] = len(combined)

    # Extract parameters
    params = extract_cosmology_params(combined)

    # Get weights for weighted statistics
    weights = params.get('weight', np.ones(len(combined)))

    # Compute statistics for each parameter
    for name, values in params.items():
        if name == 'weight':
            continue

        # Weighted statistics
        w_sum = np.sum(weights)
        mean = np.sum(weights * values) / w_sum
        var = np.sum(weights * (values - mean)**2) / w_sum
        std = np.sqrt(var)

        # Percentiles (unweighted for simplicity)
        p16, p50, p84 = np.percentile(values, [16, 50, 84])

        results['params'][name] = {
            'mean': float(mean),
            'std': float(std),
            'median': float(p50),
            'ci_low': float(p16),
            'ci_high': float(p84)
        }

    # Convert to FDT parameters
    if 'H0' in params:
        H0_mean = results['params']['H0']['mean']
        omega_m = results['params'].get('omegam', {}).get('mean')
        sigma8 = results['params'].get('sigma8', {}).get('mean')
        ns = results['params'].get('ns', {}).get('mean')
        ombh2 = results['params'].get('ombh2', {}).get('mean')
        omch2 = results['params'].get('omch2', {}).get('mean')
        w0 = results['params'].get('w', {}).get('mean', -1.0)
        wa = results['params'].get('wa', {}).get('mean', 0.0)

        fdt = convert_cosmology_to_fdt(
            H0=H0_mean,
            omega_m=omega_m,
            omega_b=ombh2,
            omega_cdm=omch2,
            sigma8=sigma8,
            ns=ns,
            w0=w0 if w0 else -1.0,
            wa=wa if wa else 0.0
        )

        results['fdt_params'] = {
            'R_universe_m': fdt.R_universe,
            'R_universe_Gly': fdt.R_universe / (9.461e15 * 1e9),  # Convert to Gly
            'alpha_initial': fdt.alpha_initial,
            'n_alpha': fdt.n_alpha,
            'alpha_baryon': fdt.alpha_baryon,
            'alpha_cdm': fdt.alpha_cdm,
            'delta_alpha_nu': fdt.delta_alpha_nu,
            'w0': fdt.w0,
            'wa': fdt.wa
        }

        # Regime analysis
        # In FDT, alpha_total represents the combined density effect
        # Using the FDT combination rule: alpha_combined = alpha1 + alpha2 - alpha1*alpha2
        # This ensures alpha_total stays in (0,1) when alpha1, alpha2 are in (0,1)
        # Physical interpretation: overlapping density regions don't simply add
        alpha_total = fdt.alpha_baryon + fdt.alpha_cdm - fdt.alpha_baryon * fdt.alpha_cdm
        alpha_total = np.clip(alpha_total, 0.0, 0.9999)  # Ensure bounded

        results['regime_analysis'] = {
            'alpha_total': alpha_total,
            'alpha_baryon': fdt.alpha_baryon,
            'alpha_cdm': fdt.alpha_cdm,
            'regime': triple_identity_regime(alpha_total),
            'shadow_boost': shadow_matter_boost(alpha_total),
            'in_shadow_regime': alpha_total > ALPHA_THRESHOLD,
            'dm_baryon_ratio_measured': fdt.alpha_cdm / fdt.alpha_baryon if fdt.alpha_baryon > 0 else None
        }

    # Consistency checks
    c2, em, diff = central_invariant_check()
    results['consistency'] = {
        'central_invariant': {
            'c_squared': c2,
            'eps0_mu0_inverse': em,
            'fractional_diff': diff,
            'verified': diff < 1e-10
        },
        'cmb_temperature': {
            'T_CMB_observed': T_CMB_FDT,
            'T_universe_predicted': T_UNIVERSE_FDT,
            'observable_fraction': T_CMB_FDT / T_UNIVERSE_FDT,
            'matches_baryon_fraction': abs(T_CMB_FDT / T_UNIVERSE_FDT - 2 * BARYON_FRACTION) < 0.01
        }
    }

    return results


def analyze_all_chains(base_dir: Path) -> List[Dict]:
    """Analyze all MCMC chain directories."""
    results = []

    mcmc_dir = base_dir / 'mcmc_chains' / 'cobaya'
    if not mcmc_dir.exists():
        print(f"MCMC chains directory not found: {mcmc_dir}")
        return results

    # Find all chain directories
    for model_dir in sorted(mcmc_dir.iterdir()):
        if not model_dir.is_dir():
            continue
        for analysis_dir in sorted(model_dir.iterdir()):
            if not analysis_dir.is_dir():
                continue
            if list(analysis_dir.glob('chain.*.txt')):
                print(f"\nAnalyzing: {model_dir.name}/{analysis_dir.name}")
                result = analyze_chain_fdt(analysis_dir)
                result['model'] = model_dir.name
                result['analysis'] = analysis_dir.name
                results.append(result)

    return results


def print_fdt_summary(results: List[Dict]):
    """Print summary of FDT analysis results."""
    print("\n" + "=" * 80)
    print("FUNDAMENTAL DENSITY THEORY (FDT) ANALYSIS OF DESI DR1 MCMC CHAINS")
    print("=" * 80)

    # Central invariant check
    c2, em, diff = central_invariant_check()
    print("\n[Central Invariant Verification]")
    print(f"  E/m = c^2 = 1/(eps0*mu0) = (d/t)^2")
    print(f"  c^2 = {c2:.10e} m^2/s^2")
    print(f"  1/(eps0*mu0) = {em:.10e} m^2/s^2")
    print(f"  Fractional difference: {diff:.2e} {'VERIFIED' if diff < 1e-10 else 'FAILED'}")

    print("\n[CMB Temperature Prediction]")
    print(f"  T_universe (FDT) = {T_UNIVERSE_FDT} K")
    print(f"  T_CMB (observed) = {T_CMB_FDT} K")
    print(f"  Observable fraction = {100*T_CMB_FDT/T_UNIVERSE_FDT:.1f}%")
    print(f"  (2 x baryon fraction = {100*2*BARYON_FRACTION:.1f}%)")

    print("\n[Maximum Force]")
    print(f"  F_max = c^4/(4G) = {F_MAX:.3e} N")

    for result in results:
        if result['n_samples'] == 0:
            continue

        print("\n" + "-" * 80)
        print(f"Model: {result.get('model', 'unknown')}")
        print(f"Analysis: {result.get('analysis', 'unknown')}")
        print(f"Samples: {result['n_samples']:,}")

        if 'params' in result and result['params']:
            print("\n[Standard Cosmological Parameters]")
            for name, stats in result['params'].items():
                if isinstance(stats, dict):
                    print(f"  {name}: {stats['mean']:.4f} +/- {stats['std']:.4f} "
                          f"(68%: [{stats['ci_low']:.4f}, {stats['ci_high']:.4f}])")

        if 'fdt_params' in result and result['fdt_params']:
            print("\n[FDT Parameters]")
            fdt = result['fdt_params']
            print(f"  R_universe = {fdt['R_universe_Gly']:.2f} Gly "
                  f"({fdt['R_universe_m']:.3e} m)")
            print(f"  alpha_initial = {fdt['alpha_initial']:.4f}")
            print(f"  n_alpha = {fdt['n_alpha']:.4f}")
            print(f"  alpha_baryon = {fdt['alpha_baryon']:.4f}")
            print(f"  alpha_cdm = {fdt['alpha_cdm']:.4f}")
            print(f"  w0 = {fdt['w0']:.3f}")
            print(f"  wa = {fdt['wa']:.3f}")

        if 'regime_analysis' in result and result['regime_analysis']:
            print("\n[Regime Analysis]")
            regime = result['regime_analysis']
            print(f"  alpha_total = {regime['alpha_total']:.4f}")
            print(f"  Regime: {regime['regime']}")
            print(f"  Shadow boost factor: {regime['shadow_boost']:.3f}")
            print(f"  In shadow regime: {regime['in_shadow_regime']}")
            if regime.get('dm_baryon_ratio_measured'):
                print(f"  DM/baryon ratio: {regime['dm_baryon_ratio_measured']:.2f} "
                      f"(expected: {DM_BARYON_RATIO})")


def main():
    """Main analysis function."""
    base_dir = Path(__file__).parent

    print("FDT Analysis of DESI DR1 MCMC Chains")
    print("=" * 50)

    # Run analysis
    results = analyze_all_chains(base_dir)

    if not results:
        print("\nNo chains analyzed. Make sure chains are downloaded to mcmc_chains/cobaya/")
        return

    # Print summary
    print_fdt_summary(results)

    # Save results to JSON
    output_file = base_dir / 'fdt_analysis_results.json'

    # Convert numpy types for JSON serialization
    def convert_types(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(i) for i in obj]
        return obj

    with open(output_file, 'w') as f:
        json.dump(convert_types(results), f, indent=2)

    print(f"\nResults saved to: {output_file}")

    return results


if __name__ == '__main__':
    main()
