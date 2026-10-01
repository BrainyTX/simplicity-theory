"""
sparc_tools.py – shared helpers for the SPARC analyses.

Loads the SPARC rotation curves, builds the baryonic mass model and provides
the gravity models that are compared against the data. All analysis scripts
in this folder (fit_sigma_crit.py, compare_models.py, rar_phase_plots.py,
plateau_test.py) import from here, so they share one data set, one set of
constants and one set of formulas.

Units
-----
radius           kpc
velocities       km/s
accelerations    (km/s)^2 / kpc      (multiply by ACC_TO_SI for m/s^2)
surface density  M_sun / pc^2
"""
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data" / "sparc" / "Rotmod_LTG"
RESULTS_DIR = REPO_ROOT / "results" / "sparc_analysis"

# --------------------------------------------------------------------------- #
# Physical constants
# --------------------------------------------------------------------------- #
G_KPC = 4.30091e-6                       # G in kpc (km/s)^2 / M_sun
ACC_TO_SI = 1e6 / 3.085677581e19         # (km/s)^2/kpc -> m/s^2
C_SI = 2.99792458e8                      # m/s
G_SI = 6.67430e-11                       # m^3 kg^-1 s^-2
M_SUN = 1.98847e30                       # kg
PC = 3.085677581e16                      # m

# --------------------------------------------------------------------------- #
# Thresholds of the Simplicity Theory
# --------------------------------------------------------------------------- #
# Upper threshold (Simplicity Codex, Sec. 45.6, quartet refinement):
#   sigma_dark = 137 - Omega_phi - 1/4 = 123.69,  Omega_phi = 2*pi/ln(phi)
# Above it the medium is "solid" (Newtonian gravity).
#
# Lower threshold: 137 - 8*Omega_phi = 32.54. Below it the medium is
# "gaseous" (time dissolves, dark-energy regime).
#
# The codex values are pure numbers, a surface density carries units. The
# lower threshold is therefore used as the dimensionless ratio
# LOW_TO_STAR_RATIO to the upper one; the absolute scale in M_sun/pc^2 comes
# from SIGMA_STAR (which agrees with c*H0/((2*pi)^2*G) for H0 ~ 70 km/s/Mpc,
# see sigma_star_from_h0).
PHI = (1 + 5 ** 0.5) / 2
OMEGA_PHI = 2 * np.pi / np.log(PHI)                       # 13.057
SIGMA_STAR = 137 - OMEGA_PHI - 0.25                       # 123.69 M_sun/pc^2
LOW_TO_STAR_RATIO = (137 - 8 * OMEGA_PHI) / SIGMA_STAR    # 0.2631
SIGMA_LOW = LOW_TO_STAR_RATIO * SIGMA_STAR                # 32.54 M_sun/pc^2

# Mass-to-light ratios at 3.6 micron. SPARC lists Vdisk and Vbul for M/L = 1;
# 0.5 (disk) and 0.7 (bulge) are the standard values (Lelli et al. 2016).
ML_DISK = 0.5
ML_BULGE = 0.7

COLUMNS = ["R", "Vobs", "eV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]


def sigma_star_from_h0(h0_km_s_mpc: float) -> float:
    """Sigma*_crit = c*H0 / ((2*pi)^2 * G) in M_sun/pc^2."""
    h0_si = h0_km_s_mpc * 1e3 / (1e6 * PC)
    sigma_kg_m2 = C_SI * h0_si / ((2 * np.pi) ** 2 * G_SI)
    return sigma_kg_m2 * PC ** 2 / M_SUN


def g_dagger(sigma: float) -> float:
    """Acceleration scale 2*pi*G*Sigma in (km/s)^2/kpc for Sigma in M_sun/pc^2."""
    return 2 * np.pi * G_KPC * sigma * 1e6


# --------------------------------------------------------------------------- #
# Data loading and mass model
# --------------------------------------------------------------------------- #
def load_galaxy(path: Path) -> pd.DataFrame:
    """Reads one SPARC *_rotmod.dat file."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=COLUMNS)
    df["galaxy"] = path.stem.replace("_rotmod", "")
    return df


def load_sparc(ml_disk: float = ML_DISK, ml_bulge: float = ML_BULGE,
               data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Loads all galaxies and adds the baryonic mass model (one row per radius)."""
    files = sorted(data_dir.glob("*_rotmod.dat"))
    if not files:
        raise FileNotFoundError(f"No *_rotmod.dat files found in {data_dir}")
    df = pd.concat([load_galaxy(p) for p in files], ignore_index=True)
    return add_mass_model(df, ml_disk, ml_bulge)


def add_mass_model(df: pd.DataFrame, ml_disk: float = ML_DISK,
                   ml_bulge: float = ML_BULGE) -> pd.DataFrame:
    """Adds baryonic and observed accelerations and the surface densities.

    SPARC convention: a negative velocity means that component pulls outward
    (e.g. a central hole in the gas disk), so contributions are added as V*|V|
    instead of V^2. Points without a positive baryonic or observed
    acceleration are dropped.
    """
    vbar2 = (df["Vgas"] * df["Vgas"].abs()
             + ml_disk * df["Vdisk"] * df["Vdisk"].abs()
             + ml_bulge * df["Vbul"] * df["Vbul"].abs())
    df = df.assign(Vbar2=vbar2)
    df = df[(df["R"] > 0) & (df["Vobs"] > 0) & (df["Vbar2"] > 0)].copy()

    df["g_bar"] = df["Vbar2"] / df["R"]
    df["g_obs"] = df["Vobs"] ** 2 / df["R"]
    # error of log10(g_obs) from the velocity error only
    df["e_log_gobs"] = 2 * df["eV"] / df["Vobs"] / np.log(10)
    # local stellar surface density (gas is not included in the SPARC files)
    df["sigma_stars"] = ml_disk * df["SBdisk"] + ml_bulge * df["SBbul"]
    # effective surface density: the Sigma that produces g_bar for a thin sheet,
    # includes gas, disk and bulge
    df["sigma_eff"] = df["g_bar"] / g_dagger(1.0)
    return df


# --------------------------------------------------------------------------- #
# Gravity models: each returns the predicted g_obs for every row of df
# --------------------------------------------------------------------------- #
def newton(df, **_):
    """Baryons only, no amplification."""
    return df["g_bar"].to_numpy()


def square_rule(df, sigma_star=SIGMA_STAR, **_):
    """Original Square-Rule of square_rule_basic.py / square_rule_smbh.py:
    g = g_bar * (1 + f),  f = 1 - (Sigma_stars / Sigma*)^2,  f in [0, 1]."""
    f = 1 - np.clip(df["sigma_stars"] / sigma_star, 0, 1) ** 2
    return (df["g_bar"] * (1 + f)).to_numpy()


def rar(df, sigma_star=SIGMA_STAR, **_):
    """Radial acceleration relation (McGaugh, Lelli & Schombert 2016) with
    g_dagger = 2*pi*G*Sigma*."""
    x = df["g_bar"] / g_dagger(sigma_star)
    return (df["g_bar"] / (1 - np.exp(-np.sqrt(x)))).to_numpy()


def mond_simple(df, sigma_star=SIGMA_STAR, **_):
    """MOND with the simple interpolation function mu(y) = y/(1+y), a0 = 2*pi*G*Sigma*."""
    a0 = g_dagger(sigma_star)
    gb = df["g_bar"]
    return (gb / 2 + np.sqrt(gb ** 2 / 4 + gb * a0)).to_numpy()


def three_phase_boost(x, x_low=LOW_TO_STAR_RATIO, gas_exponent=0.5):
    """Amplification g_obs/g_bar of the three-phase model, x = Sigma_eff / Sigma*.

    solid   x >= 1          : 1                        (Newton)
    liquid  x_low <= x < 1  : 2 - x^2                  (Square-Rule, 1 + f)
    gas     x < x_low       : (2 - x_low^2) * (x_low / x)^gas_exponent

    The gas law joins the liquid one continuously at x_low. gas_exponent = 0.5
    is what flat rotation curves require (g_obs ~ sqrt(g_bar)).
    """
    x = np.asarray(x, dtype=float)
    b_gas = (2 - x_low ** 2) * (x_low / x) ** gas_exponent
    return np.where(x >= 1, 1.0, np.where(x >= x_low, 2 - x ** 2, b_gas))


def three_phase(df, sigma_star=SIGMA_STAR, sigma_low=SIGMA_LOW, gas_exponent=0.5, **_):
    """Three-phase model (gas / liquid / solid), phase set by Sigma_eff."""
    x = df["g_bar"] / g_dagger(sigma_star)
    return (df["g_bar"] * three_phase_boost(x, sigma_low / sigma_star, gas_exponent)).to_numpy()


def three_phase_hybrid(df, sigma_star=SIGMA_STAR, sigma_low=SIGMA_LOW, gas_exponent=0.5, **_):
    """Hybrid: liquid/solid from the original Square-Rule on the local stellar
    Sigma, gas law on Sigma_eff; the stronger of the two applies."""
    x_low = sigma_low / sigma_star
    x = df["g_bar"] / g_dagger(sigma_star)
    b_square = 2 - np.clip(df["sigma_stars"] / sigma_star, 0, 1) ** 2
    b_gas = (2 - x_low ** 2) * (x_low / x) ** gas_exponent
    return (df["g_bar"] * np.maximum(b_square, b_gas)).to_numpy()


MODELS = {
    "Newton": newton,
    "Square-Rule (original)": square_rule,
    "RAR (McGaugh 2016)": rar,
    "MOND (simple mu)": mond_simple,
    "Three-phase": three_phase,
    "Three-phase hybrid": three_phase_hybrid,
}


# --------------------------------------------------------------------------- #
# Fit statistics
# --------------------------------------------------------------------------- #
# SPARC velocity errors do not contain distance, inclination and M/L
# uncertainties; this floor (in dex) is added in quadrature for global fits.
SYSTEMATIC_FLOOR_DEX = 0.1


def log_residuals(df, g_pred):
    """log10(g_obs) - log10(g_pred)."""
    return np.log10(df["g_obs"].to_numpy()) - np.log10(g_pred)


def chi2_global(df, g_pred, floor=SYSTEMATIC_FLOOR_DEX):
    """chi^2 of log residuals with velocity errors plus a systematic floor."""
    err = np.hypot(df["e_log_gobs"].to_numpy(), floor)
    return float(np.sum((log_residuals(df, g_pred) / err) ** 2))


def chi2_velocity(df, g_pred):
    """Classic per-galaxy chi^2 in velocity space with the SPARC errors."""
    v_pred = np.sqrt(g_pred * df["R"].to_numpy())
    err = np.where(df["eV"] > 0, df["eV"], 1.0)
    return float(np.sum(((df["Vobs"].to_numpy() - v_pred) / err) ** 2))


def binned_median(x, y, edges, min_count=10):
    """Median of y in bins of x; returns bin centres (geometric), medians, counts."""
    x, y = np.asarray(x), np.asarray(y)
    idx = np.digitize(x, edges) - 1
    centres, medians, counts = [], [], []
    for i in range(len(edges) - 1):
        sel = idx == i
        if sel.sum() >= min_count:
            centres.append(np.sqrt(edges[i] * edges[i + 1]))
            medians.append(np.median(y[sel]))
            counts.append(int(sel.sum()))
    return np.array(centres), np.array(medians), np.array(counts)
