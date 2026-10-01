#!/usr/bin/env python3
"""
sigma_crit.py  –  Critical surface density Σ*_crit from fundamental constants.

Usage (CLI):
    python sigma_crit.py            # mit Default-H0 = 70 km/s/Mpc
    python sigma_crit.py 73.0       # eigenes H0 in km/s/Mpc

Output:
    a0        = 1.08e-10  m/s²
    Sigma*    = 0.258     kg/m²
    Sigma*    = 123.6     Msol/pc²
"""

import sys
import math

# --- Konstanten (CODATA 2022) -----------------------------------------------
c      = 2.997_924_58e8       # m/s
G      = 6.674_30e-11         # m³ kg⁻¹ s⁻²
Msol   = 1.988_47e30          # kg
pc     = 3.085_677_581e16     # m

def sigma_crit(H0_km_s_Mpc: float = 70.0) -> None:
    """Prints a0 and Σ*_crit for a given H0 (km/s/Mpc)."""
    H0 = H0_km_s_Mpc * 1e3 / (1e6 * pc)   # -> 1/s
    a0 = c * H0 / (2*math.pi)             # m/s²
    sigma = a0 / (2*math.pi*G)            # kg/m²
    sigma_msol_pc2 = sigma * pc**2 / Msol # Msol/pc²

    print(f"H0        = {H0_km_s_Mpc:.2f}  km/s/Mpc")
    print(f"a0        = {a0:.3e}  m/s²")
    print(f"Sigma*    = {sigma:.3f}  kg/m²")
    print(f"Sigma*    = {sigma_msol_pc2:.1f}  Msol/pc²")

if __name__ == "__main__":
    h0 = float(sys.argv[1]) if len(sys.argv) > 1 else 70.0
    sigma_crit(h0)
