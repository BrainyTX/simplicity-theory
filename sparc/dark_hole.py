import numpy as np
import pandas as pd

# Deine Konstante aus der Simplicity Theory
SIGMA_STAR = 124.0  # M_sun / pc^2

def r_dh(m_msun):
    """
    Berechnet den Dark-Hole-Radius (r_DH) in Parsec
    für eine gegebene Masse in Sonnenmassen (M_sun).
    Formel: r_DH = sqrt(M / (4 * pi * SIGMA_STAR))
    """
    return np.sqrt(m_msun / (4 * np.pi * SIGMA_STAR))

# Katalog von bekannten Schwarzen Löchern
bh_catalog = pd.DataFrame({
    "Name": ["Sgr A*", "M87*", "NGC 1277", "Cygnus X-1"],
    "M_Msun": [4.3e6, 6.5e9, 1.7e10, 14.8]
})

# Berechnung des r_DH für jeden Eintrag
bh_catalog["r_DH_pc"] = r_dh(bh_catalog["M_Msun"])

# Ausgabe der formatierten Tabelle
print("--- Berechnung des Dark-Hole-Radius (r_DH) ---")
print(bh_catalog.round({"r_DH_pc": 3}))
