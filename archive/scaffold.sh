#!/usr/bin/env bash
# ============================================================
# scaffold.sh – legt ein komplettes Projekt­gerüst für
#               „simplicity-theory“ an (Ordner, Skripte, Git)
# ============================================================
set -e   # Sofort abbrechen, falls ein Befehl fehlschlägt

# -------- Metadaten (bei Bedarf anpassen) -------------------
REPO_NAME="simplicity-theory"
AUTHOR="Pieter Goldau"
EMAIL="simplicity-theory@gmx.de"
VERSION="0.1.0"

# -------- Ordnerstruktur ------------------------------------
mkdir -p "${REPO_NAME}"/{src,docs,notebooks,data,tests}
cd "${REPO_NAME}"

# -------- README.md -----------------------------------------
cat > README.md <<EOF
# ${REPO_NAME}

Mathematische Tools zur **Existential Coherence Theory**  
Version ${VERSION} – © ${AUTHOR}

## Schnelle Nutzung

\`\`\`bash
# Umgebung anlegen (optional, aber empfohlen)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Vorhersage bei H₀ = 70 km/s/Mpc prüfen
python src/sigma_crit.py --H0 70
\`\`\`
EOF

# -------- requirements.txt ----------------------------------
cat > requirements.txt <<'EOF'
numpy
astropy
matplotlib
EOF

# -------- src/sigma_crit.py ---------------------------------
cat > src/sigma_crit.py <<'EOF'
#!/usr/bin/env python
"""
sigma_crit.py
Berechnet die Grenz-Beschleunigung a₀ und die kritische
Oberflächendichte Σ*_crit nach Goldau (2025).

Formeln:
    a0       = c * H0 / (2π)
    sigma_kg = a0 / (2π G)
    sigma_ms = sigma_kg * (pc_in_m)**2 / M_SUN
"""

import argparse, math

# --- Konstanten (CODATA 2022) --------------------------------
C        = 2.99792458e8            # m  / s
G        = 6.67430e-11             # m³ / (kg s²)
M_SUN    = 1.98847e30              # kg
PC_IN_M  = 3.085677581e16          # m  (1 pc)
TWOPI    = 2.0 * math.pi

def sigma_crit(H0_km_s_Mpc: float):
    """Gibt a0 [m/s²],  Σ*_crit [kg/m²],  Σ*_crit [M_sun/pc²] zurück."""
    H0_SI = H0_km_s_Mpc * 1e3 / 3.085677581e22   # km/s/Mpc  →  1/s
    a0    = C * H0_SI / TWOPI
    sigma_kg = a0 / (TWOPI * G)
    sigma_msun = sigma_kg * (PC_IN_M**2) / M_SUN
    return a0, sigma_kg, sigma_msun

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Σ*_crit-Rechner")
    p.add_argument("--H0", type=float, default=70.0,
                   help="Hubble-Konstante in km/s/Mpc (Standard: 70)")
    args = p.parse_args()

    a0, sig_kg, sig_msun = sigma_crit(args.H0)
    print(f"H0        : {args.H0:6.2f}  km/s/Mpc")
    print(f"a0        : {a0: .3e}  m/s²")
    print(f"Σ*_crit   : {sig_kg: .3e}  kg/m²")
    print(f"Σ*_crit   : {sig_msun: .2f}   M_sun/pc²")
EOF
chmod +x src/sigma_crit.py

# -------- Quick-Notebook (Mini-Demo) -------------------------
cat > notebooks/quick_check.ipynb <<'EOF'
{
 "cells": [
  {
   "cell_type": "code",
   "metadata": {},
   "source": ["!python ../src/sigma_crit.py --H0 72"]
  }
 ],
 "metadata": {},
 "nbformat": 4,
 "nbformat_minor": 5
}
EOF

# -------- Lizenz --------------------------------------------
cat > LICENSE <<'EOF'
Creative Commons Attribution 4.0 International
<https://creativecommons.org/licenses/by/4.0/>
EOF

# -------- Git-Repo initialisieren ----------------------------
git init -q
git add .
git commit -qm "Initial scaffold for ${REPO_NAME} v${VERSION}"

echo
echo "✅  Projekt »${REPO_NAME}« wurde angelegt und initial committet."
echo "   Als Nächstes (ersetzt <URL> durch dein GitHub-Repo):"
echo "   git remote add origin <URL>"
echo "   git push -u origin main"
echo
