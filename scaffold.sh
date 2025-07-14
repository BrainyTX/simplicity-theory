#!/usr/bin/env bash
set -e                              # sofort stoppen bei Fehlern

####################################
# 0. Metadaten – einmal anpassen   #
####################################
REPO_NAME="simplicity-theory"
AUTHOR="Pieter Goldau"
EMAIL="simplicity-theory@gmx.de"
VERSION="0.1.0"

####################################
# 1. Ordnerstruktur erzeugen       #
####################################
mkdir -p "${REPO_NAME}"/{src,docs,notebooks,data,tests}
cd "${REPO_NAME}"

####################################
# 2. README.md                     #
####################################
cat > README.md <<EOF
# ${REPO_NAME}

Mathematische Tools zur **Existential Coherence Theory**  
Version ${VERSION} – © ${AUTHOR}

```bash
# Installation (virtuelle Umgebung empfohlen)
pip install -r requirements.txt

# Schnelle Prüfung der Vorhersage bei H₀ = 70 km/s/Mpc
python src/sigma_crit.py --H0 70
