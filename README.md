# simplicity-theory

This thesis proposes a revolutionary framework for understanding physical reality.

This repository contains the computational companion material for the
*Simplicity Theory*: symbolic derivations for the appendices, the
calculation of the critical surface density Σ\*<sub>crit</sub>, tests of
the model against the SPARC galaxy rotation curves, and several
illustrative simulations.

## Repository layout

```
simplicity-theory/
├── theory/                  Symbolic derivations for the paper's appendices (SymPy)
│   ├── appendix_s1_overlap_geometry.py   S1: FCC sphere overlap → unique ξ = H/R = 1/2
│   ├── appendix_s2_pi_lock_in.py         S2: π is locked into the eigen-volume
│   └── appendix_s3_lagrangian.py         S3: parameter-free Lagrangian, wave speed c₀
│
├── sigma_crit/              Critical surface density Σ*_crit = c·H₀ / ((2π)²·G)
│   ├── sigma_crit.py                     CLI calculator (Σ*_crit ≈ 124 M☉/pc² for H₀ = 70)
│   └── sigma_crit_demo.ipynb             Robustness of Σ*_crit for H₀ = 67…74 km/s/Mpc
│
├── sparc/                   Tests against the SPARC rotation-curve catalogue
│   ├── square_rule_basic.py              Square-Rule f(Σ) = 1 − (Σ/Σ*)² applied to all galaxies
│   ├── square_rule_smbh.py               Extended Square-Rule (bulge, SMBH, age shift) → per-galaxy plots
│   ├── rar_fit.ipynb                     Fits Σ_crit per galaxy, writes results/sparc_sigma_fit_results.csv
│   ├── rar_threshold_plot.py             Surface density vs. dark-matter fraction (RAR threshold figure)
│   └── dark_hole.py                      "Dark-hole" radius r_DH = √(M / 4πΣ*) for known black holes
│
├── simulations/             Illustrative simulations
│   ├── vacuum_oscillation.py             3D "Urenergie" interference: cancellation of vacuum oscillations
│   └── double_slit_transition.py         Double-slit wave→particle transition (animated GIF)
│
├── data/sparc/              Input data
│   ├── Rotmod_LTG/                       152 SPARC mass-model files (*_rotmod.dat)
│   └── SPARC_summary.csv                 μ₀ and dark-matter fraction f_DM per galaxy
│
├── results/                 Generated output
│   ├── sparc_galaxy_fits/                91 rotation-curve analyses (*_analyse.png) from square_rule_smbh.py
│   └── double_slit_transition.gif        Output of double_slit_transition.py
│
├── paper/                   LaTeX snippets for the manuscript
│   └── sparc_threshold_figure.tex
│
└── archive/                 Old helper scripts, kept for reference
    └── scaffold.sh                       Original project scaffold generator
```

## Getting started

Requires Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

All scripts locate their data relative to the repository root, so they can be
run from any directory. Outputs are written to `results/`.

```bash
# Appendix derivations
python theory/appendix_s1_overlap_geometry.py
python theory/appendix_s2_pi_lock_in.py
python theory/appendix_s3_lagrangian.py

# Critical surface density (default H0 = 70 km/s/Mpc, or pass your own)
python sigma_crit/sigma_crit.py
python sigma_crit/sigma_crit.py 73.0

# SPARC analyses
python sparc/square_rule_basic.py       # opens one plot window per galaxy
python sparc/square_rule_smbh.py        # writes results/sparc_galaxy_fits/*.png
python sparc/rar_threshold_plot.py      # writes results/RAR_threshold.pdf
python sparc/dark_hole.py

# Simulations
python simulations/vacuum_oscillation.py       # writes results/urenergie_extended_3d_interference.png
python simulations/double_slit_transition.py   # writes results/double_slit_transition.gif (takes a while)
```

The notebooks are opened with Jupyter from inside their folder:

```bash
jupyter lab
```

## Notes

- **SPARC metadata:** `sparc/square_rule_smbh.py` can additionally read the
  SPARC galaxy table (`Table1.mrt`, saved as
  `data/sparc/SPARC_Lelli2016c.txt`). This file is not included in the
  repository; without it the script still runs, but without the bulge/SMBH
  metadata. Running it will overwrite the plots in `results/sparc_galaxy_fits/`.
- **Truncated scripts:** the end of `simulations/vacuum_oscillation.py` was
  cut off in the original upload; the remaining analysis output is missing.

## Data source

The rotation curves in `data/sparc/Rotmod_LTG/` are from the SPARC database:

> F. Lelli, S. S. McGaugh, J. M. Schombert (2016),
> *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry at 3.6 μm
> and Accurate Rotation Curves*, AJ 152, 157.
> <http://astroweb.cwru.edu/SPARC/>

## License

© Pieter Goldau. Licensed under
[Creative Commons Attribution 4.0 International](LICENSE) (CC BY 4.0).
