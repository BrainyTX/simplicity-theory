"""
fit_sigma_crit.py – which Sigma*_crit do the SPARC rotation curves prefer?

1. Per galaxy: fits Sigma* to each rotation curve (68 % range from delta chi^2 = 1).
2. Global: one Sigma* for all galaxies together, with a bootstrap over galaxies.
3. Systematics: repeats the global fit for disk M/L = 0.4, 0.5, 0.6.

This replaces the fit in rar_fit.ipynb, which used M/L = 1 (overestimating the
baryons) and dropped the sign of negative gas velocities.

Usage:
    python sparc/fit_sigma_crit.py                 # RAR form (default)
    python sparc/fit_sigma_crit.py --model mond    # MOND simple-mu form
Output:
    results/sparc_analysis/sigma_crit_fits.csv
    results/sparc_analysis/sigma_crit_fits.png
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar

import sparc_tools as st

MODEL_CHOICES = {"rar": st.rar, "mond": st.mond_simple}
LOG_BOUNDS = (0.0, 3.5)          # Sigma between 1 and ~3000 M_sun/pc^2
MIN_POINTS = 5


def fit_one_galaxy(g: pd.DataFrame, model) -> dict:
    """Best Sigma* for one galaxy and its 68 % range (delta chi^2 = 1)."""
    def chi2(log_sigma):
        return st.chi2_velocity(g, model(g, sigma_star=10 ** log_sigma))

    best = minimize_scalar(chi2, bounds=LOG_BOUNDS, method="bounded")
    chi2_min = best.fun

    def edge(side):
        limit = LOG_BOUNDS[0] if side < 0 else LOG_BOUNDS[1]
        f = lambda ls: chi2(ls) - chi2_min - 1.0
        if f(limit) < 0:
            return np.nan                  # range hits the bound
        return 10 ** brentq(f, *sorted((best.x, limit)))

    n = len(g)
    return {
        "galaxy": g["galaxy"].iat[0],
        "n_points": n,
        "sigma_fit": 10 ** best.x,
        "sigma_lo": edge(-1),
        "sigma_hi": edge(+1),
        "chi2": chi2_min,
        "reduced_chi2": chi2_min / (n - 1),
        "at_bound": bool(np.isclose(best.x, LOG_BOUNDS, atol=1e-3).any()),
    }


def fit_global(df: pd.DataFrame, model) -> float:
    """One Sigma* for all points (log residuals with systematic floor)."""
    res = minimize_scalar(lambda ls: st.chi2_global(df, model(df, sigma_star=10 ** ls)),
                          bounds=LOG_BOUNDS, method="bounded")
    return 10 ** res.x


def bootstrap_global(df: pd.DataFrame, model, n_boot: int, seed: int = 1) -> np.ndarray:
    """Bootstrap over galaxies (not points): galaxies are resampled as a whole."""
    rng = np.random.default_rng(seed)
    groups = [g for _, g in df.groupby("galaxy")]
    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        out.append(fit_global(pd.concat([groups[i] for i in pick]), model))
    return np.array(out)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--model", choices=MODEL_CHOICES, default="rar")
    p.add_argument("--n-boot", type=int, default=200, help="bootstrap samples (default 200)")
    args = p.parse_args()
    model = MODEL_CHOICES[args.model]
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    df = st.load_sparc()
    print(f"Loaded {len(df)} points from {df['galaxy'].nunique()} galaxies "
          f"(M/L disk {st.ML_DISK}, bulge {st.ML_BULGE}), model: {args.model}\n")

    # 1) per galaxy
    rows = [fit_one_galaxy(g, model) for _, g in df.groupby("galaxy") if len(g) >= MIN_POINTS]
    fits = pd.DataFrame(rows).sort_values("sigma_fit")
    out_csv = st.RESULTS_DIR / "sigma_crit_fits.csv"
    fits.to_csv(out_csv, index=False, float_format="%.4g")
    free = fits[~fits["at_bound"]]
    q16, q50, q84 = np.percentile(free["sigma_fit"], [16, 50, 84])
    print(f"Per galaxy ({len(fits)} fitted, {fits['at_bound'].sum()} at the fit bound):")
    print(f"  median Sigma = {q50:.1f}   (16-84 %: {q16:.1f} - {q84:.1f}) M_sun/pc^2")

    # 2) global with bootstrap
    s_glob = fit_global(df, model)
    boot = bootstrap_global(df, model, args.n_boot)
    b16, b84 = np.percentile(boot, [16, 84])
    print(f"\nGlobal (one Sigma for all galaxies): {s_glob:.1f}   "
          f"(68 % bootstrap: {b16:.1f} - {b84:.1f}) M_sun/pc^2")

    # 3) M/L systematics
    print("\nSystematics, global fit for different disk M/L (bulge = 1.4 x disk):")
    for ml in (0.4, 0.5, 0.6):
        s = fit_global(st.load_sparc(ml_disk=ml, ml_bulge=1.4 * ml), model)
        print(f"  M/L disk {ml}: Sigma = {s:.1f}")

    h0_lo, h0_hi = st.sigma_star_from_h0(67.4), st.sigma_star_from_h0(73.0)
    print(f"\nFor comparison: codex Sigma* = {st.SIGMA_STAR:.2f}, "
          f"c*H0/(2pi)^2G = {h0_lo:.1f} (H0=67.4) ... {h0_hi:.1f} (H0=73.0)")

    # plot
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bins = np.logspace(0, 3.5, 36)
    ax.hist(free["sigma_fit"], bins=bins, color="0.6", edgecolor="0.3",
            label=f"per-galaxy fits (n={len(free)})")
    ax.axvspan(b16, b84, color="tab:green", alpha=0.25,
               label=f"global fit {s_glob:.0f} (68 %: {b16:.0f}-{b84:.0f})")
    ax.axvspan(h0_lo, h0_hi, color="tab:orange", alpha=0.25,
               label=f"c·H0/(2π)²G for H0 = 67.4-73 ({h0_lo:.0f}-{h0_hi:.0f})")
    ax.axvline(st.SIGMA_STAR, color="tab:red", lw=2,
               label=f"Σ* = 137 − Ω_φ − 1/4 = {st.SIGMA_STAR:.2f}")
    ax.axvline(137.0, color="tab:purple", ls="--", lw=1.5, label="137")
    ax.set_xscale("log")
    ax.set_xlabel("fitted Σ*  [M☉/pc²]")
    ax.set_ylabel("number of galaxies")
    ax.set_title(f"Σ* preferred by the SPARC rotation curves ({args.model.upper()} form)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out_png = st.RESULTS_DIR / "sigma_crit_fits.png"
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"\nSaved {out_csv.relative_to(st.REPO_ROOT)} and {out_png.relative_to(st.REPO_ROOT)}")


if __name__ == "__main__":
    main()
