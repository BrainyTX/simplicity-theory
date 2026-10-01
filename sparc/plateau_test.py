"""
plateau_test.py – does the plateau velocity follow v^4 = G * M * 2*pi*G*Sigma*?

Past the threshold the rotation curve stays flat (the "torque plateau"); its
height should then depend only on the total baryonic mass M and the threshold
(baryonic Tully-Fisher relation, BTFR). Prediction without a free parameter:

    RAR / MOND:         v_flat^4 = G * M * g_dagger,  g_dagger = 2*pi*G*Sigma*
    three-phase model:  v_flat^4 = c^2 * G * M * g_dagger,
                        c = (2 - x_low^2) * sqrt(x_low) = 0.990

Data: data/sparc/BTFR_Lelli2016a.mrt (Lelli, McGaugh & Schombert 2016), 118
galaxies with baryonic mass (stars with M/L 0.5 plus gas), flat velocity,
gas fraction and effective stellar surface density SBeff.

The script reports
  1. the offset of the data from the prediction for Sigma* = 123.69 and 137,
  2. the Sigma* preferred by the plateau heights (slope fixed to 4), with a
     bootstrap and the M/L systematics, separately for gas-rich galaxies
     (Fg > 0.5, hardly affected by M/L) and star-dominated ones (Fg < 0.2),
  3. the free slope (orthogonal distance regression with errors in both axes),
  4. whether the plateau depends on the "phase" of the galaxy (SBeff below
     Sigma_low, between, above Sigma*) at fixed mass.

Usage:
    python sparc/plateau_test.py
Output:
    results/sparc_analysis/plateau_test.csv
    results/sparc_analysis/plateau_test.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import odr, stats

import sparc_tools as st

PHASE_EDGES = [0, st.SIGMA_LOW, st.SIGMA_STAR, np.inf]
PHASE_NAMES = ["gas", "liquid", "solid"]
PHASE_COLOURS = ["tab:blue", "tab:cyan", "tab:purple"]


def log_mb_for_ml(t: pd.DataFrame, ml: float) -> pd.Series:
    """log10 baryonic mass for another stellar M/L (the table uses 0.5)."""
    return t["logMb"] + np.log10(t["Fg"] + (1 - t["Fg"]) * ml / 0.5)


def offset(t, log_mb, sigma, coeff=1.0):
    """Weighted mean of log10(v^4 observed / v^4 predicted) and its residuals."""
    r = 4 * t["logVf"] - (np.log10(coeff * st.G_KPC * st.g_dagger(sigma)) + log_mb)
    w = 1 / np.hypot(4 * t["e_logVf"], t["e_logMb"]) ** 2
    return float(np.sum(w * r) / w.sum()), r


def best_sigma(t, log_mb):
    """Sigma* that makes the weighted mean offset zero (slope fixed to 4)."""
    return 10 ** offset(t, log_mb, 1.0)[0]


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t = st.load_btfr()
    print(f"{len(t)} galaxies from {st.BTFR_FILE.name} (stellar M/L 0.5)\n")

    # 1) offsets from the parameter-free predictions
    x_low = st.LOW_TO_STAR_RATIO
    c_tp = (2 - x_low ** 2) * np.sqrt(x_low)
    for name, sigma, coeff in [(f"RAR/MOND, Sigma* = {st.SIGMA_STAR:.2f}", st.SIGMA_STAR, 1.0),
                               ("RAR/MOND, Sigma* = 137", 137.0, 1.0),
                               (f"three-phase, Sigma* = {st.SIGMA_STAR:.2f}", st.SIGMA_STAR, c_tp ** 2)]:
        m, r = offset(t, t["logMb"], sigma, coeff)
        print(f"{name:32s} log(v^4 observed/predicted) = {m:+.3f} dex "
              f"(factor {10 ** m:.2f} in v^4, {10 ** (m / 4):.3f} in v), scatter {np.std(r):.2f} dex")

    # 2) preferred Sigma* with bootstrap and M/L systematics
    rng = np.random.default_rng(3)
    s_best = best_sigma(t, t["logMb"])
    boot = [best_sigma(t.iloc[i], t["logMb"].iloc[i])
            for i in (rng.integers(0, len(t), len(t)) for _ in range(2000))]
    print(f"\nSigma* preferred by the plateau heights (slope 4): {s_best:.1f} "
          f"(68 %: {np.percentile(boot, 16):.1f} - {np.percentile(boot, 84):.1f})")
    gas_rich, star_rich = t["Fg"] > 0.5, t["Fg"] < 0.2
    print(f"M/L systematics      all ({len(t)})   gas-rich ({gas_rich.sum()})   star-dominated ({star_rich.sum()})")
    for ml in (0.3, 0.4, 0.5, 0.6):
        lm = log_mb_for_ml(t, ml)
        print(f"  M/L {ml}:          {best_sigma(t, lm):7.1f}   {best_sigma(t[gas_rich], lm[gas_rich]):12.1f}"
              f"   {best_sigma(t[star_rich], lm[star_rich]):14.1f}")

    # 3) free slope: log Vf = a + b log Mb, errors in both axes
    line = lambda beta, x: beta[0] + beta[1] * x
    fit = odr.ODR(odr.RealData(t["logMb"], t["logVf"], sx=t["e_logMb"], sy=t["e_logVf"]),
                  odr.Model(line), beta0=[0.0, 0.25]).run()
    b, e_b = fit.beta[1], fit.sd_beta[1]
    print(f"\nFree slope of log Vf vs log Mb: {b:.3f} +/- {e_b:.3f} "
          f"(plateau rule 0.25, i.e. Mb ~ Vf^{1 / b:.2f} instead of Vf^4; "
          f"{abs(b - 0.25) / e_b:.1f} sigma from 0.25)")

    # 4) does the phase of the galaxy matter at fixed mass?
    t["resid_logVf"] = t["logVf"] - line(fit.beta, t["logMb"])
    t["phase"] = pd.cut(t["SBeff"], PHASE_EDGES, labels=PHASE_NAMES)
    rho, p = stats.spearmanr(np.log10(t["SBeff"]), t["resid_logVf"])
    print(f"Residual (at fixed mass) vs. SBeff: Spearman {rho:+.2f}, p = {p:.2f}")
    print(t.groupby("phase", observed=True)["resid_logVf"]
          .agg(n="size", mean="mean", sem="sem").round(4).to_string())
    t.to_csv(st.RESULTS_DIR / "plateau_test.csv", index=False, float_format="%.4g")

    # plot
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(14, 5.8))
    for name, col in zip(PHASE_NAMES, PHASE_COLOURS):
        s = t[t["phase"] == name]
        ax.errorbar(s["logMb"], s["logVf"], xerr=s["e_logMb"], yerr=s["e_logVf"], fmt="o", ms=4,
                    color=col, ecolor=col, alpha=0.8, elinewidth=0.8, label=f"{name} (SBeff, n={len(s)})")
    m = np.linspace(t["logMb"].min() - 0.2, t["logMb"].max() + 0.2, 50)
    for label, sigma, ls, col in [(f"Σ* = {st.SIGMA_STAR:.2f}", st.SIGMA_STAR, "--", "tab:red"),
                                  ("Σ* = 137", 137.0, ":", "k"),
                                  (f"best Σ* = {s_best:.0f} (slope 4)", s_best, "-", "0.5")]:
        ax.plot(m, 0.25 * (np.log10(st.G_KPC * st.g_dagger(sigma)) + m), ls=ls, color=col, lw=2,
                label=f"v⁴ = G·M·2πGΣ*, {label}")
    ax.plot(m, line(fit.beta, m), color="tab:green", lw=1.5, label=f"free fit, slope {b:.3f}")
    ax.set_xlabel("log baryonic mass Mb  [M☉]")
    ax.set_ylabel("log plateau velocity Vf  [km/s]")
    ax.set_title("Plateau test (baryonic Tully-Fisher), Lelli et al. 2016")
    ax.legend(fontsize=7.5, loc="upper left")

    for (lo, hi), col in zip(zip(PHASE_EDGES[:-1], PHASE_EDGES[1:]), PHASE_COLOURS):
        axr.axvspan(max(lo, 1), min(hi, 1e4), color=col, alpha=0.08)
    axr.axvline(st.SIGMA_LOW, color="tab:blue", ls="--", lw=1)
    axr.axvline(st.SIGMA_STAR, color="tab:purple", ls="--", lw=1)
    axr.errorbar(t["SBeff"], t["resid_logVf"], yerr=t["e_logVf"], fmt="o", ms=4, color="0.35",
                 ecolor="0.7", elinewidth=0.8)
    axr.axhline(0, color="k", lw=0.8)
    axr.set_xscale("log")
    axr.set_xlabel("effective stellar surface density SBeff  [M☉/pc²]")
    axr.set_ylabel("plateau residual at fixed mass  [dex in Vf]")
    axr.set_title(f"Does the phase of the galaxy matter?  Spearman {rho:+.2f}, p = {p:.2f}")
    fig.tight_layout()
    out = st.RESULTS_DIR / "plateau_test.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\nSaved {out.relative_to(st.REPO_ROOT)} and plateau_test.csv")


if __name__ == "__main__":
    main()
