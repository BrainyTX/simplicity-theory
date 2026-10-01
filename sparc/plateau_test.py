"""
plateau_test.py – does the plateau velocity follow v^4 = G * M * 2*pi*G*Sigma*?

Past the threshold the rotation curve stays flat (the "torque plateau"); its
height should then depend only on the total baryonic mass M and the threshold
(baryonic Tully-Fisher relation). Prediction without any free parameter:

    RAR / MOND:         v_flat^4 = G * M * g_dagger,  g_dagger = 2*pi*G*Sigma*
    three-phase model:  v_flat^4 = c^2 * G * M * g_dagger,
                        c = (2 - x_low^2) * sqrt(x_low) = 0.990

Plateau galaxies: at least 6 points, last three points vary by < 10 %.
M is estimated from the baryonic rotation curve at the last measured radius
(M = V_bar^2 R / G). This is rough: gas beyond the last point is missing.
The SPARC galaxy table (SPARC_Lelli2016c.mrt, luminosity and HI mass) would
give better masses.

Usage:
    python sparc/plateau_test.py
Output:
    results/sparc_analysis/plateau_test.csv
    results/sparc_analysis/plateau_test.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sparc_tools as st

FLATNESS = 0.10
MIN_POINTS = 6


def plateau_galaxies(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gal, g in df.groupby("galaxy"):
        g = g.sort_values("R")
        if len(g) < MIN_POINTS:
            continue
        last = g.tail(3)
        v_flat = last["Vobs"].mean()
        if abs(last["Vobs"].iat[-1] - last["Vobs"].iat[0]) / v_flat >= FLATNESS:
            continue
        rows.append({
            "galaxy": gal,
            "v_flat": v_flat,
            "e_v_flat": last["eV"].mean(),
            "M_bar": g["Vbar2"].iat[-1] * g["R"].iat[-1] / st.G_KPC,
            "R_last": g["R"].iat[-1],
        })
    return pd.DataFrame(rows)


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t = plateau_galaxies(st.load_sparc())
    t.to_csv(st.RESULTS_DIR / "plateau_test.csv", index=False, float_format="%.4g")
    print(f"{len(t)} galaxies with a flat plateau\n")

    x_low = st.LOW_TO_STAR_RATIO
    c_tp = (2 - x_low ** 2) * np.sqrt(x_low)
    predictions = {
        f"RAR/MOND, Σ* = {st.SIGMA_STAR:.2f}": (st.SIGMA_STAR, 1.0, "tab:red", "--"),
        "RAR/MOND, Σ* = 137": (137.0, 1.0, "tab:purple", ":"),
        f"three-phase, Σ* = {st.SIGMA_STAR:.2f}": (st.SIGMA_STAR, c_tp ** 2, "tab:blue", "-"),
    }
    log_v4 = np.log10(t["v_flat"] ** 4)
    for name, (sigma, coeff, *_) in predictions.items():
        r = log_v4 - np.log10(coeff * st.G_KPC * t["M_bar"] * st.g_dagger(sigma))
        print(f"{name:32s} median log(v^4 observed / predicted) = {np.median(r):+.3f} dex "
              f"(factor {10 ** np.median(r):.2f} in v^4, {10 ** (np.median(r) / 4):.3f} in v), "
              f"scatter {np.std(r):.2f} dex")

    g_best = 10 ** np.median(log_v4 - np.log10(st.G_KPC * t["M_bar"]))
    slope, _ = np.polyfit(np.log10(t["M_bar"]), np.log10(t["v_flat"]), 1)
    print(f"\nSigma* preferred by the plateau heights: {g_best / st.g_dagger(1.0):.1f} M_sun/pc^2")
    print(f"Slope of log v_flat vs log M: {slope:.3f} (plateau rule expects 0.25)")

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.errorbar(t["M_bar"], t["v_flat"], yerr=t["e_v_flat"], fmt="o", ms=4, color="0.35",
                ecolor="0.7", label=f"SPARC plateau galaxies (n={len(t)})")
    m = np.logspace(np.log10(t["M_bar"].min()) - 0.2, np.log10(t["M_bar"].max()) + 0.2, 100)
    for name, (sigma, coeff, col, ls) in predictions.items():
        ax.plot(m, (coeff * st.G_KPC * m * st.g_dagger(sigma)) ** 0.25, color=col, ls=ls, lw=2,
                label=f"v⁴ = {'' if coeff == 1 else f'{coeff:.3f}·'}G·M·2πGΣ*  ({name})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("baryonic mass M_bar  [M☉]  (from V_bar at the last radius)")
    ax.set_ylabel("plateau velocity v_flat  [km/s]")
    ax.set_title("Plateau test (baryonic Tully-Fisher): no free parameter")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = st.RESULTS_DIR / "plateau_test.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\nSaved {out.relative_to(st.REPO_ROOT)} and plateau_test.csv")


if __name__ == "__main__":
    main()
