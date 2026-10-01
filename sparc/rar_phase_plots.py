"""
rar_phase_plots.py – the two overview figures for all SPARC data points.

1. Radial acceleration relation (RAR): observed vs. baryonic acceleration for
   every measured radius of every galaxy, with the phase boundaries
   g = 2*pi*G*Sigma_low and g = 2*pi*G*Sigma* and the model curves.
2. Phase diagram: required amplification g_obs/g_bar vs. surface density,
   once for Sigma_eff = g_bar/(2 pi G) (includes gas) and once for the local
   stellar Sigma (what the original Square-Rule uses).

No model is fitted here: Sigma* = 123.69 and Sigma_low = 32.54 are fixed.

Usage:
    python sparc/rar_phase_plots.py
Output:
    results/sparc_analysis/rar.png
    results/sparc_analysis/phase_diagram.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sparc_tools as st

PHASES = [  # (label, lower Sigma, upper Sigma, colour)
    ("gas", 0.0, st.SIGMA_LOW, "tab:blue"),
    ("liquid", st.SIGMA_LOW, st.SIGMA_STAR, "tab:cyan"),
    ("solid", st.SIGMA_STAR, np.inf, "tab:purple"),
]


def shade_phases(ax, scale=1.0, lo=1e-30, hi=1e30, label_y=None):
    """Shades the three phases along the x axis; scale converts Sigma to x units."""
    for name, a, b, col in PHASES:
        x0, x1 = max(a * scale, lo), min(b * scale, hi)
        ax.axvspan(x0, x1, color=col, alpha=0.07)
        if label_y is not None:
            ax.text(np.sqrt(x0 * x1) if np.isfinite(x1) else x0 * 3, label_y, name,
                    color=col, ha="center", fontsize=10, weight="bold")
    for s, col in [(st.SIGMA_LOW, "tab:blue"), (st.SIGMA_STAR, "tab:purple")]:
        ax.axvline(s * scale, color=col, ls="--", lw=1)


def rar_figure(df):
    si = st.ACC_TO_SI
    gb, go = df["g_bar"].to_numpy() * si, df["g_obs"].to_numpy() * si
    g_curve = np.logspace(np.log10(gb.min()), np.log10(gb.max()), 400)
    curve_df = pd.DataFrame({"g_bar": g_curve / si})
    rar_curve = st.rar(curve_df) * si
    tp_curve = curve_df["g_bar"] * st.three_phase_boost(curve_df["g_bar"] / st.g_dagger(st.SIGMA_STAR)) * si
    edges = np.logspace(np.log10(gb.min()), np.log10(gb.max()), 25)
    c, m, _ = st.binned_median(gb, np.log10(go), edges)
    scale = st.g_dagger(1.0) * si          # Sigma [M_sun/pc^2] -> g [m/s^2]

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8.5, 9), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    shade_phases(ax, scale, gb.min(), gb.max() * 2, label_y=go.max() * 0.4)
    ax.scatter(gb, go, s=3, color="0.45", alpha=0.25, label=f"SPARC ({len(df)} points, "
               f"{df['galaxy'].nunique()} galaxies)")
    ax.plot(c, 10 ** m, "ko", ms=5, label="median per bin")
    ax.plot(g_curve, g_curve, color="0.3", lw=1, ls=":", label="Newton (g_obs = g_bar)")
    ax.plot(g_curve, rar_curve, color="tab:red", lw=2, ls="--", label="RAR, g† = 2πGΣ*")
    ax.plot(g_curve, tp_curve, color="tab:blue", lw=2, label="three-phase model")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("observed acceleration g_obs  [m/s²]")
    ax.set_title(f"Radial acceleration relation  (Σ* = {st.SIGMA_STAR:.2f}, "
                 f"Σ_low = {st.SIGMA_LOW:.2f} M☉/pc²)")
    ax.legend(fontsize=8, loc="lower right")

    r = np.log10(go) - np.log10(st.rar(df) * si)
    shade_phases(axr, scale, gb.min(), gb.max() * 2)
    axr.scatter(gb, r, s=3, color="0.45", alpha=0.25)
    cr, mr, _ = st.binned_median(gb, r, edges)
    axr.plot(cr, mr, "ko-", ms=4)
    axr.axhline(0, color="tab:red", ls="--")
    axr.set_ylim(-0.8, 0.8)
    axr.set_xlabel("baryonic acceleration g_bar  [m/s²]")
    axr.set_ylabel("log(g_obs / RAR)")
    fig.tight_layout()
    return fig


def phase_figure(df):
    boost = (df["g_obs"] / df["g_bar"]).to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    panels = [
        ("sigma_eff", "effective surface density Σ_eff = g_bar / 2πG  [M☉/pc²]",
         "phase set by Σ_eff (gas, disk and bulge)"),
        ("sigma_stars", "local stellar surface density Σ_stars  [M☉/pc²]",
         "phase set by the local stellar Σ (original Square-Rule)"),
    ]
    for ax, (col, xlabel, title) in zip(axes, panels):
        sigma = df[col].to_numpy()
        ok = sigma > 0
        s, b = sigma[ok], boost[ok]
        edges = np.logspace(np.log10(s.min()), np.log10(s.max()), 22)
        c, m, _ = st.binned_median(s, b, edges, min_count=15)
        shade_phases(ax, 1.0, s.min(), s.max() * 1.5, label_y=0.62)
        ax.scatter(s, b, s=3, color="0.45", alpha=0.2, label=f"SPARC points (n={ok.sum()})")
        ax.plot(c, m, "ko-", lw=2, ms=4, label="median per bin")
        ax.axhline(1, color="k", lw=0.8)
        ax.axhline(2, color="tab:green", ls=":", label="maximum of the original Square-Rule")
        if col == "sigma_eff":
            xs = np.logspace(np.log10(s.min()), np.log10(s.max()), 400)
            ax.plot(xs, 1 / (1 - np.exp(-np.sqrt(xs / st.SIGMA_STAR))), "r--", lw=2, label="RAR")
            ax.plot(xs, st.three_phase_boost(xs / st.SIGMA_STAR), color="tab:blue", lw=2,
                    label="three-phase model")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(0.55, 40)
        ax.set_xlim(s.min(), s.max() * 1.5)
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="upper right")
    axes[0].set_ylabel("required amplification  g_obs / g_bar")
    fig.suptitle("Phase diagram: how much must gravity be amplified at a given density?")
    fig.tight_layout()
    return fig


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df = st.load_sparc()
    for name, fig in [("rar.png", rar_figure(df)), ("phase_diagram.png", phase_figure(df))]:
        out = st.RESULTS_DIR / name
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"Saved {out.relative_to(st.REPO_ROOT)}")


if __name__ == "__main__":
    main()
