"""
compare_models.py – how well do the gravity models describe the SPARC data?

All models use the same data, the same mass model (M/L 0.5 / 0.7) and the same
Sigma* = 123.69 M_sun/pc^2; none has a free parameter per galaxy.

    Newton                 baryons only
    Square-Rule (original) g = g_bar (1 + f), f = 1 - (Sigma_stars/Sigma*)^2
    RAR (McGaugh 2016)     g = g_bar / (1 - exp(-sqrt(g_bar/g_dagger)))
    MOND (simple mu)       g = g_bar/2 + sqrt(g_bar^2/4 + g_bar a0)
    Three-phase            solid (Newton) / liquid (Square-Rule) / gas,
                           thresholds Sigma* = 123.69 and Sigma_low = 32.54,
                           phase set by Sigma_eff = g_bar / (2 pi G)
    Three-phase hybrid     Square-Rule on the local stellar Sigma for
                           liquid/solid, gas law on Sigma_eff

Additionally:
  * the gas exponent of the three-phase model is fitted globally
    (0.5 = flat rotation curves),
  * a kink test checks whether the data change slope at the phase boundaries
    more than at arbitrary other places (sharp phase transition or smooth?).

Usage:
    python sparc/compare_models.py
Output:
    results/sparc_analysis/model_comparison.csv      (summary per model)
    results/sparc_analysis/model_comparison_galaxies.csv  (chi^2 per galaxy)
    results/sparc_analysis/model_residuals.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import sparc_tools as st


def summarize(df, models):
    rows, per_gal = [], {}
    for name, g_pred in models.items():
        r = st.log_residuals(df, g_pred)
        rows.append({
            "model": name,
            "median_residual_dex": np.median(r),
            "scatter_dex": np.std(r),
            "chi2_global": st.chi2_global(df, g_pred),
        })
        per_gal[name] = {gal: st.chi2_velocity(g, g_pred[g.index.map(df.index.get_loc)])
                         for gal, g in df.groupby("galaxy")}
    summary = pd.DataFrame(rows)
    galaxies = pd.DataFrame(per_gal)
    galaxies.index.name = "galaxy"
    # "best model" is counted among the base models only (fitted variants excluded)
    best = galaxies[[m for m in st.MODELS if m in galaxies]].idxmin(axis=1).value_counts()
    is_base = summary["model"].isin(st.MODELS)
    summary["best_for_n_galaxies"] = (summary["model"].map(best)
                                      .where(~is_base, summary["model"].map(best).fillna(0))
                                      .astype("Int64"))
    summary["median_galaxy_chi2"] = summary["model"].map(galaxies.median())
    return summary, galaxies


def kink_test(df, n_boot=200, seed=0):
    """Slope change of the RAR residuals at the two phase boundaries.

    A sharp phase transition shows up as a kink in log(g_obs) vs log(x).
    Because the residuals are not perfectly straight, the same fit is repeated
    with the boundaries shifted to other positions (control): only a kink that
    is clearly larger there than elsewhere points to a real phase boundary.
    """
    x = (df["g_bar"] / st.g_dagger(st.SIGMA_STAR)).to_numpy()
    lx = np.log10(x)
    r = st.log_residuals(df, st.rar(df))
    b_low, b_star = np.log10(st.LOW_TO_STAR_RATIO), 0.0

    def kinks(b1, b2, idx=slice(None)):
        X = np.column_stack([np.ones_like(lx), lx,
                             np.maximum(0, lx - b1), np.maximum(0, lx - b2)])[idx]
        c, *_ = np.linalg.lstsq(X, r[idx], rcond=None)
        return c[2:]

    k = kinks(b_low, b_star)
    rng = np.random.default_rng(seed)
    groups = list(df.groupby("galaxy").indices.values())
    boot = np.array([kinks(b_low, b_star, np.concatenate(
        [groups[i] for i in rng.integers(0, len(groups), len(groups))])) for _ in range(n_boot)])
    shifts = [s for s in np.linspace(-1.5, 1.0, 26) if abs(s) > 0.05]
    control = np.abs([kinks(b_low + s, b_star + s) for s in shifts])
    return k, np.percentile(boot, [2.5, 97.5], axis=0), np.median(control, axis=0)


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df = st.load_sparc()
    print(f"{len(df)} points, {df['galaxy'].nunique()} galaxies, Sigma* = {st.SIGMA_STAR:.2f}, "
          f"Sigma_low = {st.SIGMA_LOW:.2f} M_sun/pc^2\n")

    models = {name: f(df) for name, f in st.MODELS.items()}

    # three-phase model with fitted gas exponent
    res = minimize_scalar(lambda n: st.chi2_global(df, st.three_phase(df, gas_exponent=n)),
                          bounds=(0.1, 1.5), method="bounded")
    n_fit = res.x
    models[f"Three-phase (gas exponent {n_fit:.2f})"] = st.three_phase(df, gas_exponent=n_fit)

    summary, galaxies = summarize(df, models)
    pd.set_option("display.width", 140)
    print(summary.to_string(index=False, float_format="%.3f"))
    summary.to_csv(st.RESULTS_DIR / "model_comparison.csv", index=False, float_format="%.4g")
    galaxies.to_csv(st.RESULTS_DIR / "model_comparison_galaxies.csv", float_format="%.4g")

    deep = (2 - st.LOW_TO_STAR_RATIO ** 2) * np.sqrt(st.LOW_TO_STAR_RATIO)
    print(f"\nThree-phase deep gas limit: g = {deep:.3f} * sqrt(g_bar * g_dagger)  "
          f"(MOND/RAR: 1.000)")

    # residuals by regime
    x = (df["g_bar"] / st.g_dagger(st.SIGMA_STAR)).to_numpy()
    edges = [x.min(), st.LOW_TO_STAR_RATIO / 10, st.LOW_TO_STAR_RATIO, 0.6, 1, 2, 5, x.max() * 1.01]
    regimes = pd.cut(x, edges)
    table = pd.DataFrame({name: st.log_residuals(df, g) for name, g in models.items()
                          if name in ("Square-Rule (original)", "RAR (McGaugh 2016)",
                                      "Three-phase", "Three-phase hybrid")})
    print("\nMedian residual log(g_obs/g_pred) [dex] by x = Sigma_eff/Sigma*  "
          "(positive = model too weak):")
    print(table.groupby(regimes, observed=True).median()
          .assign(n=table.groupby(regimes, observed=True).size()).round(3).to_string())

    # kink test
    k, ci, ctrl = kink_test(df)
    print("\nKink test (slope change of RAR residuals at the phase boundaries):")
    for i, name in enumerate(["lower boundary 32.5", "upper boundary 123.7"]):
        sig = not (ci[0, i] < 0 < ci[1, i])
        print(f"  {name}: {k[i]:+.3f}  95 % [{ci[0, i]:+.3f}, {ci[1, i]:+.3f}]  "
              f"bootstrap-significant: {sig};  typical |kink| at shifted positions: {ctrl[i]:.3f}")

    # plot: binned residuals vs x for each model
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    bin_edges = np.logspace(np.log10(x.min()), np.log10(x.max()), 22)
    styles = {"Newton": ("0.5", ":"), "Square-Rule (original)": ("tab:green", "-."),
              "RAR (McGaugh 2016)": ("tab:red", "--"), "MOND (simple mu)": ("tab:orange", ":"),
              "Three-phase": ("tab:blue", "-"), "Three-phase hybrid": ("tab:purple", "-")}
    for name, (color, ls) in styles.items():
        c, m, _ = st.binned_median(x, st.log_residuals(df, models[name]), bin_edges)
        ax.plot(c * st.SIGMA_STAR, m, color=color, ls=ls, lw=2, marker="o", ms=3, label=name)
    ax.axhline(0, color="k", lw=0.8)
    for lo, hi, col, txt in [(None, st.SIGMA_LOW, "tab:blue", "gas"),
                             (st.SIGMA_LOW, st.SIGMA_STAR, "tab:cyan", "liquid"),
                             (st.SIGMA_STAR, None, "tab:purple", "solid")]:
        ax.axvspan(lo or 1e-3, hi or 1e5, color=col, alpha=0.07)
    ax.axvline(st.SIGMA_LOW, color="tab:blue", ls="--", lw=1)
    ax.axvline(st.SIGMA_STAR, color="tab:purple", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(x.min() * st.SIGMA_STAR, x.max() * st.SIGMA_STAR)
    ax.set_ylim(-0.3, 0.8)
    ax.set_xlabel("effective surface density Σ_eff = g_bar / 2πG  [M☉/pc²]")
    ax.set_ylabel("median log(g_obs / g_model)  [dex]")
    ax.set_title("Where do the models miss the data? (0 = perfect, >0 = model too weak)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = st.RESULTS_DIR / "model_residuals.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"\nSaved results to {st.RESULTS_DIR.relative_to(st.REPO_ROOT)}/")


if __name__ == "__main__":
    main()
