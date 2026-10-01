"""
dsi_test.py – pre-registered DSI screen test on the SPARC RAR residuals.

Implements the protocol of Meta-Ledger Scaffold Part I, Appendix B
(DOI 10.5281/zenodo.19413799) on the SPARC data:

  observable    RAR residual log10(g_obs / g_RAR) in dex
  carriers      primary (pre-registered): xi = ln(g_bar / g_dagger)
                further channels: ln(Sigma_stars / Sigma*), ln(R / kpc)
  detrending    (declared) per-galaxy mean removed, then a global quadratic
                in the carrier coordinate
  frequency     fixed a priori: omega_phi = 2*pi/ln(phi) = 13.057, plus 2*omega_phi
  nulls         (1) surrogates: each galaxy gets a random phase (keeps the
                    structure inside a galaxy, destroys a common phase)
                (2) look-elsewhere: same amplitude at other frequencies
  coherence     Rayleigh test of the per-galaxy phases (galaxies spanning
                >= 2 phi-steps)
  control       deliberately wrong carrier: the linear axis g_bar / g_dagger,
                tested with the same number of cycles across its range

Two readings of the phase are tested per channel:
  common phase  (Part I, Appendix B) detection only if the amplitude is
                significant against the surrogates AND the phases are
                coherent across galaxies
  domain-local  (Part IV, Step 20 C4) every domain (galaxy) has its own
                anchor/phase, only the step ln(phi) is fixed: the power summed
                incoherently over galaxies at omega_phi must exceed the power
                at the other frequencies between omega_phi/2 and 2 omega_phi.
                The sensitivity is calibrated by injecting modulations with a
                random phase per galaxy.
Otherwise label H. Phases of different carriers are not compared: that needs
carrier offsets declared in advance (Mapping level), which the theory does not
fix yet.

Usage:
    python sparc/dsi_test.py
Output:
    results/sparc_analysis/dsi_test.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sparc_tools as st

N_SURROGATES = 5000
ALPHA = 0.05


def detrend(r, x, galaxies):
    """Declared detrending: per-galaxy mean removed, then a global quadratic in x."""
    r = r - (pd.Series(r).groupby(galaxies).transform("mean")).to_numpy()
    return r - np.polyval(np.polyfit(x, r, 2), x)


def galaxy_sums(r, x, omega, groups):
    """Complex amplitude sum_points r * exp(-i omega x) for every galaxy."""
    e = r * np.exp(-1j * omega * x)
    return np.array([e[idx].sum() for idx in groups])


def amplitude(z, n_points):
    return 2 * abs(z.sum()) / n_points


def surrogate_null(z, n_points, rng):
    phases = np.exp(1j * rng.uniform(0, 2 * np.pi, (N_SURROGATES, len(z))))
    return 2 * np.abs((z * phases).sum(axis=1)) / n_points


def domain_power(r, x, starts, omegas):
    """Phase-free power: sum over galaxies of |sum r exp(-i omega x)|^2.
    Data must be sorted by galaxy; starts are the first index of each galaxy."""
    e = r[:, None] * np.exp(-1j * np.outer(x, omegas))
    return (np.abs(np.add.reduceat(e, starts, axis=0)) ** 2).sum(axis=0)


def domain_local_test(r, x, galaxies, rng, omega, step, n_inject=20):
    """Part IV, Step 20 C4: DSI with a free phase per galaxy.

    Returns (galaxies used, power ratio at omega vs. median of the band,
    fraction of band frequencies with at least the same power, amplitude in dex
    detected in >= 90 % of injections)."""
    order = np.argsort(galaxies, kind="stable")
    r, x, galaxies = r[order], x[order], galaxies[order]
    starts = np.flatnonzero(np.r_[True, galaxies[1:] != galaxies[:-1]])
    sizes = np.diff(np.r_[starts, len(x)])
    spans = np.maximum.reduceat(x, starts) - np.minimum.reduceat(x, starts)
    keep_gal = spans / step >= 2
    keep = np.repeat(keep_gal, sizes)
    r, x, galaxies = r[keep], x[keep], galaxies[keep]
    starts = np.flatnonzero(np.r_[True, galaxies[1:] != galaxies[:-1]])
    gal_index = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(x)]))

    resolution = 2 * np.pi / np.median(spans[keep_gal])
    omegas = np.linspace(omega / 2, 2 * omega, 300)
    band = omegas[(np.abs(omegas - omega) > resolution) & (np.abs(omegas - 2 * omega) > resolution)]
    test_omegas = np.r_[omega, band]

    def evaluate(rr):
        power = domain_power(rr, x, starts, test_omegas)
        return (power[1:] >= power[0]).mean(), power[0] / np.median(power[1:])

    p, ratio = evaluate(r)
    a90 = np.nan
    for a in (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08):
        hits = sum(evaluate(r + a * np.cos(omega * x + rng.uniform(0, 2 * np.pi, len(starts))[gal_index]))[0]
                   < ALPHA for _ in range(n_inject))
        if hits >= 0.9 * n_inject:
            a90 = a
            break
    return len(starts), ratio, p, a90


def test_channel(name, r, x, galaxies, rng, omega=st.OMEGA_PHI, step=np.log(st.PHI)):
    """Runs the Appendix-B test on one carrier; returns spectrum and verdict."""
    r = detrend(r, x, galaxies)
    groups = list(pd.Series(np.arange(len(r))).groupby(galaxies).indices.values())
    n = len(r)
    span = np.array([np.ptp(x[idx]) / step for idx in groups])
    print(f"\n== {name}: {n} points, {np.ptp(x) / step:.1f} steps in total, "
          f"median {np.median(span):.1f} per galaxy")
    detected = False
    for label, w in [("fundamental", omega), ("harmonic 2x", 2 * omega)]:
        z = galaxy_sums(r, x, w, groups)
        a = amplitude(z, n)
        null = surrogate_null(z, n, rng)
        p_sur = (null >= a).mean()
        sel = span >= 2
        rayleigh_r = abs(np.exp(1j * np.angle(z[sel])).mean())
        p_ray = np.exp(-sel.sum() * rayleigh_r ** 2)
        hit = p_sur < ALPHA and p_ray < ALPHA
        if label == "fundamental":
            detected = hit
        print(f"   {label}: amplitude {a:.4f} dex ({100 * (10 ** a - 1):.1f} %), "
              f"surrogates p = {p_sur:.3f} (95 % threshold {np.quantile(null, 0.95):.4f} dex), "
              f"Rayleigh over {sel.sum()} galaxies p = {p_ray:.3f}")
    omegas = np.linspace(3, 40, 400) * omega / st.OMEGA_PHI
    spectrum = np.array([amplitude(galaxy_sums(r, x, w, groups), n) for w in omegas])
    p_le = (spectrum[np.abs(omegas - omega) > omega / 13] >= amplitude(galaxy_sums(r, x, omega, groups), n)).mean()
    print(f"   look-elsewhere: {p_le:.2f} of other frequencies reach the same amplitude"
          f"  ->  common phase: {'DETECTION' if detected else 'label H'}")
    n_gal, ratio, p_dom, a90 = domain_local_test(r, x, galaxies, rng, omega, step)
    detected_dom = p_dom < ALPHA
    print(f"   domain-local phases ({n_gal} galaxies): power at omega / band median = {ratio:.2f}, "
          f"{p_dom:.2f} of band frequencies reach it; sensitivity (90 % of injections): "
          f"{a90:.3f} dex ({100 * (10 ** a90 - 1):.1f} %)  ->  "
          f"{'DETECTION' if detected_dom else 'label H'}")
    return omegas / omega * st.OMEGA_PHI, spectrum, detected or detected_dom


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    df = st.load_sparc()
    r_all = st.log_residuals(df, st.rar(df))
    gal = df["galaxy"].to_numpy()
    x_gbar = (df["g_bar"] / st.g_dagger(st.SIGMA_STAR)).to_numpy()
    stars = (df["sigma_stars"] > 0).to_numpy()

    channels = {
        "ln(g_bar / g_dagger)  [primary]": (r_all, np.log(x_gbar), gal),
        "ln(Sigma_stars / Sigma*)": (r_all[stars], np.log(df["sigma_stars"].to_numpy()[stars] / st.SIGMA_STAR), gal[stars]),
        "ln(R / kpc)": (r_all, np.log(df["R"].to_numpy()), gal),
    }
    results = {name: test_channel(name, *args, rng) for name, args in channels.items()}

    # wrong carrier: linear axis with as many cycles over its range as the primary carrier
    xi = np.log(x_gbar)
    omega_lin = st.OMEGA_PHI * np.ptp(xi) / np.ptp(x_gbar)
    print("\nControl, wrong carrier (linear axis g_bar/g_dagger):", end="")
    _, _, wrong = test_channel("linear g_bar / g_dagger  [control]", r_all, x_gbar, gal, rng,
                               omega=omega_lin, step=2 * np.pi / omega_lin)

    primary = results["ln(g_bar / g_dagger)  [primary]"][2]
    verdict = "DETECTION" if primary else "no detection -> label H (horizon)"
    print(f"\nVerdict on the pre-registered carrier (common or domain-local phase): {verdict}")
    if wrong:
        print("WARNING: the wrong carrier also shows a signal -> any detection is inadmissible")

    fig, ax = plt.subplots(figsize=(9, 4.8))
    styles = ["-", "--", ":"]
    for (name, (om, spec, _)), ls in zip(results.items(), styles):
        ax.plot(om, spec, ls=ls, lw=1.3, label=name)
    ax.axvline(st.OMEGA_PHI, color="tab:red", lw=1.5, alpha=0.6, label="2π/ln φ = 13.057 (pre-registered)")
    ax.axvline(2 * st.OMEGA_PHI, color="tab:orange", ls=":", alpha=0.8, label="harmonic 26.11")
    ax.set_xlabel("log-frequency ω  (per unit of the log carrier)")
    ax.set_ylabel("amplitude of RAR residuals  [dex]")
    ax.set_title(f"DSI screen test on SPARC (Appendix B): {verdict}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = st.RESULTS_DIR / "dsi_test.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out.relative_to(st.REPO_ROOT)}")


if __name__ == "__main__":
    main()
