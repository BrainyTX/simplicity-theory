"""
dsi_test.py – pre-registered DSI screen test on the SPARC RAR residuals.

Implements the protocol of Meta-Ledger Scaffold Part I, Appendix B
(DOI 10.5281/zenodo.19413799) on the SPARC data:

  carrier       mu = g_bar / g_dagger (dimensionless), xi = ln(mu)
  observable    RAR residual log10(g_obs / g_RAR) in dex
  detrending    (declared) per-galaxy mean removed, then a global quadratic in xi
  frequency     fixed a priori: omega_phi = 2*pi/ln(phi) = 13.057, plus 2*omega_phi
  nulls         (1) surrogates: each galaxy gets a random phase (keeps the
                    structure inside a galaxy, destroys a common phase)
                (2) look-elsewhere: same amplitude at other frequencies
  coherence     Rayleigh test of the per-galaxy phases (galaxies spanning
                >= 2 phi-steps) and phase difference of two random halves
  control       deliberately wrong carrier xi = ln(R)

Verdict: detection only if the amplitude is significant against the
surrogates AND the phases are coherent across galaxies; otherwise label H.

Usage:
    python sparc/dsi_test.py
Output:
    results/sparc_analysis/dsi_test.png
"""
import matplotlib.pyplot as plt
import numpy as np

import sparc_tools as st

N_SURROGATES = 5000
ALPHA = 0.05


def detrended_residuals(df, xi):
    r = st.log_residuals(df, st.rar(df))
    r = r - df.assign(r=r).groupby("galaxy")["r"].transform("mean").to_numpy()
    return r - np.polyval(np.polyfit(xi, r, 2), xi)


def galaxy_sums(r, x, omega, groups):
    """Complex amplitude sum_points r * exp(-i omega x) for every galaxy."""
    e = r * np.exp(-1j * omega * x)
    return np.array([e[idx].sum() for idx in groups])


def amplitude(z, n_points):
    return 2 * abs(z.sum()) / n_points


def surrogate_null(z, n_points, rng):
    phases = np.exp(1j * rng.uniform(0, 2 * np.pi, (N_SURROGATES, len(z))))
    return 2 * np.abs((z * phases).sum(axis=1)) / n_points


def main():
    st.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    df = st.load_sparc()
    xi = np.log(df["g_bar"] / st.g_dagger(st.SIGMA_STAR)).to_numpy()
    r = detrended_residuals(df, xi)
    groups = list(df.groupby("galaxy").indices.values())
    n = len(r)
    ln_phi = np.log(st.PHI)
    span = np.array([np.ptp(xi[idx]) / ln_phi for idx in groups])
    print(f"Carrier xi = ln(g_bar/g_dagger): {np.ptp(xi) / ln_phi:.1f} phi-steps in total, "
          f"median {np.median(span):.1f} per galaxy; residual scatter {np.std(r):.3f} dex\n")

    omegas = np.linspace(3, 40, 600)
    spectrum = np.array([amplitude(galaxy_sums(r, xi, w, groups), n) for w in omegas])

    detected = {}
    for label, w in [("fundamental 2pi/ln(phi)", st.OMEGA_PHI), ("harmonic 2x", 2 * st.OMEGA_PHI)]:
        z = galaxy_sums(r, xi, w, groups)
        a = amplitude(z, n)
        null = surrogate_null(z, n, rng)
        p_sur = (null >= a).mean()
        p_le = (spectrum[np.abs(omegas - w) > 1.0] >= a).mean()
        sel = span >= 2
        rayleigh_r = abs(np.exp(1j * np.angle(z[sel])).mean())
        p_ray = np.exp(-sel.sum() * rayleigh_r ** 2)
        idx = rng.permutation(len(z))
        dphi = np.degrees(np.angle(z[idx[: len(z) // 2]].sum() / z[idx[len(z) // 2:]].sum()))
        detected[label] = p_sur < ALPHA and p_ray < ALPHA
        print(f"{label}: amplitude {a:.4f} dex ({100 * (10 ** a - 1):.1f} %), "
              f"phase {np.degrees(np.angle(z.sum())):+.0f} deg")
        print(f"   surrogates p = {p_sur:.3f} (95 % detection threshold {np.quantile(null, 0.95):.4f} dex), "
              f"look-elsewhere {p_le:.2f}")
        print(f"   phase coherence: Rayleigh over {sel.sum()} galaxies R = {rayleigh_r:.3f}, p = {p_ray:.3f}; "
              f"two random halves differ by {dphi:+.0f} deg")

    z_wrong = galaxy_sums(r, np.log(df["R"].to_numpy()), st.OMEGA_PHI, groups)
    a_wrong = amplitude(z_wrong, n)
    p_wrong = (surrogate_null(z_wrong, n, rng) >= a_wrong).mean()
    print(f"\nControl, wrong carrier ln(R): amplitude {a_wrong:.4f} dex, p = {p_wrong:.3f}")
    verdict = "DETECTION" if detected["fundamental 2pi/ln(phi)"] else "no detection -> label H (horizon)"
    print(f"\nVerdict (Appendix B, all-or-nothing): {verdict}")

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(omegas, spectrum, color="0.3", lw=1.2, label="amplitude of RAR residuals")
    ax.axvline(st.OMEGA_PHI, color="tab:red", ls="--", label="2π/ln φ = 13.057 (pre-registered)")
    ax.axvline(2 * st.OMEGA_PHI, color="tab:orange", ls=":", label="harmonic 26.11")
    ax.set_xlabel("log-frequency ω  (per unit of ln(g_bar/g†))")
    ax.set_ylabel("amplitude  [dex]")
    ax.set_title(f"DSI screen test on SPARC (Appendix B): {verdict}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = st.RESULTS_DIR / "dsi_test.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out.relative_to(st.REPO_ROOT)}")


if __name__ == "__main__":
    main()
