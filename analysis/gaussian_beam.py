"""HeNe laser: beam waist from two beam-radius methods, and resonator stability.

Data: data/hene_knife_edge.csv (radii from erf fits to the knife-edge scans),
data/hene_clip_level.csv (90 % / 10 % clip positions), data/hene_stability.csv.
"""

import numpy as np
import pandas as pd

from _common import DATA, plt, save, write_results
from optics import fitting as ft
from optics import gaussian_beam as gb

LAM = gb.HENE_WAVELENGTH
CLIP_POS_ERR = 20e-6                     # blade position reading, per position
CLIP_W_ERR = np.sqrt(2) * CLIP_POS_ERR / gb.clip_level_factor()

ke = pd.read_csv(DATA / "hene_knife_edge.csv")
cl = pd.read_csv(DATA / "hene_clip_level.csv")
cl["D_c_um"] = cl.x10_um - cl.x90_um
cl["w_um"] = gb.clip_width_to_radius(cl.D_c_um)


def waist_fit(z, w, w_err=None):
    model = lambda z, w0: gb.beam_radius(z, w0, LAM)
    return ft.fit(model, z, w, p0=[2e-4], sigma=w_err)


def dataset(method, res):
    if method == "knife edge":
        d = ke[ke.resonator == res]
        return d.z_cm.values * 1e-2, d.w_um.values * 1e-6, d.w_err_um.values * 1e-6
    d = cl[cl.resonator == res]
    return d.z_cm.values * 1e-2, d.w_um.values * 1e-6, np.full(len(d), CLIP_W_ERR)


# --- waist fits -------------------------------------------------------------
rows = []
fits = {}
for res in (1, 2):
    for method in ("knife edge", "clip level"):
        z, w, s = dataset(method, res)
        weighted = waist_fit(z, w, s)
        unweighted = waist_fit(z, w)
        w0, sw0 = weighted.params[0], weighted.errors[0]
        # quoted errors only cover reading precision; if the points scatter
        # more than that (chi2/dof > 1), scale the error up (Birge ratio)
        birge = max(1.0, np.sqrt(weighted.reduced_chi2))
        fits[(res, method)] = (w0, sw0 * birge)
        rows.append(dict(
            resonator=res, method=method,
            w0_weighted=w0 * 1e6, w0_err_stat=sw0 * 1e6, chi2_dof=weighted.reduced_chi2,
            w0_err_scaled=sw0 * birge * 1e6,
            w0_unweighted=unweighted.params[0] * 1e6, w0_err_unweighted=unweighted.errors[0] * 1e6,
            zR_cm=gb.rayleigh_length(w0, LAM) * 100,
            theta_mrad=gb.divergence(w0, LAM) * 1e3,
        ))
table = pd.DataFrame(rows)

# --- figure: beam radius vs distance -----------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
z_plot = np.linspace(0, 2.3, 300)
for ax, res in zip(axes, (1, 2)):
    for method, marker, color in (("knife edge", "o", "C0"), ("clip level", "s", "C1")):
        z, w, s = dataset(method, res)
        w0, _ = fits[(res, method)]
        ax.errorbar(z * 100, w * 1e6, s * 1e6, fmt=marker, color=color, ms=5, capsize=3,
                    label=f"{method}: w0 = {w0*1e6:.1f} µm")
        ax.plot(z_plot * 100, gb.beam_radius(z_plot, w0, LAM) * 1e6, color=color, lw=1)
    w0_ke = fits[(res, "knife edge")][0]
    zr = gb.rayleigh_length(w0_ke, LAM)
    ax.axvspan(0, zr * 100, color="grey", alpha=0.15, label="Rayleigh range (knife edge)")
    ax.set_title(f"Resonator {res}" + ("" if res == 1 else "  (mirror moved +10 cm)"))
    ax.set_xlabel("distance from plane mirror z (cm)")
    ax.legend(loc="upper left")
axes[0].set_ylabel("1/e² beam radius w (µm)")
save(fig, "beam_radius_vs_z.png")

# --- figure: why the far field only measures the divergence ------------------
# Allowing M^2 != 1 makes w0 and M^2 degenerate when all points are far field:
# only the product theta = M^2 lambda / (pi w0) is constrained.
z, w, s = dataset("knife edge", 1)
w0_grid = np.linspace(120e-6, 330e-6, 211)
m2_grid = np.linspace(0.3, 2.0, 1701)
chi2_fixed = np.array([np.sum(((gb.beam_radius(z, w0, LAM) - w) / s) ** 2) for w0 in w0_grid])
chi2_profile = np.array([
    min(np.sum(((gb.beam_radius(z, w0, LAM, m2) - w) / s) ** 2) for m2 in m2_grid)
    for w0 in w0_grid
])
best_m2 = np.array([
    m2_grid[np.argmin([np.sum(((gb.beam_radius(z, w0, LAM, m2) - w) / s) ** 2)
                       for m2 in m2_grid])]
    for w0 in w0_grid
])
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
ax = axes[0]
ax.plot(w0_grid * 1e6, chi2_fixed - chi2_fixed.min(), label="M² fixed to 1")
ax.plot(w0_grid * 1e6, chi2_profile - chi2_fixed.min(), label="M² free (minimised over M²)")
ax.set_ylim(-1, 30)
ax.set_xlabel("waist w0 (µm)")
ax.set_ylabel("Δχ²")
ax.set_title("How well is w0 determined?")
ax.legend()
ax = axes[1]
ax.plot(w0_grid * 1e6, best_m2, color="C1")
ax.axhline(1.0, color="grey", lw=0.8)
ax.set_xlabel("waist w0 (µm)")
ax.set_ylabel("best-fit M²")
ax.set_title(r"Free M² just rescales: $\theta = M^2\lambda/(\pi w_0)$ fixed")
fig.suptitle("Knife-edge data, resonator 1: all points are in the far field (z > 4 zR)")
save(fig, "waist_m2_degeneracy.png")

# --- resonator stability -----------------------------------------------------
st = pd.read_csv(DATA / "hene_stability.csv")
U_BG, R2_NOMINAL = 0.25, 70.0
st["U_net"] = st.U_V - U_BG
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(st.L_cm, st.U_net, "o-", ms=5)
ax.axvline(R2_NOMINAL, color="C3", ls="--", lw=1, label="nominal limit L = R2 = 70 cm")
ax.axvline(st.L_cm.max(), color="k", ls=":", lw=1,
           label=f"last lasing point L = {st.L_cm.max():.1f} cm")
ax.axvspan(60, R2_NOMINAL, color="C2", alpha=0.08)
ax.text(64, 0.25, "stable if\nR2 = 70 cm", fontsize=9, color="C2")
ax.set_xlim(63, 78)
ax.set_xlabel("mirror separation L (cm)")
ax.set_ylabel("photodiode signal above background (V)")
ax.set_title("Output power near the stability limit")
ax.legend(loc="upper right")
save(fig, "resonator_stability.png")

# --- results -----------------------------------------------------------------
fmt = lambda r: (f"| {r.resonator} | {r.method} | {r.w0_weighted:.1f} ± {r.w0_err_stat:.1f} | "
                 f"{r.chi2_dof:.1f} | {r.w0_weighted:.1f} ± {r.w0_err_scaled:.1f} | "
                 f"{r.w0_unweighted:.1f} ± {r.w0_err_unweighted:.1f} | {r.zR_cm:.1f} | {r.theta_mrad:.3f} |")
lines = [
    "# Gaussian beam results (HeNe)",
    "",
    f"Clip-level factor D_c / w = {gb.clip_level_factor():.4f} (eps = 0.1); "
    f"clip-level radius error {CLIP_W_ERR*1e6:.1f} µm.",
    "",
    "| Resonator | Method | w0 weighted (µm) | χ²/dof | w0, error scaled by √(χ²/dof) (µm) "
    "| w0 unweighted (µm) | zR (cm) | θ (mrad) |",
    "|---|---|---|---|---|---|---|---|",
    *[fmt(r) for r in table.itertuples()],
    "",
    f"Stability: laser still oscillates at L = {st.L_cm.max():.1f} cm > nominal R2 = 70 cm, "
    f"so the real R2 is at least {st.L_cm.max():.1f} cm.",
]
write_results("gaussian_beam.md", lines)
