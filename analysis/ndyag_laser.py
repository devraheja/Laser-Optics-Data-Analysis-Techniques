"""Diode-pumped Nd:YAG laser: pump-diode calibration, laser threshold and slope efficiency.

Data: data/yag_pump_diode.csv (diode current vs optical power),
data/yag_laser_output.csv (diode current vs 1064 nm output behind the OC100).
No lasing was observed for I <= 275 mA.
"""

import numpy as np
import pandas as pd

from _common import DATA, plt, save, write_results
from optics import fitting as ft

rng = np.random.default_rng(1064)
pump = pd.read_csv(DATA / "yag_pump_diode.csv")
laser = pd.read_csv(DATA / "yag_laser_output.csv")
I_NO_LASING = 275.0  # mA, highest current with zero output

# --- 1. pump diode -------------------------------------------------------------
above = pump[pump.I_mA > 200]
line_fit = ft.fit(ft.line, above.I_mA, above.P_mW, [1.0, -200.0])
s, b = line_fit.params
I_th, sI_full, sI_nocov = ft.propagate_ratio(s, b, line_fit.cov)
hinge_fit = ft.fit(ft.hinge, above.I_mA, above.P_mW, [1.0, 200.0])  # fits I_th directly
boot = ft.residual_bootstrap(ft.hinge, above.I_mA.values, above.P_mW.values,
                             [1.0, 200.0], rng=rng)
I_boot_lo, I_boot_hi = np.percentile(boot[:, 1], [16, 84])

# sensitivity to where the "linear" region is taken to start
starts = [225, 250, 300, 350, 400]
window = []
for st in starts:
    d = pump[pump.I_mA >= st]
    f = ft.fit(ft.hinge, d.I_mA, d.P_mW, [1.0, 200.0])
    window.append((st, len(d), *f.params, *f.errors))

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
I_plot = np.linspace(190, 580, 100)
ax = axes[0]
ax.plot(pump.I_mA, pump.P_mW, "o", ms=5, label="measured")
ax.plot(I_plot, ft.line(I_plot, s, b), lw=1,
        label=f"fit I > 200 mA: s = {s:.3f} W/A")
ax.axvline(I_th, color="k", ls=":", lw=1, label=f"I_th = {I_th:.1f} mA")
ax.set_xlabel("diode current I (mA)")
ax.set_ylabel("pump power P (mW)")
ax.set_title("Pump diode characteristic")
ax.legend()
ax = axes[1]
ax.semilogy(pump.I_mA[1:], pump.P_mW[1:], "o", ms=5)
ax.axvline(I_th, color="k", ls=":", lw=1)
ax.text(40, 3, "spontaneous\nemission", fontsize=9)
ax.text(300, 3, "stimulated\nemission", fontsize=9)
ax.set_xlabel("diode current I (mA)")
ax.set_ylabel("pump power P (mW), log scale")
ax.set_title("Below threshold: factor ~8 jump at 200→225 mA")
ax = axes[2]
st_, n_, _, ith_, _, sith_ = zip(*window)
ax.errorbar(st_, ith_, sith_, fmt="o", capsize=3)
for x, y, n in zip(st_, ith_, n_):
    ax.annotate(f"{n} pts", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
ax.set_xlabel("lowest current included in the fit (mA)")
ax.set_ylabel("fitted threshold I_th (mA)")
ax.set_title("Sensitivity to the fit window")
save(fig, "ndyag_pump_diode.png")

# --- 2. laser characteristic ---------------------------------------------------
# pump power from the diode calibration; the calibration uncertainty is
# propagated by resampling (s, b) from their joint distribution
P_pump = s * laser.I_mA.values + b
P_out = laser.P_out_mW.values
lin = ft.fit(ft.hinge, P_pump, P_out, [0.1, 100.0])   # P_out = eta (P_pump - P_th)
eta, P_th = lin.params

draws = rng.multivariate_normal(line_fit.params, line_fit.cov, size=4000)
mc = np.array([ft.fit(ft.hinge, sd * laser.I_mA.values + bd, P_out, [0.1, 100.0]).params
               for sd, bd in draws])
sP_th_calib = mc[:, 1].std()
sP_th_total = np.hypot(lin.errors[1], sP_th_calib)
s_eta_total = np.hypot(lin.errors[0], mc[:, 0].std())

# The linear fit puts the threshold above the lowest lasing point, which is
# impossible. A quadratic captures the curvature visible in the residuals.
quad = np.polyfit(P_pump, P_out, 2)
roots = np.roots(quad)
P_th_quad = float(roots[(roots > 0) & np.isreal(roots)].real.min())
P_no_lasing = s * I_NO_LASING + b
P_lowest_lasing = P_pump[0]

fig, axes = plt.subplots(2, 1, figsize=(6.5, 6), sharex=True,
                         gridspec_kw={"height_ratios": [3, 1]})
x = np.linspace(60, 410, 200)
ax = axes[0]
ax.plot(P_pump, P_out, "o", ms=5, label="measured")
ax.plot(x, np.clip(ft.hinge(x, eta, P_th), 0, None), lw=1,
        label=f"linear: η = {eta:.3f}, P_th = {P_th:.0f} mW")
ax.plot(x, np.clip(np.polyval(quad, x), 0, None), lw=1, ls="--",
        label=f"quadratic: P_th = {P_th_quad:.0f} mW")
ax.axvspan(P_no_lasing, P_lowest_lasing, color="C2", alpha=0.15,
           label="threshold must lie here (data)")
ax.set_ylabel("laser output P_out (mW)")
ax.set_title("Nd:YAG laser characteristic (OC100)")
ax.legend(loc="upper left")
ax = axes[1]
ax.axhline(0, color="grey", lw=0.8)
ax.plot(P_pump, lin.residuals, "o-", ms=4, label="linear fit residuals")
ax.plot(P_pump, P_out - np.polyval(quad, P_pump), "s--", ms=4, label="quadratic")
ax.set_xlabel("pump power P_pump (mW)")
ax.set_ylabel("residual (mW)")
ax.legend(fontsize=8)
save(fig, "ndyag_laser_characteristic.png")

# --- results -------------------------------------------------------------------
lines = [
    "# Nd:YAG results",
    "",
    "## Pump diode (fit to I > 200 mA)",
    f"- slope s = {s:.4f} ± {line_fit.errors[0]:.4f} W/A, intercept b = {b:.2f} ± {line_fit.errors[1]:.2f} mW",
    f"- threshold I_th = -b/s = {I_th:.1f} mA",
    f"  - error ignoring the s-b correlation: ± {sI_nocov:.2f} mA",
    f"  - error with full covariance: ± {sI_full:.2f} mA "
    f"(fitting I_th directly gives ± {hinge_fit.errors[1]:.2f} mA)",
    f"  - residual bootstrap 68% interval: [{I_boot_lo:.1f}, {I_boot_hi:.1f}] mA",
    "- fit-window sensitivity:",
    "",
    "| lowest current (mA) | points | I_th (mA) |",
    "|---|---|---|",
    *[f"| {w[0]} | {w[1]} | {w[3]:.1f} ± {w[5]:.1f} |" for w in window],
    "",
    "## Laser characteristic",
    f"- linear fit: slope efficiency η = {eta:.3f} ± {s_eta_total:.3f}, "
    f"P_th = {P_th:.1f} ± {sP_th_total:.1f} mW "
    f"(fit {lin.errors[1]:.1f}, pump calibration {sP_th_calib:.1f})",
    f"- quadratic fit: P_th = {P_th_quad:.1f} mW",
    f"- data bracket: no lasing at {P_no_lasing:.1f} mW, lasing at {P_lowest_lasing:.1f} mW",
    f"- residual signs of the linear fit: "
    + " ".join("+" if r > 0 else "-" for r in lin.residuals),
]
write_results("ndyag_laser.md", lines)
