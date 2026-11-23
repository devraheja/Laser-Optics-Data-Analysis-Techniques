"""Michelson interferometer with a tilted mirror as a Fourier-transform spectrometer.

Real data (from the lab report):
  data/michelson_calibration.csv  FFT peak index of two interference filters
  data/michelson_slit_peaks.csv   peak index of eight prism-selected halogen bands
  data/michelson_double_slit.csv  peak pairs for five double-slit configurations
  measured widths of the white-light spectrum (26 bins) and fringe envelope (41 px)

The raw camera images were not kept, so the full processing chain
(interferogram -> FFT -> spectrum, Hilbert envelope) is demonstrated on
interferograms simulated with the measured parameters.
"""

import numpy as np
import pandas as pd

from _common import DATA, plt, save, write_results
from optics import interferometry as fts

C = fts.C
N_PIX = 1280
PEAK_ERR = 0.5           # bins: peak position read to the nearest FFT bin
SPECTRAL_FWHM_BINS = 26  # white-light spectrum, measured
ENVELOPE_FWHM_PX = 41    # white-light fringe envelope, measured
BARYCENTER_NM = 735      # white-light barycenter, measured
LINE_FWHM_BINS = 4.5     # width of a single prism-selected band (4-5 bins measured)

cal = pd.read_csv(DATA / "michelson_calibration.csv")
bands = pd.read_csv(DATA / "michelson_slit_peaks.csv")
dslit = pd.read_csv(DATA / "michelson_double_slit.csv")
rng = np.random.default_rng(632)

# --- 1. calibration ----------------------------------------------------------
# Fringe count across the frame is proportional to optical frequency, so the
# physical model is nu = a k with no intercept. A straight line with an
# intercept through two points has zero degrees of freedom and an
# intercept uncertainty of tens of THz, which is what inflated the original
# wavelength errors.
nu_cal = C / (cal.wavelength_nm.values * 1e-9)
k_cal = cal.peak_index.values

m_line, c_line = np.polyfit(k_cal, nu_cal, 1)

a_i = nu_cal / k_cal
sa_i = a_i * PEAK_ERR / k_cal
w = 1.0 / sa_i**2
a = np.sum(w * a_i) / np.sum(w)
sa = 1.0 / np.sqrt(np.sum(w))
chi2_cal = float(np.sum(w * (a_i - a) ** 2))      # 1 degree of freedom

opd_per_px = C / (a * N_PIX)                       # optical path difference per pixel

k = bands.peak_index.values.astype(float)
lam_prop = C / (a * k) * 1e9
slam_prop = lam_prop * np.hypot(sa / a, PEAK_ERR / k)
lam_line = C / (m_line * k + c_line) * 1e9
# report errors (Table 1), quoted for comparison
REPORT_ERR = {"violet": 41.3, "blue": 47.3, "green": 54.4, "yellow": 62.9,
              "orange": 67.0, "red": 75.6, "IR2": 112.7, "IR1": 118.6}
report_err = np.array([REPORT_ERR[b] for b in bands.band])

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
kk = np.linspace(0, 150, 100)
ax.plot(kk, a * kk / 1e12, label=f"ν = a·k,  a = {a/1e12:.3f} ± {sa/1e12:.3f} THz/bin")
ax.plot(kk, (m_line * kk + c_line) / 1e12, "--",
        label=f"two-point line, intercept {c_line/1e12:.0f} THz")
ax.errorbar(k_cal, nu_cal / 1e12, xerr=PEAK_ERR, fmt="o", color="k", capsize=3,
            label="green / red filters")
ax.set_xlim(0, 150)
ax.set_ylim(0, 800)
ax.set_xlabel("FFT peak index k (fringes per frame)")
ax.set_ylabel("optical frequency ν (THz)")
ax.set_title("Calibration")
ax.legend(loc="upper left")

ax = axes[1]
order = np.argsort(lam_prop)
y = np.arange(len(k))
ax.errorbar(lam_line[order], y + 0.15, xerr=report_err[order], fmt="s", color="C1",
            capsize=3, label="line with intercept (report errors)")
ax.errorbar(lam_prop[order], y - 0.15, xerr=slam_prop[order], fmt="o", color="C0",
            capsize=3, label="proportional calibration")
ax.set_yticks(y)
ax.set_yticklabels(bands.band.values[order])
ax.set_xlabel("reconstructed wavelength (nm)")
ax.set_title("Prism-selected halogen bands")
ax.legend(loc="lower right", fontsize=8)
save(fig, "fts_calibration.png")

# --- 2. white light: simulated interferogram with the measured spectrum ------
k_axis = np.arange(1, N_PIX // 2)
k0 = C / (BARYCENTER_NM * 1e-9) / a
S_true = fts.gaussian_band(k_axis, k0, SPECTRAL_FWHM_BINS)
row = fts.interferogram(S_true, k_axis, N_PIX, zpd=640, noise=0.01, rng=rng)
k_fft, S_fft = fts.spectrum_from_interferogram(row)
env = fts.envelope(row)
p = np.arange(N_PIX)

sel = (k_fft > 40) & (k_fft < 140)
sim_spec_fwhm = fts.fwhm(k_fft[sel], S_fft[sel])
sim_env_fwhm = fts.fwhm(p, env)
predicted_env = fts.gaussian_envelope_width(SPECTRAL_FWHM_BINS, N_PIX)

fig, axes = plt.subplots(2, 1, figsize=(8, 6))
ax = axes[0]
ax.plot(p, row, lw=0.6, label="simulated camera row")
ax.plot(p, env + row.mean(), "k--", lw=1, label="Hilbert envelope")
ax.set_xlim(500, 780)
ax.set_xlabel("pixel")
ax.set_ylabel("intensity (norm.)")
ax.set_title(f"White light: envelope FWHM {sim_env_fwhm:.1f} px "
             f"(Gaussian prediction {predicted_env:.1f} px, measured {ENVELOPE_FWHM_PX} px)")
ax.legend(loc="upper right")
ax = axes[1]
ax.plot(k_fft, S_fft / S_fft[sel].max(), label="|FFT| of the row")
ax.plot(k_axis, S_true, "k--", lw=1, label="input spectrum")
ax.set_xlim(0, 200)
ax.set_xlabel("FFT index k")
ax.set_ylabel("spectral density (norm.)")
sec = ax.secondary_xaxis("top", functions=(lambda x: x * a / 1e12, lambda f: f * 1e12 / a))
sec.set_xlabel("ν (THz)")
ax.set_title(f"Recovered FWHM {sim_spec_fwhm:.1f} bins (input {SPECTRAL_FWHM_BINS})", pad=28)
ax.legend()
save(fig, "fts_white_light.png")

# --- 3. double slit: beats and resolution ------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
k1, k2 = 102, 113  # measured, w = 500 µm, d = 2 mm
S2 = fts.gaussian_band(k_axis, k1, LINE_FWHM_BINS) + fts.gaussian_band(k_axis, k2, LINE_FWHM_BINS)
row2 = fts.interferogram(S2, k_axis, N_PIX, zpd=640, noise=0.01, rng=rng)
ax = axes[0]
ax.plot(p, row2, lw=0.5)
ax.plot(p, fts.envelope(row2) + row2.mean(), "k--", lw=1)
ax.set_xlabel("pixel")
ax.set_ylabel("intensity (norm.)")
ax.set_title(f"Simulated beat pattern, peaks at k = {k1} and {k2}\n"
             f"beat period N/Δk = {N_PIX/(k2-k1):.0f} px")

# two equal Gaussian lines show a dip between them once their separation
# exceeds ~0.85 FWHM; scan the separation and record the dip depth
seps = np.linspace(0, 18, 181)
dips = []
for d_k in seps:
    kk = np.linspace(80, 130, 2001)
    prof = fts.gaussian_band(kk, 105 - d_k / 2, LINE_FWHM_BINS) + \
        fts.gaussian_band(kk, 105 + d_k / 2, LINE_FWHM_BINS)
    centre = prof[np.argmin(np.abs(kk - 105))]
    dips.append(1 - centre / prof.max())
ax = axes[1]
ax.plot(seps, dips)
crit = seps[np.argmax(np.array(dips) > 1e-3)]
ax.axvline(crit, color="k", ls=":", lw=1, label=f"dip appears at Δk ≈ {crit:.1f} bins")
measured = dslit.dropna()
for _, r in measured.iterrows():
    dk = r.k2 - r.k1
    ax.plot(dk, np.interp(dk, seps, dips), "o", ms=6,
            label=f"w={r.slit_width_um:.0f} µm, d={r.slit_gap_um/1000:.1f} mm: Δk={dk:.0f}")
# d = 0.6 mm: estimate from the measured d = 2 mm spacing (~10.5 bins / 2 mm)
dk_06 = 10.5 / 2000 * 600
ax.axvspan(dk_06 - 0.5, dk_06 + 0.5, color="C3", alpha=0.2,
           label=f"d = 0.6 mm (estimated Δk ≈ {dk_06:.1f}): unresolved")
ax.set_xlabel("peak separation Δk (bins)")
ax.set_ylabel("relative dip between peaks")
ax.set_title(f"Resolution for lines of FWHM {LINE_FWHM_BINS} bins")
ax.legend(fontsize=7.5, loc="lower right")
save(fig, "fts_double_slit.png")

# --- results -------------------------------------------------------------------
lines = [
    "# Fourier-transform spectroscopy results",
    "",
    "## Calibration",
    f"- two-point line with intercept: ν = {m_line/1e12:.3f} k + ({c_line/1e12:.1f}) THz "
    "(0 degrees of freedom)",
    f"- proportional model ν = a k: a = {a/1e12:.4f} ± {sa/1e12:.4f} THz/bin, "
    f"χ² = {chi2_cal:.2f} for 1 dof (peak positions ± {PEAK_ERR} bin)",
    f"- implied optical path difference per pixel: {opd_per_px*1e9:.2f} nm "
    f"(mirror displacement {opd_per_px*1e9/2:.2f} nm/px; report: 23.4 nm/px)",
    "",
    "| band | k | λ, proportional (nm) | λ, line (nm) | report error (nm) |",
    "|---|---|---|---|---|",
    *[f"| {b} | {int(kk_)} | {lp:.1f} ± {sl:.1f} | {ll:.1f} | ± {re:.1f} |"
      for b, kk_, lp, sl, ll, re in zip(bands.band, k, lam_prop, slam_prop, lam_line, report_err)],
    "",
    "## White light (Wiener-Khinchin check on measured widths)",
    f"- spectral FWHM {SPECTRAL_FWHM_BINS} bins = {SPECTRAL_FWHM_BINS*a/1e12:.0f} THz",
    f"- predicted envelope FWHM for a Gaussian spectrum: {predicted_env:.1f} px; "
    f"measured {ENVELOPE_FWHM_PX} px ({(ENVELOPE_FWHM_PX/predicted_env-1)*100:+.0f} %)",
    f"- coherence length from the envelope: {ENVELOPE_FWHM_PX*opd_per_px*1e6:.2f} µm "
    "(FWHM of the path-difference envelope)",
    f"- simulation with the measured spectrum: spectrum FWHM {sim_spec_fwhm:.1f} bins, "
    f"envelope FWHM {sim_env_fwhm:.1f} px",
    "",
    "## Double slit",
    f"- two equal lines of FWHM {LINE_FWHM_BINS} bins are resolved for Δk > {crit:.1f} bins",
    f"- d = 0.6 mm corresponds to Δk ≈ {dk_06:.1f} bins (scaled from d = 2 mm) -> unresolved, "
    "as observed",
]
write_results("fourier_spectroscopy.md", lines)
