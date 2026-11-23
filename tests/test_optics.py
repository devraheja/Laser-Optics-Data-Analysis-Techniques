import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from optics import fitting as ft
from optics import gaussian_beam as gb
from optics import interferometry as fts

DATA = ROOT / "data"


# --- Gaussian beam -------------------------------------------------------------

def test_clip_level_factor():
    # 10 % / 90 % clip width is 1.2816 times the 1/e^2 radius
    assert gb.clip_level_factor(0.1) == pytest.approx(1.2816, abs=1e-4)


def test_knife_edge_signal_limits():
    U0, x0, w = 1.7, 1000.0, 800.0
    assert gb.knife_edge_signal(x0, U0, x0, w) == pytest.approx(U0 / 2)
    assert gb.knife_edge_signal(x0 - 5 * w, U0, x0, w) == pytest.approx(U0, rel=1e-6)
    assert gb.knife_edge_signal(x0 + 5 * w, U0, x0, w) == pytest.approx(0, abs=1e-6)
    # 10-90 % width of the erf curve equals the clip width
    x = np.linspace(x0 - 3 * w, x0 + 3 * w, 200001)
    U = gb.knife_edge_signal(x, U0, x0, w)
    width = x[np.argmin(abs(U - 0.1 * U0))] - x[np.argmin(abs(U - 0.9 * U0))]
    assert width == pytest.approx(gb.clip_level_factor(0.1) * w, rel=1e-4)


def test_beam_radius_at_rayleigh_length():
    w0 = 210e-6
    zr = gb.rayleigh_length(w0)
    assert gb.beam_radius(zr, w0) == pytest.approx(np.sqrt(2) * w0)
    assert gb.beam_radius(100 * zr, w0) == pytest.approx(gb.divergence(w0) * 100 * zr, rel=1e-4)


def test_waist_fit_recovers_synthetic_beam():
    rng = np.random.default_rng(0)
    w0 = 205e-6
    z = np.array([0.935, 1.745, 2.06])
    w = gb.beam_radius(z, w0) * (1 + 0.002 * rng.standard_normal(3))
    f = ft.fit(lambda z, a: gb.beam_radius(z, a), z, w, [2e-4], sigma=0.002 * w)
    assert f.params[0] == pytest.approx(w0, rel=5e-3)


# --- fitting helpers -----------------------------------------------------------

def test_threshold_error_with_covariance_matches_direct_fit():
    d = pd.read_csv(DATA / "yag_pump_diode.csv")
    d = d[d.I_mA > 200]
    lin = ft.fit(ft.line, d.I_mA, d.P_mW, [1, -200])
    hinge = ft.fit(ft.hinge, d.I_mA, d.P_mW, [1, 200])
    x0, s_full, s_nocov = ft.propagate_ratio(*lin.params, lin.cov)
    assert x0 == pytest.approx(hinge.params[1])
    assert s_full == pytest.approx(hinge.errors[1], rel=1e-4)
    assert s_full < s_nocov  # slope and intercept are anti-correlated


def test_report_values_reproduced():
    # pump diode: s = 1.066 W/A, I_th = 199.8 mA (lab report)
    d = pd.read_csv(DATA / "yag_pump_diode.csv")
    d = d[d.I_mA > 200]
    s, b = ft.fit(ft.line, d.I_mA, d.P_mW, [1, -200]).params
    assert s == pytest.approx(1.066, abs=1e-3)
    assert -b / s == pytest.approx(199.8, abs=0.1)
    # laser: slope efficiency 0.120, threshold 116.5 mW (linear model)
    L = pd.read_csv(DATA / "yag_laser_output.csv")
    eta, P_th = ft.fit(ft.hinge, s * L.I_mA + b, L.P_out_mW, [0.1, 100]).params
    assert eta == pytest.approx(0.120, abs=1e-3)
    assert P_th == pytest.approx(116.5, abs=0.2)


def test_clip_level_radii_match_report():
    cl = pd.read_csv(DATA / "hene_clip_level.csv")
    w = gb.clip_width_to_radius(cl.x10_um - cl.x90_um)
    assert np.allclose(w, [905, 1654, 1904, 921, 1779, 2044], atol=1.0)


# --- interferometry ------------------------------------------------------------

def test_fft_recovers_line_position():
    k = np.arange(1, 640)
    row = fts.interferogram(fts.gaussian_band(k, 113, 1.0), k, zpd=300.0)
    k_fft, S = fts.spectrum_from_interferogram(row)
    assert k_fft[np.argmax(S)] == 113


def test_fwhm_of_gaussian():
    x = np.linspace(-50, 50, 10001)
    assert fts.fwhm(x, fts.gaussian_band(x, 3.0, 12.0)) == pytest.approx(12.0, rel=1e-4)


def test_envelope_width_is_fourier_partner_of_spectral_width():
    k = np.arange(1, 640)
    row = fts.interferogram(fts.gaussian_band(k, 80, 26.0), k, n_pix=1280, zpd=640)
    env = fts.envelope(row)
    measured = fts.fwhm(np.arange(1280), env)
    assert measured == pytest.approx(fts.gaussian_envelope_width(26.0, 1280), rel=0.03)


def test_proportional_calibration_is_consistent():
    cal = pd.read_csv(DATA / "michelson_calibration.csv")
    a_i = fts.C / (cal.wavelength_nm * 1e-9) / cal.peak_index
    # the two filters agree on THz per bin within the +-0.5 bin peak resolution
    assert abs(a_i[0] - a_i[1]) / a_i.mean() < 0.5 / 92 + 0.5 / 113
