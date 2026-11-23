"""Fourier-transform spectroscopy on a camera interferogram.

With one interferometer mirror tilted, the optical path difference grows
linearly across the camera. A spectral component then appears as a cosine
fringe whose number of periods across the N-pixel frame, k, is proportional
to its optical frequency. The FFT of one camera row is therefore the
spectrum, with the FFT bin index k as frequency axis.
"""

import numpy as np
from scipy.signal import hilbert

C = 299_792_458.0  # m/s


def gaussian_band(k, k0, fwhm, amplitude=1.0):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((k - k0) / sigma) ** 2)


def interferogram(spectrum, k, n_pix=1280, zpd=640.0, noise=0.0, rng=None):
    """Camera row for a spectrum S(k) given on the grid k (bins = fringes per frame).

    I(p) = sum_k S(k) [1 + cos(2 pi k (p - zpd) / N)] + noise.
    The zero-path-difference position zpd sets where white-light fringes appear.
    """
    rng = np.random.default_rng(rng)
    p = np.arange(n_pix)
    phase = 2.0 * np.pi * np.outer(k, p - zpd) / n_pix
    signal = (spectrum[:, None] * (1.0 + np.cos(phase))).sum(axis=0)
    signal = signal / signal.max()
    return signal + noise * rng.standard_normal(n_pix)


def spectrum_from_interferogram(row):
    """|FFT| of the mean-subtracted row for bins k = 0 ... N/2."""
    row = np.asarray(row, float)
    mag = np.abs(np.fft.rfft(row - row.mean()))
    return np.arange(len(mag)), mag


def envelope(row):
    """Fringe envelope from the analytic signal (Hilbert transform)."""
    row = np.asarray(row, float)
    return np.abs(hilbert(row - row.mean()))


def fwhm(x, y):
    """Full width at half maximum of a single peak, with linear interpolation
    of the two half-maximum crossings."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    i_max = int(np.argmax(y))
    half = y[i_max] / 2.0
    left = i_max
    while left > 0 and y[left] > half:
        left -= 1
    right = i_max
    while right < len(y) - 1 and y[right] > half:
        right += 1
    xl = np.interp(half, [y[left], y[left + 1]], [x[left], x[left + 1]])
    xr = np.interp(half, [y[right], y[right - 1]], [x[right], x[right - 1]])
    return xr - xl


def barycenter(x, weights):
    weights = np.asarray(weights, float)
    return float(np.sum(weights * x) / np.sum(weights))


def gaussian_envelope_width(spectral_fwhm_bins, n_pix):
    """FWHM in pixels of the fringe envelope for a Gaussian spectrum of the given
    FWHM in bins. The envelope is the Fourier transform of the spectrum, and for
    Gaussians FWHM_k * FWHM_p = 4 ln 2 N / pi."""
    return 4.0 * np.log(2.0) * n_pix / (np.pi * spectral_fwhm_bins)
