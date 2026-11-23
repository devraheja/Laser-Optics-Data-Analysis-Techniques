"""Gaussian beam relations and the two beam-radius measurement methods."""

import numpy as np
from scipy.special import erf, erfinv

HENE_WAVELENGTH = 632.8e-9  # m


def rayleigh_length(w0, wavelength=HENE_WAVELENGTH, m2=1.0):
    return np.pi * w0**2 / (m2 * wavelength)


def beam_radius(z, w0, wavelength=HENE_WAVELENGTH, m2=1.0):
    """1/e^2 radius at distance z from the waist (waist at z = 0)."""
    return w0 * np.sqrt(1.0 + (z / rayleigh_length(w0, wavelength, m2)) ** 2)


def divergence(w0, wavelength=HENE_WAVELENGTH, m2=1.0):
    """Far-field half angle theta = M^2 lambda / (pi w0)."""
    return m2 * wavelength / (np.pi * w0)


def knife_edge_signal(x, U0, x0, w):
    """Power passing a blade at position x that moves into the beam.

    The detector integrates the Gaussian intensity over the uncovered half
    plane, which gives an error function. w is the 1/e^2 radius.
    """
    return 0.5 * U0 * (1.0 - erf(np.sqrt(2.0) * (x - x0) / w))


def clip_level_factor(eps=0.1):
    """Ratio D_c / w between the clip width at level eps and the 1/e^2 radius.

    Between the two clip positions the fraction 1 - 2 eps of the power passes,
    so erf(D_c / (2 sqrt(2) sigma)) = 1 - 2 eps with w = 2 sigma.
    """
    return np.sqrt(2.0) * erfinv(1.0 - 2.0 * eps)


def clip_width_to_radius(D_c, eps=0.1):
    return D_c / clip_level_factor(eps)
