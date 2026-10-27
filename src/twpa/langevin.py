"""Constant-pump TLS/3WM model extracted from sim_tls_and_spine.py."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from .models import positive


@dataclass
class LangevinResult:
    gain: float
    added_noise: float
    transfer: np.ndarray
    bath_covariance: np.ndarray


def constant_langevin(kp, ks, ki, alpha_s, alpha_i, bias, cells, *, bath_occupation=(0., 0.)):
    """Forward photon-normalized gain and input-referred symmetrized noise.

    This is the source's rotating-frame transfer matrix and Green-function
    noise integral. Matrix exponentials integrate the bath covariance exactly
    for constant coefficients, replacing its 2001-point trapezoidal rule.
    Input idler is vacuum; internal baths may have specified occupations.
    """
    for value in (kp, ks, ki, cells):
        positive(value, 'wavenumbers/cells')
    for value in (alpha_s, alpha_i, *bath_occupation):
        positive(value, 'attenuation/occupation', True)
    if len(bath_occupation) != 2:
        raise ValueError('two bath occupations are required')
    db = kp-ks-ki+bias.xi*bias.pump_current**2/8*(kp-2*ks-2*ki)
    g = bias.epsilon*bias.pump_current*np.sqrt(ks*ki)/4
    m = np.array([[-alpha_s-1j*db/2, 1j*g], [-1j*g, -alpha_i+1j*db/2]])
    u = expm(m*cells)
    diffusion = np.diag(2*np.array([alpha_s, alpha_i])*(np.asarray(bath_occupation)+.5))
    generator = np.zeros((5, 5), complex)
    generator[:4, :4] = np.kron(np.eye(2), m)+np.kron(m.conj(), np.eye(2))
    generator[:4, 4] = diffusion.ravel(order='F')
    v = expm(generator*cells)[:4, 4].reshape((2, 2), order='F')
    gain = float(abs(u[0, 0])**2)
    noise = float((.5*abs(u[0, 1])**2+v[0, 0].real)/gain)
    return LangevinResult(gain, noise, u, v)
