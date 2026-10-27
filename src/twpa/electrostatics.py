"""Rectangular finite-volume electrostatics for supplied dielectric/electrode maps."""
from dataclasses import dataclass
import warnings
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve, MatrixRankWarning
from scipy.constants import epsilon_0


@dataclass
class ElectrostaticResult:
    potential: np.ndarray
    capacitance_per_length: float  # F/m along the invariant third dimension.
    participation: float


def electrostatic(relative_permittivity, fixed_potential, oxide, *, dx, dy, voltage=1., periodic_x=False):
    """NaN denotes free potential; other entries are Dirichlet electrodes.

    Exterior faces are insulating except optional periodic x boundaries.
    Permittivity is cell centered, face conductance uses harmonic averaging.
    Convert C' to a finger capacitance with the appropriate invariant length.
    The caller supplies the oxide geometry; refine the mesh to convergence.
    """
    er = np.asarray(relative_permittivity, float)
    fixed = np.asarray(fixed_potential, float)
    oxide = np.asarray(oxide, bool)
    if er.ndim != 2 or min(er.shape) < 2 or er.shape != fixed.shape or er.shape != oxide.shape:
        raise ValueError('maps must have identical 2D shapes, at least 2 by 2')
    if np.any(~np.isfinite(er)) or np.any(er <= 0) or np.any(np.isinf(fixed)):
        raise ValueError('invalid permittivity or potential map')
    if not all(np.isfinite(v) and v > 0 for v in (dx, dy, voltage)):
        raise ValueError('dx, dy and voltage must be positive and finite')
    boundary = np.isfinite(fixed)
    if len(np.unique(fixed[boundary])) < 2:
        raise ValueError('at least two electrode voltages are required')
    if not np.isclose(np.ptp(fixed[boundary]), voltage):
        raise ValueError('voltage must equal electrode potential difference')
    ny, nx = er.shape
    a = lil_matrix((er.size, er.size)); b = np.zeros(er.size)
    faces = []
    for y in range(ny):
        for x in range(nx):
            i = y*nx+x
            neighbors = []
            if x+1 < nx:
                neighbors.append((y, x+1, dy/dx))
            elif periodic_x:
                neighbors.append((y, 0, dy/dx))
            if y+1 < ny:
                neighbors.append((y+1, x, dx/dy))
            for yy, xx, aspect in neighbors:
                j = yy*nx+xx
                e1, e2 = er[y, x], er[yy, xx]
                g = epsilon_0*2*e1*e2/(e1+e2)*aspect
                # Energy divides between half cells in proportion to 1/epsilon.
                fraction = (oxide[y, x]/e1+oxide[yy, xx]/e2)/(1/e1+1/e2)
                faces.append((i, j, g, fraction))
                a[i, i] += g; a[j, j] += g
                a[i, j] -= g; a[j, i] -= g
    for i in np.flatnonzero(boundary):
        a.rows[i] = [i]; a.data[i] = [1.]; b[i] = fixed.ravel()[i]
    with warnings.catch_warnings():
        warnings.simplefilter('error', MatrixRankWarning)
        phi = spsolve(a.tocsr(), b)
    energy = sum(.5*g*(phi[i]-phi[j])**2 for i, j, g, _ in faces)
    oxide_energy = sum(.5*g*(phi[i]-phi[j])**2*p for i, j, g, p in faces)
    if not np.isfinite(energy) or energy <= 0:
        raise ValueError('invalid electrostatic energy')
    return ElectrostaticResult(phi.reshape(er.shape), 2*energy/voltage**2, oxide_energy/energy)
