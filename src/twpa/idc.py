"""Parameterized IDC edge Laplace solver ported from generate_reorganized_figures.py."""
from dataclasses import dataclass, replace
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.constants import epsilon_0
from .models import positive


@dataclass(frozen=True)
class IDCGeometry:
    finger_width: float = .5e-6
    gap: float = .25e-6
    metal_thickness: float = .06e-6
    substrate_permittivity: float = 11.45
    oxide_permittivity: float = 3.9
    substrate_depth: float = .2e-6
    air_height: float = .2e-6
    mesh: float = 2.5e-9

    def __post_init__(self):
        for name, value in vars(self).items():
            positive(value, name)
        if self.air_height <= self.metal_thickness:
            raise ValueError('air_height must exceed metal_thickness')


@dataclass
class IDCResult:
    x: np.ndarray
    y: np.ndarray
    potential: np.ndarray
    energy_density: np.ndarray
    energy_per_length: float
    participation: float


def idc_fields(geometry=IDCGeometry(), oxide_thickness=60e-9):
    """Source boundary conditions: +0.5 V finger, 0 V gap-center symmetry plane.

    Neumann left/top/bottom boundaries. Mesh and domain must converge for the
    intended geometry. Output energy covers this HALF-gap edge domain only.
    Use coated/bare energy ratios to scale calibrated full-cell capacitances.
    """
    positive(oxide_thickness, 'oxide_thickness', True)
    g = geometry
    edge = g.finger_width/2
    if oxide_thickness+g.metal_thickness >= g.air_height:
        raise ValueError('oxide reaches air-domain boundary; increase air_height')
    lx = edge+g.gap/2
    nx = int(round(lx/g.mesh))+1
    ny = int(round((g.substrate_depth+g.air_height)/g.mesh))+1
    if nx < 3 or ny < 3:
        raise ValueError('mesh too coarse')
    x = np.linspace(0, lx, nx)
    y = np.linspace(-g.substrate_depth, g.air_height, ny)
    dx, dy = np.diff(x), np.diff(y)
    xc, yc = (x[:-1]+x[1:])/2, (y[:-1]+y[1:])/2
    xx, yy = np.meshgrid(xc, yc, indexing='ij')
    er = np.ones_like(xx); er[yy < 0] = g.substrate_permittivity
    gap_ox = (yy >= 0) & (yy <= oxide_thickness) & (xx >= edge)
    top_ox = (yy >= g.metal_thickness) & (yy <= g.metal_thickness+oxide_thickness) & (xx <= edge)
    side_ox = (yy >= 0) & (yy <= g.metal_thickness+oxide_thickness) & (xx >= edge) & (xx <= edge+oxide_thickness)
    oxide = gap_ox | top_ox | side_ox
    if oxide_thickness == 0:
        oxide[:] = False
    er[oxide] = g.oxide_permittivity
    a = lil_matrix((nx*ny, nx*ny)); b = np.zeros(nx*ny)
    metal = np.zeros((nx, ny), bool)
    ix = int(round(edge/(lx/(nx-1))))
    jm = np.where((y >= -1e-12) & (y <= g.metal_thickness+1e-12))[0]
    metal[:ix+1, jm] = True
    for i in range(nx):
        for j in range(ny):
            idx = i*ny+j
            if metal[i, j]:
                a[idx, idx] = 1; b[idx] = .5
            elif i == nx-1:
                a[idx, idx] = 1
            elif j == 0 or j == ny-1:
                a[idx, idx] = 1; a[idx, i*ny+(1 if j == 0 else ny-2)] = -1
            else:
                ds, dn = dy[j-1], dy[j]
                if i == 0:
                    de = dx[0]
                    ce = (er[0, j-1]*ds+er[0, j]*dn)/de
                    cs, cn = er[0, j-1]*de/ds, er[0, j]*de/dn
                    a[idx, idx] = -(ce+cs+cn)
                    a[idx, ny+j] = ce
                else:
                    dw, de = dx[i-1], dx[i]
                    cw = .5*(er[i-1, j-1]*ds+er[i-1, j]*dn)/dw
                    ce = .5*(er[i, j-1]*ds+er[i, j]*dn)/de
                    cs = .5*(er[i-1, j-1]*dw+er[i, j-1]*de)/ds
                    cn = .5*(er[i-1, j]*dw+er[i, j]*de)/dn
                    a[idx, idx] = -(cw+ce+cs+cn)
                    a[idx, (i-1)*ny+j] = cw; a[idx, (i+1)*ny+j] = ce
                a[idx, idx-1] = cs; a[idx, idx+1] = cn
    v = spsolve(a.tocsc(), b).reshape(nx, ny)
    ex = -(v[1:, :-1]+v[1:, 1:]-v[:-1, :-1]-v[:-1, 1:])/(2*dx[:, None])
    ey = -(v[:-1, 1:]+v[1:, 1:]-v[:-1, :-1]-v[1:, :-1])/(2*dy[None, :])
    energy = .5*epsilon_0*er*(ex**2+ey**2)
    in_metal = (xx <= edge) & (yy >= 0) & (yy <= g.metal_thickness)
    energy[in_metal] = np.nan
    area = dx[:, None]*dy[None, :]
    total = np.nansum(energy*area)
    participation = np.nansum(np.where(oxide, energy*area, 0))/total
    return IDCResult(x, y, v, energy, float(total), float(participation))


def oxide_design(design, geometry=IDCGeometry(), oxide_thickness=60e-9):
    """Wire the source's Laplace energy-ratio prescription into the circuit.

    The same edge participation/capacitance ratio applies to every section.
    For unlike cross-sections, extract and assign each Cell separately.
    """
    bare = idc_fields(geometry, 0.)
    coated = bare if oxide_thickness == 0 else idc_fields(geometry, oxide_thickness)
    ratio = coated.energy_per_length/bare.energy_per_length
    def coat(section):
        cell = replace(section.cell, capacitance=section.cell.capacitance*ratio,
                       participation=coated.participation)
        return replace(section, cell=cell)
    result = replace(design, sections=tuple(coat(s) for s in design.sections),
                     tail=tuple(coat(s) for s in design.tail))
    return result, bare, coated
