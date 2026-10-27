"""Electrical designs and operating conditions, in SI units."""
from dataclasses import dataclass
from math import isfinite, pi
from scipy.constants import epsilon_0


def positive(value, name, zero=False):
    if not isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")


@dataclass(frozen=True)
class Cell:
    inductance: float
    capacitance: float  # Per finger; the circuit contains TWO fingers.
    stub_inductance: float = 0.0
    participation: float = 0.0
    topology: str = 'shunt-series'

    def __post_init__(self):
        positive(self.inductance, 'inductance')
        positive(self.capacitance, 'capacitance')
        positive(self.stub_inductance, 'stub_inductance', True)
        if not 0 <= self.participation <= 1:
            raise ValueError('participation must lie in [0, 1]')
        if self.topology not in ('shunt-series', 'symmetric'):
            raise ValueError('unknown cell topology')


@dataclass(frozen=True)
class Section:
    cell: Cell
    count: int

    def __post_init__(self):
        if not isinstance(self.count, int) or self.count < 0:
            raise ValueError('count must be a nonnegative integer')


@dataclass(frozen=True)
class Design:
    sections: tuple[Section, ...]
    repeats: int
    tail: tuple[Section, ...] = ()
    name: str = 'custom'

    def __post_init__(self):
        if not isinstance(self.repeats, int) or self.repeats < 1:
            raise ValueError('repeats must be a positive integer')
        if self.period_cells < 1:
            raise ValueError('period must contain at least one cell')

    @property
    def period_cells(self):
        return sum(s.count for s in self.sections)

    @property
    def total_cells(self):
        return self.period_cells * self.repeats + sum(s.count for s in self.tail)


@dataclass(frozen=True)
class Nonlinearity:
    """Phenomenological current-envelope 3WM and Kerr coefficients [1/A, 1/A^2].

    Supply extracted coefficients to avoid imposing a kinetic-inductance law.
    This does not derive Josephson-array or other microscopic device physics.
    """
    three_wave: float
    kerr: float

    def __post_init__(self):
        if not isfinite(self.three_wave) or not isfinite(self.kerr):
            raise ValueError('nonlinear coefficients must be finite')


@dataclass(frozen=True)
class Bias:
    impedance: float = 50.0
    scaling_current: float = 1e-3
    dc_current: float = .2e-3
    pump_current: float = 1e-3 / 60  # Peak current, not RMS.
    kerr_factor: float = 1.0  # TLS source uses 2; MATLAB/nonlinear sources use 1.
    nonlinearity: Nonlinearity | None = None

    def __post_init__(self):
        positive(self.impedance, 'impedance')
        positive(self.scaling_current, 'scaling_current')
        positive(self.dc_current, 'dc_current', True)
        positive(self.pump_current, 'pump_current', True)
        positive(self.kerr_factor, 'kerr_factor', True)

    @property
    def epsilon(self):
        if self.nonlinearity is not None:
            return self.nonlinearity.three_wave
        return 2*self.dc_current/(self.scaling_current**2+self.dc_current**2)

    @property
    def xi(self):
        if self.nonlinearity is not None:
            return self.nonlinearity.kerr
        return self.kerr_factor/(self.scaling_current**2+self.dc_current**2)


@dataclass(frozen=True)
class TLS:
    loss_tangent: float = 0.0  # Intrinsic low-power, zero-temperature value.
    temperature: float = .02
    critical_field: float = 1.0  # V/m
    field_per_amp: float = 0.0  # Peak local oxide E / peak pump current, V/m/A.
    linewidth: float = 1e6  # Hz, phenomenological saturation-hole HWHM.

    def __post_init__(self):
        for name in ('loss_tangent', 'temperature', 'field_per_amp'):
            positive(getattr(self, name), name, True)
        positive(self.critical_field, 'critical_field')
        positive(self.linewidth, 'linewidth')

    @classmethod
    def from_density(cls, density, dipole, relative_permittivity, **kwargs):
        """P0 [1/J/m^3], dipole [C m]; isotropic dipole orientation."""
        for value in (density, dipole, relative_permittivity):
            positive(value, 'density/dipole/permittivity')
        return cls(pi*density*dipole**2/(3*epsilon_0*relative_permittivity), **kwargs)
