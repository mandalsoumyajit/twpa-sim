"""Physical RCSJ resonator decks and classical reduced-model parameters."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, elementary_charge as e, physical_constants
from scipy.optimize import brentq
from .jpa import Resonator
from .models import positive

PHI0 = physical_constants['mag. flux quantum'][0]
PHIBAR = PHI0/(2*np.pi)


@dataclass(frozen=True)
class Tone:
    frequency: float
    power: float
    phase: float = 0.  # Radians in sin(omega*t+phase).

    def __post_init__(self):
        positive(self.frequency, 'tone frequency')
        positive(self.power, 'tone power', True)
        if not np.isfinite(self.phase):
            raise ValueError('tone phase must be finite')


@dataclass(frozen=True)
class Transient:
    step: float = 1e-12
    stop: float = 6e-6
    start: float = 2e-6
    output_step: float = 20e-12

    def __post_init__(self):
        for value in (self.step, self.stop, self.output_step):
            positive(value, 'transient times')
        if not np.isfinite(self.start) or not 0 <= self.start < self.stop or self.output_step < self.step:
            raise ValueError('require 0 <= start < stop and output_step >= step')


@dataclass(frozen=True)
class ParallelJPA:
    critical_current: float
    shunt_capacitance: float
    junction_capacitance: float = 1e-12
    junction_resistance: float = 1e9
    port_impedance: float = 50.
    internal_resistance: float = 5000.


    def __post_init__(self):
        for name, value in vars(self).items():
            positive(value, name)

    @property
    def capacitance(self):
        return self.shunt_capacitance+self.junction_capacitance

    @property
    def inductance(self):
        return PHIBAR/self.critical_current

    def reduced_model(self):
        c = self.capacitance
        return Resonator(1/(2*np.pi*np.sqrt(self.inductance*c)), 1/(self.port_impedance*c),
                         (1/self.internal_resistance+1/self.junction_resistance)/c,
                         -e**2/(2*c*hbar))

    def components(self):
        return [f'.model JJMOD jj(RTYPE=0 IC={self.critical_current:.17g} '
                f'C={self.junction_capacitance:.17g} RN={self.junction_resistance:.17g} '
                f'R0={self.junction_resistance:.17g})',
                'BJJ PORT 0 JJMOD', f'CSH PORT 0 {self.shunt_capacitance:.17g}',
                f'RINT PORT 0 {self.internal_resistance:.17g}']


@dataclass(frozen=True)
class SQUIDJPA:
    critical_current: float  # Each of two identical junctions.
    shunt_capacitance: float
    flux_bias: float = .2  # External flux / flux quantum.
    arm_inductance: float = 1e-12
    coil_inductance: float = 1e-9
    coupling: float = .2  # +/- mutual coupling to the two arms.
    junction_capacitance: float = .5e-12
    junction_resistance: float = 1e9
    port_impedance: float = 50.
    internal_resistance: float = 5000.
    bias_ramp: float = 100e-9

    def __post_init__(self):
        for name in ('critical_current', 'shunt_capacitance', 'arm_inductance', 'coil_inductance',
                     'junction_capacitance', 'junction_resistance', 'port_impedance', 'internal_resistance', 'bias_ramp'):
            positive(getattr(self, name), name)
        if not np.isfinite(self.flux_bias) or not abs(self.flux_bias) < .45:
            raise ValueError('initial SQUID adapter requires |flux_bias| < 0.45')
        if not 0 < self.coupling < 1/np.sqrt(2):
            raise ValueError('coupling must give a positive-definite three-inductor matrix')

    @property
    def mutual_inductance(self):
        return self.coupling*np.sqrt(self.arm_inductance*self.coil_inductance)

    def junction_phase(self, flux=None):
        bias = self.flux_bias if flux is None else flux
        beta = self.arm_inductance*self.critical_current/PHIBAR
        return brentq(lambda x: x+beta*np.sin(x)-np.pi*bias, -np.pi/2, np.pi/2)

    def _mode(self, flux=None):
        lj = PHIBAR/(self.critical_current*np.cos(self.junction_phase(flux)))
        def admittance(w):
            yj = 1/self.junction_resistance+1j*w*self.junction_capacitance+1/(1j*w*lj)
            return 1j*w*self.shunt_capacitance+2/(1j*w*self.arm_inductance+1/yj)
        guess = 1/np.sqrt((lj+self.arm_inductance)/2*(self.shunt_capacitance+2*self.junction_capacitance))
        w = brentq(lambda w: admittance(w).imag, guess*.8, guess*1.2)
        dw = w*1e-5
        ceff = (admittance(w+dw).imag-admittance(w-dw).imag)/(4*dw)
        return w, ceff, admittance(w).real

    def reduced_model(self):
        w, c, conductance = self._mode()
        # Kerr is not used by flux_response; leave it zero for this linear adapter.
        return Resonator(w/(2*np.pi), 1/(self.port_impedance*c),
                         (1/self.internal_resistance+conductance)/c)

    def frequency_sensitivity(self):
        """d omega_mode / d(flux/Phi0), angular frequency per flux quantum."""
        d = 1e-5
        return (self._mode(self.flux_bias+d)[0]-self._mode(self.flux_bias-d)[0])/(2*d)

    def pair_drive(self, flux_amplitude, phase=0.):
        positive(flux_amplitude, 'flux_amplitude', True)
        if abs(self.flux_bias)+flux_amplitude >= .45 or not np.isfinite(phase):
            raise ValueError('flux drive exceeds the supported bias branch or phase is nonfinite')
        derivative = self.frequency_sensitivity()
        # Voltage-normalized port envelope, sine flux pump. No fitting to pumped gain.
        return -.5*derivative*flux_amplitude*np.exp(-1j*phase)

    def components(self, flux_pump=None):
        lines = [f'.model JJMOD jj(RTYPE=0 IC={self.critical_current:.17g} '
                 f'C={self.junction_capacitance:.17g} RN={self.junction_resistance:.17g} '
                 f'R0={self.junction_resistance:.17g})',
                 'BJ1 PORT ARM1 JJMOD', 'BJ2 PORT ARM2 JJMOD',
                 f'LARM1 ARM1 0 {self.arm_inductance:.17g}',
                 f'LARM2 ARM2 0 {self.arm_inductance:.17g}',
                 f'LCOIL COIL 0 {self.coil_inductance:.17g}',
                 f'KONE LCOIL LARM1 {self.coupling:.17g}',
                 f'KTWO LCOIL LARM2 {-self.coupling:.17g}',
                 f'IBIAS 0 COIL PWL(0 0 {self.bias_ramp:.17g} '
                 f'{self.flux_bias*PHI0/(2*self.mutual_inductance):.17g})',
                 f'CSH PORT 0 {self.shunt_capacitance:.17g}',
                 f'RINT PORT 0 {self.internal_resistance:.17g}']
        if flux_pump is not None:
            frequency, amplitude, phase = flux_pump
            lines.append(f'IFLUX 0 COIL SIN(0 {amplitude*PHI0/(2*self.mutual_inductance):.17g} '
                         f'{frequency:.17g} 0 0 {phase:.17g})')
        return lines


def jpa_netlist(circuit, tones, transient=Transient(), *, flux_pump=None):
    """Self-contained physical JJ deck with Thevenin port and explicit damping.

    Port voltage source has twice the incident-wave peak voltage. The two JJ
    capacitances in a SQUID and the single-JJ C are part of the RCSJ models;
    CSH is additional capacitance. No thermal-noise sources are enabled.
    """
    if flux_pump is not None and not isinstance(circuit, SQUIDJPA):
        raise ValueError('flux pumping requires a SQUID circuit')
    if flux_pump is not None:
        frequency, amplitude, phase = flux_pump
        positive(frequency, 'flux pump frequency')
        circuit.pair_drive(amplitude, phase)  # Validate supported flux excursion.
    tones = tuple(tones)
    if not tones:
        raise ValueError('provide at least one tone, even if power is zero')
    lines = ['* TWPA Sim: classical Josephson resonator circuit']
    for i, tone in enumerate(tones):
        nxt = f'SRC{i+1}' if i+1<len(tones) else '0'
        peak = np.sqrt(8*circuit.port_impedance*tone.power)
        lines.append(f'VDRIVE{i} SRC{i} {nxt} SIN(0 {peak:.17g} '
                     f'{tone.frequency:.17g} 0 0 {tone.phase:.17g})')
    lines.append(f'RPORT SRC0 PORT {circuit.port_impedance:.17g}')
    lines += circuit.components(flux_pump) if isinstance(circuit, SQUIDJPA) else circuit.components()
    junction_traces = 'P(BJ1) P(BJ2) I(LCOIL)' if isinstance(circuit, SQUIDJPA) else 'P(BJJ)'
    lines += [f'.tran {transient.step:.17g} {transient.stop:.17g} {transient.start:.17g} '
              f'{transient.output_step:.17g}', f'.print V(PORT) V(SRC0) {junction_traces}', '.end']
    return '\n'.join(lines)+'\n'
