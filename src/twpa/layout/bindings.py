"""Constrained electrical-to-layout sizing (no inverse EM synthesis)."""
import math

from .generate import JosephsonLayout
from .process import positive


def from_circuit(circuit, process, *, capacitor_width=12., relative_tolerance=.01, **geometry):
    """Size a bare ParallelJPA or uniform 4WM Chain from Ic, Cj and shunt C.

    Dimensions are rounded to the explicit process grid. Incompatible junction
    area/current/capacitance targets are rejected, including after rounding.
    No inductance, resonant loading, loss resistor or port coupling is synthesized.
    Use geometric design inputs for SQUID/SNAIL loops pending EM extraction.
    """
    from ..jpa_circuit import ParallelJPA
    from ..jtwpa import Chain
    positive(relative_tolerance,'relative_tolerance')
    if relative_tolerance > .1:
        raise ValueError('Sizing tolerance must be <=10%')
    if isinstance(circuit,ParallelJPA):
        family='jpa'; cells=1; shunt=circuit.shunt_capacitance
    elif isinstance(circuit,Chain):
        if circuit.resonator is not None:
            raise ValueError('Resonant phase-matching loads need additional layout primitives')
        family='jtwpa4wm'; cells=circuit.cells; shunt=circuit.ground_capacitance
        geometry.setdefault('cell_pitch',circuit.cell_length*1e6)
    else:
        raise TypeError('Automatic sizing currently supports ParallelJPA and uniform 4WM Chain')
    positive(capacitor_width,'capacitor_width'); process.coordinate(capacitor_width)
    jc=process.parameter('critical_current_density','uA/um2')*1e-6
    cj=process.parameter('junction_specific_capacitance','fF/um2')*1e-15
    cs=process.parameter('capacitor_specific_capacitance','fF/um2')*1e-15
    def snap(x): return round(x/process.grid)*process.grid
    width=snap(math.sqrt(circuit.critical_current/jc))
    height=snap(shunt/(cs*capacitor_width))
    achieved=(width**2*jc,width**2*cj,capacitor_width*height*cs)
    targets=(circuit.critical_current,circuit.junction_capacitance,shunt)
    for name,value,target in zip(('critical_current','junction_capacitance','shunt_capacitance'),achieved,targets):
        if abs(value/target-1)>relative_tolerance:
            raise ValueError(f'{name} incompatible with process/geometry: target={target:g}, realized={value:g}; '
                             'Ic and Cj share the same fabricated junction area')
    reserved={'family','cells','junction_width','capacitor_width','capacitor_height'}
    if reserved & geometry.keys():
        raise ValueError('Cannot override an electrically sized dimension')
    layout=JosephsonLayout(family=family,cells=cells,junction_width=width,
                           capacitor_width=capacitor_width,capacitor_height=height,**geometry)
    layout.validate(process)
    return layout
