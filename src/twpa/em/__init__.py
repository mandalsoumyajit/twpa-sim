"""Geometry-derived passive extraction and circuit back-annotation.

The first supported device is the single-junction resonant JPA. openEMS
runs in its own Python environment; no openEMS imports occur in this package.
"""
from .stack import EMStack, load_stack
from .openems import OpenEMS
from .geometry import capacitor_job
from .backannotation import certify, back_annotate_jpa, read_extraction, load_circuit
from .fasthenry import FastHenry
from .magnetic import MagneticStack, load_magnetic_stack, jpa_magnetic_job
from .magnetic_annotation import certify_magnetic, annotate_coupled_jpa as back_annotate_network_jpa
from .coupled_jpa import CoupledJPA, GalvanicJPA

__all__ = ['EMStack', 'load_stack', 'OpenEMS', 'capacitor_job', 'certify',
           'back_annotate_jpa', 'read_extraction', 'load_circuit', 'FastHenry',
           'MagneticStack','load_magnetic_stack','jpa_magnetic_job','certify_magnetic',
           'back_annotate_network_jpa','CoupledJPA','GalvanicJPA']
