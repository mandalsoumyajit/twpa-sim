"""Geometry-derived passive extraction and circuit back-annotation.

Finite-terminal Josephson networks retain nonuniform devices and signed
London coupling. KIT uses a separately declared calibration-backed periodic
reduction. openEMS runs in its own environment; it is not imported here.
"""
from .stack import EMStack, load_stack
from .openems import OpenEMS
from .geometry import capacitor_job
from .backannotation import certify, back_annotate_jpa, read_extraction, load_circuit
from .fasthenry import FastHenry
from .magnetic import MagneticStack, load_magnetic_stack, jpa_magnetic_job
from .magnetic_annotation import certify_magnetic, annotate_coupled_jpa as back_annotate_network_jpa
from .coupled_jpa import CoupledJPA, GalvanicJPA
from .magnetic import layout_magnetic_job
from .layout_annotation import back_annotate_layout
from .network import LayoutNetwork
from .backannotation import reuse_capacitor
from .kit_annotation import KITPeriodicModel,kit_calibration,back_annotate_kit

__all__ = ['EMStack', 'load_stack', 'OpenEMS', 'capacitor_job', 'certify',
           'back_annotate_jpa', 'read_extraction', 'load_circuit', 'FastHenry',
           'MagneticStack','load_magnetic_stack','jpa_magnetic_job','certify_magnetic',
           'back_annotate_network_jpa','CoupledJPA','GalvanicJPA','layout_magnetic_job',
           'back_annotate_layout','LayoutNetwork','reuse_capacitor','KITPeriodicModel','kit_calibration','back_annotate_kit']
