"""Optional KLayout mask geometry, process contracts and research DRC.

Install ``twpa-sim[layout]`` for geometry operations. Process loading itself
uses only the standard library. Layouts require EM extraction/calibration
before they can represent the numerical amplifier models quantitatively.
"""
from .process import Process, load_process

__all__ = ['Process', 'load_process']
