"""Compare two flux-biased loop elements as three-wave transmission lines.

An rf-SQUID has one junction in an inductive loop. A SNAIL has a small junction
opposite several larger series junctions. Applied flux makes each potential
asymmetric, allowing a 1 GHz pump to couple a 0.45 GHz signal to a 0.55 GHz
idler. The script compares synthetic 128-cell, dispersion-engineered circuits,
not identical fabricated devices or published gain curves.

Gain is a transmitted signal power ratio. Idler conversion is a photon-flux
ratio, and added/limit are input-referred noise in symmetrized quanta. The
derived linear critical-current scale chooses a reference available pump
power; it does not prescribe the physical current through every junction.
Begin with examples/howto/06_rf_squid.py and 07_snail.py. See
docs/three-wave-josephson.md and its Zorin/Frattini/Roudsari references.
"""
from twpa import three_wave as tw

for kind in ['rf-squid', 'snail']:
    chain = tw.example_chain(kind, 128, dispersion_engineered=True)
    power = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
    response = tw.three_wave_response(chain, 1e9, power, .45e9, power*1e-6)
    noise = tw.quantum_noise(chain, 1e9, power, .45e9)
    print(kind, 'gain:', response.gain, 'idler:', response.idler_conversion,
          'added quanta:', noise.added_noise, 'limit:', noise.quantum_limit)
