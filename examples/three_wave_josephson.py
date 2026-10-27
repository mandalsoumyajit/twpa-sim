"""Synthetic rf-SQUID and SNAIL 3WM gain and quantum noise."""
from twpa import three_wave as tw

for kind in ['rf-squid', 'snail']:
    chain = tw.example_chain(kind, 128, dispersion_engineered=True)
    power = .5*chain.port_impedance*(.18*chain.linear_chain.critical_current)**2
    response = tw.three_wave_response(chain, 1e9, power, .45e9, power*1e-6)
    noise = tw.quantum_noise(chain, 1e9, power, .45e9)
    print(kind, 'gain:', response.gain, 'idler:', response.idler_conversion,
          'added quanta:', noise.added_noise, 'limit:', noise.quantum_limit)
