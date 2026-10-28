"""A checked KIT mask reduced with explicitly declared research calibration.

This example does not run full-layout EM. It demonstrates which information
comes from geometry, which is supplied, and what must be recalibrated after
changing the design. All electrical quantities are in SI units.
"""
import argparse
from pathlib import Path
import numpy as np
from twpa.layout import load_process
from twpa.layout.generate import KITLayout,generate
from twpa.models import Bias,TLS
from twpa.em import kit_calibration,back_annotate_kit,load_circuit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);gds=a.output/'kit.gds'
    # Two rows make a small mask with one complete tapered turn. The default
    # finger widths/lengths retain the documented reference geometry.
    generate(KITLayout(rows=2,periods_per_row=1),load_process(),gds)
    # Unloaded/loaded capacitance is PER FINGER; each ladder cell has two.
    # These reference values are demonstration calibration, not a new EM result.
    # The stub term is the original sheet-squares finger correction / 3.
    cal=kit_calibration(gds,capacitance_per_finger_F=(22.32e-15/2,8.38e-15/2),
        stub_inductance_H=(20e-12*(75.5/.5)/3,20e-12*(27.5/.5)/3),
        participation=(.01,.01),
        temperature_K=.02,source='Research reference C/stub values; assumed 1% TLS participation for this demonstration',
        turn_reduction='sheet_squares_matched_periodic')
    # I_dc/I_star=.2 matches the turn design. Sheet L is taken explicitly from
    # the process profile; the electrical model uses L(I)=L(0)*(1+(I/I_star)^2).
    bias=Bias(scaling_current=1e-3,dc_current=.2e-3,pump_current=1e-3/60)
    model,record=back_annotate_kit(gds,cal,bias,output=a.output/'backannotation.json')
    restored=load_circuit(a.output/'backannotation.json')
    # A frequency grid begins at DC so Bloch phase branches can be followed.
    # First inspect passive reflection/transmission before choosing a pump.
    result=restored.linear(np.linspace(0,8e9,801))
    print('Cells:',model.design.total_cells,'period repeats:',model.design.repeats)
    print('Turn sheet-squares error:',record['checks']['turn_sheet_relative_error'])
    print('Transmission at 4 GHz (dB):',result.s21_db[400])
    print(record['scope'])
    # This short mask is an integration example, not a useful high-gain line.
    # TLS participation above and this loss tangent are explicit assumptions.
    # With field_per_amp=0, no pump-saturation field calibration is invented.
    tls=TLS(loss_tangent=1e-3,temperature=.02,field_per_amp=0.)
    coupled=model.coupled([2e9],tls,frequency=np.linspace(0,8e9,801),pump_frequency=6e9)
    print('Coupled TLS/amplifier gain (dB):',float(coupled.gain_db[0]))
    print('Input-referred added quanta:',float(coupled.added_noise[0]))


if __name__=='__main__':main()
