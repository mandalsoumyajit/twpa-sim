"""Capture portable numerical evidence from completed EM and JoSIM runs."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import numpy as np


def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))


def capture(em_root,final,fine,coarse,short,output,figure=None):
    extraction=load(em_root/'extraction-certified.json')
    results=[load(p/'summary.json') for p in (final,fine,coarse,short)]
    final_data=deepcopy(results[0])
    em={key:extraction[key] for key in ('source','stack','stack_sha256','geometry','settings','simulator',
         'mesh','frequency_hz','capacitance_by_frequency_F','capacitance_F',
         'parallel_C_fit_max_relative_residual','port_tail_relative','negative_conductance_relative',
         'fit_accepted','convergence','scope','artifact_sha256')}
    em['worker_sha256']=extraction['run']['worker_sha256']
    em['source_artifact_sha256']=em.pop('artifact_sha256')
    def difference(a,b):
        return dict(signal_gain_db=float(max(abs(10*np.log10(np.asarray(a['pumped']['measured']['gain'])/
                                                                          np.asarray(b['pumped']['measured']['gain']))))),
                    idler_conversion_db=float(max(abs(10*np.log10(np.asarray(a['pumped']['measured']['idler_conversion'])/
                                                                                np.asarray(b['pumped']['measured']['idler_conversion']))))))
    sensitivity=dict(timestep_final_vs_fine=difference(results[0],results[1]),
                     timestep_fine_vs_coarse=difference(results[1],results[2]),
                     window_coarse_vs_short=difference(results[2],results[3]))
    rows=[dict(transient=r['transient'],passed=r['passed'],
               passive_error=r['passive']['max_complex_s11_error'],
               gain_error_db=r['pumped']['gain_error_db'],idler_error_db=r['pumped']['idler_error_db']) for r in results]
    data=dict(schema_version=1,date='2026-10-28',extraction=em,response=final_data,
              refinement_cases=rows,sensitivity=sensitivity,
              scope='Component-level PEC capacitor extraction and classical annotated-circuit validation. Not full-layout or fabricated-device validation.')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    if figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        freq=np.asarray(em['frequency_hz'])/1e9
        ax[0,0].plot(freq,np.asarray(em['capacitance_by_frequency_F'])*1e15,'o-',label='openEMS')
        ax[0,0].axhline(em['source']['nominal_capacitance_F']*1e15,color='0.5',ls='--',label='Nominal annotation estimate')
        ax[0,0].set(xlabel='Frequency (GHz)',ylabel='Capacitance (fF)',title='Actual electrode coupon');ax[0,0].legend(fontsize=8)
        ax[0,1].bar(['Nominal','EM annotated'],np.array([final_data['nominal_resonance_Hz'],final_data['extracted_resonance_Hz']])/1e9,
                    color=['#9da6b0','#196a9d'])
        ax[0,1].set(ylabel='Resonance (GHz)',title='Circuit changes without gain fitting')
        passive=final_data['passive']; f=np.asarray(passive['measured']['signal_frequency'])/1e9
        for part,style in [('real','-'),('imag','--')]:
            ax[1,0].plot(f,passive['model_'+part],style,label='Model '+part)
            ax[1,0].plot(f,passive['measured']['signal_'+part],'o',fillstyle='none',label='JoSIM '+part)
        ax[1,0].set(xlabel='Probe (GHz)',ylabel='Complex S11',title='Passive reflection');ax[1,0].legend(fontsize=8)
        pumped=final_data['pumped']; fs=(np.asarray(pumped['measured']['signal_frequency'])-pumped['pump_frequency_Hz'])/1e6
        for key,actual,label in [('model_gain','gain','Signal'),('model_idler_conversion','idler_conversion','Idler')]:
            ax[1,1].plot(fs,10*np.log10(pumped[key]),'-',label=label+' model')
            ax[1,1].plot(fs,10*np.log10(pumped['measured'][actual]),'o',fillstyle='none',label=label+' JoSIM')
        ax[1,1].set(xlabel='Probe offset from pump (MHz)',ylabel='Power / photon conversion (dB)',title='Pumped response');ax[1,1].legend(fontsize=8)
        fig.suptitle('JPA: GDS capacitor → openEMS → annotated RCSJ circuit',fontsize=13)
        for axis in ax.flat: axis.grid(alpha=.2)
        figure.parent.mkdir(parents=True,exist_ok=True);fig.savefig(figure,dpi=150);plt.close(fig)
    print(json.dumps(dict(extraction_F=em['capacitance_F'],sensitivity=sensitivity,refinement_cases=rows),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('em-root','final','fine','coarse','short','output'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--figure',type=Path)
    a=p.parse_args(); capture(a.em_root,a.final,a.fine,a.coarse,a.short,a.output,a.figure)
