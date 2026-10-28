"""A magnetically coupled JPA with extracted interconnect and pickup loop.

Five-state full-sine circuit; no fitted external linewidth or added lumped
kinetic inductance. The 2x2 London matrix already contains both kinds of
stored inductive energy, the explicit return, and signed mutual coupling.
"""
from dataclasses import dataclass,asdict
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp

from ..jpa_circuit import PHIBAR
from ..models import positive


@dataclass(frozen=True)
class CoupledJPA:
    critical_current: float
    shunt_capacitance: float
    junction_capacitance: float
    interconnect_inductance: float
    pickup_inductance: float
    mutual_inductance: float
    port_impedance: float = 50.
    internal_resistance: float = 1e8
    junction_resistance: float = 1e9

    def __post_init__(self):
        for name,value in asdict(self).items():
            if name!='mutual_inductance':positive(value,name)
        if not np.isfinite(self.mutual_inductance) or np.linalg.eigvalsh(self.inductance_matrix)[0]<=0:
            raise ValueError('A symmetric positive-definite signed inductance matrix is required')

    @property
    def inductance_matrix(self):
        return np.array([[self.interconnect_inductance,self.mutual_inductance],
                         [self.mutual_inductance,self.pickup_inductance]])

    @property
    def plasma_omega(self):
        return np.sqrt(self.critical_current/(PHIBAR*self.junction_capacitance))

    def linear_system(self):
        """Scaled time tau=omega_J*t; state=(phase,Vj,Vc,Itank,Ipickup)."""
        w=self.plasma_omega; inv=np.linalg.inv(self.inductance_matrix);lj=PHIBAR/self.critical_current
        a=np.zeros((5,5));a[0,1]=1;a[1,0]=-1
        a[1,1]=-1/(self.junction_resistance*self.junction_capacitance*w);a[1,3]=1
        a[2,2]=-1/(self.internal_resistance*self.shunt_capacitance*w)
        a[2,3]=-self.junction_capacitance/self.shunt_capacitance
        a[3:,1]=-inv[:,0]*lj;a[3:,2]=inv[:,0]*lj
        a[3:,4]=-inv[:,1]*self.port_impedance/w
        source=np.zeros(5);source[3:]=inv[:,1]/(self.critical_current*w)
        return a,source

    def modes(self):
        a,b=self.linear_system();values,vectors=np.linalg.eig(a)
        modes=[]
        for k,v in enumerate(values):
            if v.imag>1e-9:
                modes.append(dict(frequency=float(v.imag*self.plasma_omega/(2*np.pi)),
                                  kappa=float(-2*v.real*self.plasma_omega),
                                  phase_participation=float(abs(vectors[0,k])/np.linalg.norm(vectors[:,k]))))
        return sorted(modes,key=lambda m:m['frequency'])

    def passive(self,frequencies):
        a,b=self.linear_system();result=[]
        for f in np.atleast_1d(frequencies):
            x=np.linalg.solve(-1j*2*np.pi*f/self.plasma_omega*np.eye(5)-a,b)
            result.append(1-2*self.port_impedance*self.critical_current*(self.port_incidence@x[3:]))
        return np.array(result)

    @property
    def port_incidence(self):return np.array([0.,1.])

    def rhs_scaled(self,state,source_voltage=0.):
        x=np.asarray(state);a,b=self.linear_system();value=a@x+b*source_voltage
        value[1]-=np.sin(x[0])-x[0]
        return value

    def energy(self,state):
        phi,uj,uc,i,k=np.asarray(state)
        ej=PHIBAR*self.critical_current;lj=PHIBAR/self.critical_current
        currents=np.array([i,k])
        return ej*(1-np.cos(phi)+.5*uj**2+.5*self.shunt_capacitance/self.junction_capacitance*uc**2+
                   .5*currents@self.inductance_matrix@currents/lj)

    def components(self):
        k=self.mutual_inductance/np.sqrt(self.interconnect_inductance*self.pickup_inductance)
        lines=[f'.model JJMOD jj(RTYPE=0 IC={self.critical_current:.17g} C={self.junction_capacitance:.17g} RN={self.junction_resistance:.17g} R0={self.junction_resistance:.17g})',
               'BJJ J 0 JJMOD',f'CSH C 0 {self.shunt_capacitance:.17g}',
               f'RINT C 0 {self.internal_resistance:.17g}',
               f'LTANK C J {self.interconnect_inductance:.17g}',
               f'LPICK PORT 0 {self.pickup_inductance:.17g}']
        if k:lines.append(f'KEM LTANK LPICK {k:.17g}')
        return lines


@dataclass(frozen=True)
class GalvanicJPA:
    """Existing bare JPA driven at its actual signal/ground reference patches.

    The 2x2 matrix comprises JJ-feed and capacitor-feed loops referred to the
    finite external electrodes, retaining shared-return mutual inductance.
    """
    critical_current: float
    shunt_capacitance: float
    junction_capacitance: float
    junction_feed_inductance: float
    capacitor_feed_inductance: float
    mutual_inductance: float
    port_impedance: float = 5000.
    internal_resistance: float = 1e8
    junction_resistance: float = 1e9

    def __post_init__(self):
        for name,value in asdict(self).items():
            if name!='mutual_inductance':positive(value,name)
        if not np.isfinite(self.mutual_inductance) or np.linalg.eigvalsh(self.inductance_matrix)[0]<=0:raise ValueError('Positive-definite feed matrix required')

    @property
    def inductance_matrix(self):return np.array([[self.junction_feed_inductance,self.mutual_inductance],[self.mutual_inductance,self.capacitor_feed_inductance]])

    @property
    def port_incidence(self):return np.ones(2)

    plasma_omega=CoupledJPA.plasma_omega
    modes=CoupledJPA.modes
    passive=CoupledJPA.passive
    energy=CoupledJPA.energy
    rhs_scaled=CoupledJPA.rhs_scaled

    def linear_system(self):
        w=self.plasma_omega;inv=np.linalg.inv(self.inductance_matrix);lj=PHIBAR/self.critical_current
        a=np.zeros((5,5));a[0,1]=1;a[1,0]=-1;a[1,1]=-1/(self.junction_resistance*self.junction_capacitance*w);a[1,3]=1
        a[2,2]=-1/(self.internal_resistance*self.shunt_capacitance*w);a[2,4]=self.junction_capacitance/self.shunt_capacitance
        a[3:,1]=-inv[:,0]*lj;a[3:,2]=-inv[:,1]*lj
        a[3:,3:]=-inv@np.ones((2,2))*self.port_impedance/w
        b=np.zeros(5);b[3:]=inv@np.ones(2)/(self.critical_current*w)
        return a,b

    def components(self):
        k=self.mutual_inductance/np.sqrt(self.junction_feed_inductance*self.capacitor_feed_inductance)
        lines=[f'.model JJMOD jj(RTYPE=0 IC={self.critical_current:.17g} C={self.junction_capacitance:.17g} RN={self.junction_resistance:.17g} R0={self.junction_resistance:.17g})',
               'BJJ J 0 JJMOD',f'CSH C 0 {self.shunt_capacitance:.17g}',f'RINT C 0 {self.internal_resistance:.17g}',
               f'LJFEED PORT J {self.junction_feed_inductance:.17g}',f'LCFEED PORT C {self.capacitor_feed_inductance:.17g}']
        if k:lines.append(f'KEM LJFEED LCFEED {k:.17g}')
        return lines


@dataclass
class PeriodicPump:
    circuit: CoupledJPA | GalvanicJPA
    frequency: float
    power: float
    harmonics: np.ndarray
    state: np.ndarray
    residual: float
    phase_peak: float
    floquet_radius: float | None = None


def pump(circuit,frequency,power,*,harmonics=3,initial=None,quadrature=512,stability=True):
    positive(frequency,'pump frequency');positive(power,'pump power',True)
    if (type(harmonics) is not int or harmonics<1 or type(quadrature) is not int or
        quadrature<max(64,8*harmonics)):raise ValueError('Invalid harmonic/quadrature controls')
    modes=np.arange(1,2*harmonics,2);theta=2*np.pi*np.arange(quadrature)/quadrature
    basis=np.exp(-1j*np.outer(theta,modes));project=2*basis.conj().T/quadrature
    a,b=circuit.linear_system();mat=np.array([-1j*2*np.pi*frequency*k/circuit.plasma_omega*np.eye(5)-a for k in modes])
    source=np.zeros((len(modes),5),complex);source[0]=b*1j*np.sqrt(8*circuit.port_impedance*power)
    def unpack(v):return (v[:len(v)//2]+1j*v[len(v)//2:]).reshape(len(modes),5)
    def pack(x):return np.r_[x.real.ravel(),x.imag.ravel()]
    if initial is None:guess=np.array([np.linalg.solve(m,s) for m,s in zip(mat,source)])
    else:
        guess=np.zeros_like(source)
        for k,h in enumerate(initial.harmonics):
            if h in modes:guess[np.where(modes==h)[0][0]]=initial.state[k]
    def residual(v):
        x=unpack(v);phase=(basis@x[:,0]).real
        result=np.einsum('kij,kj->ki',mat,x)-source
        result[:,1]+=project@(np.sin(phase)-phase)
        return pack(result)
    fit=root(residual,pack(guess),tol=1e-10)
    x=unpack(fit.x);err=float(np.max(abs(residual(fit.x))))
    if not np.isfinite(err) or not np.all(np.isfinite(x)) or err>1e-8:raise RuntimeError(f'Periodic pump did not converge: {err:g}')
    phase=float(np.max(abs((basis@x[:,0]).real)))
    if phase>1:raise ValueError('Junction phase exceeds the supported un-switched phase guard')
    result=PeriodicPump(circuit,frequency,power,modes,x,err,phase)
    if stability:
        omega=2*np.pi*frequency/circuit.plasma_omega
        def ode(t,y):
            ph=float(np.real(np.exp(-1j*modes*omega*t)@x[:,0]))
            jac=a.copy();jac[1,0]=-np.cos(ph)
            return (jac@y.reshape(5,5)).ravel()
        solution=solve_ivp(ode,(0,2*np.pi/omega),np.eye(5).ravel(),method='DOP853',rtol=1e-9,atol=1e-11)
        if not solution.success:raise RuntimeError('Floquet stability integration failed')
        result.floquet_radius=float(max(abs(np.linalg.eigvals(solution.y[:,-1].reshape(5,5)))))
        if result.floquet_radius>=1-1e-8:raise ValueError(f'Periodic pump is unstable/marginal: {result.floquet_radius:g}')
    return result


def response(point,frequencies,*,sideband_order=2,quadrature=1024):
    c=point.circuit;a,b=c.linear_system();fs=np.atleast_1d(np.asarray(frequencies,float))
    if fs.ndim!=1 or len(fs)==0 or np.any(~np.isfinite(fs)) or np.any(fs<=0) or np.any(2*point.frequency-fs<=0):raise ValueError('Positive finite signal and idler frequencies required')
    if type(sideband_order) is not int or sideband_order<1:raise ValueError('sideband_order must be positive integer')
    if type(quadrature) is not int or quadrature<max(64,8*len(point.harmonics),8*sideband_order):raise ValueError('Invalid Fourier quadrature')
    theta=2*np.pi*np.arange(quadrature)/quadrature
    phase=np.real(np.exp(-1j*np.outer(theta,point.harmonics))@point.state[:,0]);cosine=np.cos(phase)-1
    bands=np.arange(-sideband_order,sideband_order+1);n=len(bands);signal=[];idler=[]
    for f in fs:
        freqs=f+2*bands*point.frequency
        matrix=np.zeros((5*n,5*n),complex);rhs=np.zeros(5*n,complex)
        for i,band in enumerate(bands):
            matrix[5*i:5*i+5,5*i:5*i+5]=-1j*2*np.pi*freqs[i]/c.plasma_omega*np.eye(5)-a
            if band==0:rhs[5*i:5*i+5]=b*2
            for j,other in enumerate(bands):
                coeff=np.mean(cosine*np.exp(2j*(band-other)*theta))
                matrix[5*i+1,5*j]+=coeff
        x=np.linalg.solve(matrix,rhs).reshape(n,5)
        signal.append(1-c.port_impedance*c.critical_current*(c.port_incidence@x[bands==0,3:][0]))
        vi=np.conj(-c.port_impedance*c.critical_current*(c.port_incidence@x[bands==-1,3:][0]))
        idler.append(vi*np.sqrt(f/(2*point.frequency-f)))
    return dict(signal_frequency=fs.tolist(),signal=np.array(signal),idler=np.array(idler),
                gain=abs(np.array(signal))**2,idler_conversion=abs(np.array(idler))**2,
                sideband_order=sideband_order)
