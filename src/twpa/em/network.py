"""Explicit extracted terminal networks, without uniform-cell reduction.

Node flux divided by PHIBAR is the coordinate. All metal star branches retain
their signed inductance matrix; JJ sine currents, individual capacitances and
losses enter nodal KCL. Applied loop flux is an explicit external assumption.
"""
from dataclasses import dataclass
from copy import deepcopy
import numpy as np
from scipy.optimize import root
from ..jpa_circuit import PHIBAR,PHI0


@dataclass
class LayoutNetwork:
    data: dict

    def __post_init__(self):
        self.data=deepcopy(self.data)
        d=self.data;self.nodes=list(d['nodes']);self.reference=d['reference']
        if len(set(self.nodes))!=len(self.nodes) or self.reference not in self.nodes:raise ValueError('Invalid network nodes')
        self.active=[n for n in self.nodes if n!=self.reference];self.index={n:i for i,n in enumerate(self.active)}
        self.n=len(self.active)
        if self.n>200:raise ValueError('Local network envelope is 200 active nodes')
        self.L=np.asarray(d['inductance_H'],float)
        self.BL=self.incidence(d['metal_branches']);self.BJ=self.incidence(d['junctions'])
        if self.L.shape!=(self.BL.shape[1],)*2 or not np.all(np.isfinite(self.L)):raise ValueError('Invalid inductance dimensions')
        if not np.allclose(self.L,self.L.T,rtol=1e-8,atol=1e-20) or np.linalg.eigvalsh(self.L)[0]<=0:raise ValueError('Inductance matrix must be reciprocal and positive definite')
        self.invL=np.linalg.inv(self.L);self.K=self.BL@self.invL@self.BL.T
        self.ic=np.array([j['critical_current'] for j in d['junctions']],float)
        self.scale=float(max(self.ic))
        self.C=np.zeros((self.n,self.n));self.G=np.zeros_like(self.C)
        for family in ('junctions','capacitors'):
            for item in d[family]:
                b=self.incidence([item])[:,0];c=float(item['capacitance']);r=float(item['resistance'])
                if not np.isfinite(c) or c<=0 or not np.isfinite(r) or r<=0:raise ValueError('Positive finite C/R required')
                self.C+=c*np.outer(b,b);self.G+=np.outer(b,b)/r
        if np.any(~np.isfinite(self.ic)) or np.any(self.ic<=0):raise ValueError('Positive finite junction currents required')
        self.BP=self.incidence(d['ports']);self.R=np.array([p['impedance'] for p in d['ports']],float)
        if np.any(~np.isfinite(self.R)) or np.any(self.R<=0):raise ValueError('Positive finite port resistances required')
        self.G+=(self.BP/self.R)@self.BP.T
        self.cycles=self.cycle_basis()
        self.CL=self.cycles[:,:self.BL.shape[1]]
        # DC phases have one arbitrary offset per disconnected superconducting
        # component. Remove those gauges only in the static equations.
        graph=np.hstack((self.BL,self.BJ));adj=abs(graph)@abs(graph).T
        remaining=set(range(self.n));free=set(remaining)
        grounded=set()
        for edge in d['metal_branches']+d['junctions']:
            if self.reference in (edge['positive'],edge['negative']):
                other=edge['negative'] if edge['positive']==self.reference else edge['positive']
                grounded.add(self.index[other])
        while remaining:
            group={next(iter(remaining))};todo=list(group)
            while todo:
                node=todo.pop()
                for other in np.flatnonzero(adj[node]):
                    if int(other) not in group:group.add(int(other));todo.append(int(other))
            if not group & grounded:free.remove(min(group))
            remaining-=group
        self.dc_free=np.array(sorted(free),int)

    def incidence(self,edges):
        b=np.zeros((len(self.active),len(edges)))
        for k,e in enumerate(edges):
            for key,sign in [('positive',1),('negative',-1)]:
                n=e[key]
                if n not in self.nodes:raise ValueError('Unknown terminal')
                if n!=self.reference:b[self.index[n],k]+=sign
        return b

    def cycle_basis(self):
        """Fundamental superconducting loops, oriented by the closing JJ."""
        edges=self.data['metal_branches']+self.data['junctions'];tree={n:[] for n in self.nodes};cycles=[]
        for k,e in enumerate(edges):
            a,b=e['positive'],e['negative'];todo=[a];previous={a:None}
            while todo and b not in previous:
                node=todo.pop()
                for other,index,orientation in tree[node]:
                    if other not in previous:previous[other]=(node,index,orientation);todo.append(other)
            if b not in previous:
                tree[a].append((b,k,1));tree[b].append((a,k,-1));continue
            row=np.zeros(len(edges));row[k]=1;node=b
            while node!=a:
                old,index,orientation=previous[node];row[index]=-orientation;node=old
            cycles.append(row)
        return np.array(cycles).reshape(-1,len(edges))

    def flux_vector(self,flux_bias=0.):
        """Bias is one signed applied flux/Phi0 per fundamental physical loop."""
        value=np.asarray(flux_bias,float)
        if len(self.CL)==0 and np.any(value):raise ValueError('No superconducting flux loop in this circuit')
        if value.ndim==0:value=np.full(len(self.CL),float(value))
        if value.shape!=(len(self.CL),) or np.any(~np.isfinite(value)):raise ValueError('One finite flux bias per loop required')
        if len(self.CL)==0:return np.zeros(len(self.L))
        return -self.CL.T@np.linalg.solve(self.CL@self.CL.T,value*PHI0)

    def equilibrium(self,flux_bias=0.,initial=None):
        flux=self.flux_vector(flux_bias);free=self.dc_free
        base=np.zeros(self.n) if initial is None else np.array(initial,copy=True)
        def residual(x):
            theta=base.copy();theta[free]=x
            return (PHIBAR*self.K@theta-self.BL@self.invL@flux+self.BJ@(self.ic*np.sin(self.BJ.T@theta)))[free]/self.scale
        fit=root(residual,base[free],tol=1e-10)
        base[free]=fit.x
        if not np.all(np.isfinite(base)) or max(abs(residual(fit.x)),default=0)>1e-8:raise RuntimeError('Static flux equilibrium failed')
        hessian=PHIBAR*self.K+(self.BJ*(self.ic*np.cos(self.BJ.T@base)))@self.BJ.T
        if np.linalg.eigvalsh(hessian[np.ix_(free,free)])[0]<=0:raise ValueError('Static equilibrium is unstable/marginal')
        return base

    def linear_matrix(self,frequency,theta):
        w=2*np.pi*frequency
        return PHIBAR*(self.K-w*w*self.C-1j*w*self.G)+(self.BJ*(self.ic*np.cos(self.BJ.T@theta)))@self.BJ.T

    def passive(self,frequencies,*,flux_bias=0.,input_port=0):
        theta=self.equilibrium(flux_bias);results=[]
        for f in np.atleast_1d(frequencies):
            if not np.isfinite(f) or f<=0:raise ValueError('Positive finite frequencies required')
            rhs=2*self.BP[:,input_port]/np.sqrt(self.R[input_port])
            x=np.linalg.solve(self.linear_matrix(f,theta),rhs)
            v=-1j*2*np.pi*f*PHIBAR*(self.BP.T@x)
            s=v/np.sqrt(self.R);s[input_port]-=1;results.append(s)
        return np.array(results)

    def periodic(self,frequency,power,*,flux_bias=0.,flux_pump=0.,harmonics=4,quadrature=256,initial=None):
        if not np.isfinite(frequency) or frequency<=0 or not np.isfinite(power) or power<0:raise ValueError('Invalid drive')
        if type(harmonics) is not int or harmonics<1 or type(quadrature) is not int or quadrature<max(64,8*harmonics):raise ValueError('Invalid spectral controls')
        static=self.equilibrium(flux_bias);flux=self.flux_vector(flux_bias);free=self.dc_free
        hs=np.arange(1,harmonics+1);time=2*np.pi*np.arange(quadrature)/quadrature
        basis=np.exp(-1j*np.outer(time,hs));project=2*basis.conj().T/quadrature
        matrices=np.array([PHIBAR*(self.K-(2*np.pi*frequency*h)**2*self.C-1j*2*np.pi*frequency*h*self.G) for h in hs])
        source=np.zeros((harmonics,self.n),complex);source[0]=self.BP[:,0]*1j*np.sqrt(8*self.R[0]*power)/self.R[0]
        source[0]+=1j*self.BL@self.invL@self.flux_vector(flux_pump)
        linear=np.array([np.linalg.solve(self.linear_matrix(h*frequency,static),s) for h,s in zip(hs,source)])
        def pack(dc,x):return np.r_[dc[free],x.real.ravel(),x.imag.ravel()]
        def unpack(v):
            dc=static.copy();dc[free]=v[:len(free)];q=v[len(free):];n=harmonics*self.n
            return dc,(q[:n]+1j*q[n:]).reshape(harmonics,self.n)
        def residual(v):
            dc,x=unpack(v);delta=self.BJ.T@dc+(basis@(x@self.BJ)).real
            current=np.sin(delta)*self.ic
            dc_res=(PHIBAR*self.K@dc-self.BL@self.invL@flux+self.BJ@np.mean(current,axis=0))/self.scale
            ac=(np.einsum('hij,hj->hi',matrices,x)+(project@current)@self.BJ.T-source)/self.scale
            return pack(dc_res,ac)
        guess=pack(static,linear) if initial is None else pack(initial['dc'],initial['state'])
        fit=root(residual,guess,tol=1e-9);dc,x=unpack(fit.x);err=float(max(abs(residual(fit.x))))
        if not np.isfinite(err) or err>1e-7:raise RuntimeError(f'Periodic network residual {err:g}')
        excursion=float(np.max(np.sum(abs(x@self.BJ),axis=0)))
        if excursion>1.:raise ValueError('AC junction phase guard exceeded')
        return dict(frequency=frequency,power=power,flux_bias=np.asarray(flux_bias).tolist(),flux_pump=np.asarray(flux_pump).tolist(),dc=dc,state=x,
                    harmonics=harmonics,quadrature=quadrature,residual=err,maximum_phase=excursion)

    def response(self,point,frequencies,*,sidebands=2,quadrature=512):
        if type(sidebands) is not int or sidebands<1:raise ValueError('Positive sideband count required')
        if type(quadrature) is not int or quadrature<max(64,8*(point['harmonics']+sidebands)):
            raise ValueError('Insufficient sideband quadrature')
        fp=point['frequency'];theta=2*np.pi*np.arange(quadrature)/quadrature
        delta=self.BJ.T@point['dc']+(np.exp(-1j*np.outer(theta,np.arange(1,point['harmonics']+1)))@(point['state']@self.BJ)).real
        bands=np.arange(-sidebands,sidebands+1);coefficient={k:np.mean(np.cos(delta)*np.exp(1j*k*theta)[:,None],axis=0)*self.ic for k in range(-2*sidebands,2*sidebands+1)}
        out=[]
        for f in np.atleast_1d(frequencies):
            fs=f+bands*fp
            if np.any(abs(fs)<1.) or not np.isfinite(f) or f<=0 or f>=2*fp:raise ValueError('Supported distinct positive signal/idler frequencies required')
            a=np.zeros((self.n*len(bands),)*2,complex);rhs=np.zeros(self.n*len(bands),complex)
            for i,k in enumerate(bands):
                s=slice(i*self.n,(i+1)*self.n);w=2*np.pi*fs[i]
                a[s,s]=PHIBAR*(self.K-w*w*self.C-1j*w*self.G)
                if k==0:rhs[s]=2*self.BP[:,0]/np.sqrt(self.R[0])
                for j,l in enumerate(bands):a[s,j*self.n:(j+1)*self.n]+=(self.BJ*coefficient[int(k-l)])@self.BJ.T
            x=np.linalg.solve(a,rhs).reshape(len(bands),self.n)
            voltage=-1j*2*np.pi*fs[:,None]*PHIBAR*(x@self.BP)
            waves=voltage/np.sqrt(self.R);waves[bands==0,0]-=1
            out.append(waves)
        return dict(frequency=np.atleast_1d(frequencies),bands=bands,waves=np.array(out))

    def netlist(self,tones,transient,*,flux_bias=0.,flux_pump=0.,flux_frequency=None,bias_ramp=20e-9):
        """Physical RCSJ/linear-metal deck; ideal external flux coils are explicit.

        External coils provide declared loop flux, not an extracted coil mask.
        Their mutuals and self L form a positive-definite augmented matrix.
        """
        from ..jpa_circuit import voltage_sources
        def node(name):return '0' if name==self.reference else 'T'+str(self.index[name])
        lines=['* twpa-sim geometry-derived finite terminal network']+voltage_sources(tones,float(self.R[0]))
        for i,p in enumerate(self.data['ports']):
            source='SRC0' if i==0 else node(p['negative'])
            lines.append(f'RPORT{i} {source} {node(p["positive"])} {self.R[i]:.17g}')
        for i,j in enumerate(self.data['junctions']):
            lines += [f'.model JJ{i} jj(RTYPE=0 IC={j["critical_current"]:.17g} C={j["capacitance"]:.17g} RN={j["resistance"]:.17g} R0={j["resistance"]:.17g})',
                      f'BJJ{i} {node(j["positive"])} {node(j["negative"])} JJ{i}']
        for i,c in enumerate(self.data['capacitors']):
            lines += [f'CC{i} {node(c["positive"])} {node(c["negative"])} {c["capacitance"]:.17g}',
                      f'RC{i} {node(c["positive"])} {node(c["negative"])} {c["resistance"]:.17g}']
        for i,b in enumerate(self.data['metal_branches']):lines.append(f'LM{i} {node(b["positive"])} {node(b["negative"])} {self.L[i,i]:.17g}')
        for i in range(len(self.L)):
            for j in range(i):
                k=self.L[i,j]/np.sqrt(self.L[i,i]*self.L[j,j])
                if k:lines.append(f'KM{i}_{j} LM{i} LM{j} {k:.17g}')
        nl=len(self.CL)
        if nl:
            amplitudes=[];dc=np.broadcast_to(np.asarray(flux_bias,float),(nl,));ac=np.broadcast_to(np.asarray(flux_pump,float),(nl,))
            if np.any(ac) and (flux_frequency is None or flux_frequency<=0):raise ValueError('Positive physical flux pump frequency required')
            # One milliampere per Phi0 is an external source normalization.
            basis=np.array([self.flux_vector(np.eye(nl)[i]) for i in range(nl)]).T/1e-3
            coil=max(1e-9,2*float(np.linalg.eigvalsh(basis.T@self.invL@basis)[-1]))
            for j in range(nl):
                grid=np.linspace(0,1,129);shape=10*grid**3-15*grid**4+6*grid**5
                points=' '.join(f'{t*bias_ramp:.17g} {v*dc[j]*1e-3:.17g}' for t,v in zip(grid,shape))
                lines += [f'IFLUX{j} 0 COIL{j} PWL({points})',f'LBIAS{j} COIL{j} 0 {coil:.17g}']
                if ac[j]:
                    # Parallel prescribed AC current; the inductor integrates
                    # exactly the same flux vector as Python's forcing term.
                    lines += [f'IPUMP{j} 0 COIL{j} SIN(0 {ac[j]*1e-3:.17g} {flux_frequency:.17g} 0 0 0)']
                for i in range(len(self.L)):
                    k=basis[i,j]/np.sqrt(self.L[i,i]*coil)
                    if k:lines.append(f'KBIAS{i}_{j} LM{i} LBIAS{j} {k:.17g}')
        else:
            if np.any(np.asarray(flux_bias)) or np.any(np.asarray(flux_pump)):raise ValueError('No physical flux loop in this layout')
        traces=['V(SRC0)']+[f'V({node(n)})' for n in self.active]+[f'P(BJJ{i})' for i in range(len(self.ic))]
        lines += [f'.tran {transient.step:.17g} {transient.stop:.17g} {transient.start:.17g} {transient.output_step:.17g}',
                  '.print '+' '.join(traces),'.end']
        return '\n'.join(lines)+'\n'
