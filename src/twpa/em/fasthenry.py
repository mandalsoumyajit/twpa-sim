"""Whiteley superconducting FastHenry runner and strict signed-matrix parsing."""
from dataclasses import dataclass
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import numpy as np

from ..backends.josim import wsl_path
from ..layout.generate import digest
from ..layout.process import positive


def parse_impedance(text,ports):
    mapping={int(i):(a.casefold(),b.casefold()) for i,a,b in re.findall(r'Row\s+(\d+):\s+(\S+)\s+to\s+(\S+)',text,re.I)}
    expected={i+1:(p['positive'].casefold(),p['negative'].casefold()) for i,p in enumerate(ports)}
    if mapping!=expected:raise ValueError('FastHenry port order/orientation mismatch')
    headers=list(re.finditer(r'Impedance matrix for frequency\s*=\s*([\deE+.-]+)\s+(\d+)\s+x\s+(\d+)',text))
    if not headers or 'ADMITTANCE matrix' in text:raise ValueError('Missing impedance matrices')
    samples=[]
    for k,m in enumerate(headers):
        f=float(m[1]);n=int(m[2]);cols=int(m[3])
        words=text[m.end():headers[k+1].start() if k+1<len(headers) else len(text)].split()
        if n!=len(ports) or cols!=n or len(words)!=2*n*n or f<=0:raise ValueError('Incomplete/invalid impedance matrix')
        z=np.array([complex(float(words[2*i]),float(words[2*i+1].rstrip('j'))) for i in range(n*n)]).reshape(n,n)
        if not np.all(np.isfinite(z)):raise ValueError('Nonfinite impedance')
        samples.append((f,z))
    if np.any(np.diff([f for f,z in samples])<=0):raise ValueError('Frequency samples not increasing')
    return samples


@dataclass(frozen=True)
class FastHenry:
    executable: str
    wsl_distribution: str | None = None
    binary_sha256: str | None = None
    timeout: float = 180.
    memory_GiB: float = 4.
    solver: str = 'auto'
    threads: int = 2

    @classmethod
    def from_installation(cls,path=None,**overrides):
        if path is None:
            if os.name=='nt':path=Path(os.environ['LOCALAPPDATA'])/'SeeQC/tools/fasthenry-wr-oriented-field.json'
            else:path=Path('/usr/local/share/seeqc/fasthenry-wr-oriented-field.json')
        data=json.loads(Path(path).read_text())
        if data.get('variant')!='Whiteley superconducting FastHenry':raise ValueError('A London-aware Whiteley installation is required')
        if os.name=='nt':
            return cls(data['linux_binary'],data['wsl_distro'],data['binary_sha256'],**overrides)
        return cls(data.get('linux_binary') or data['binary_path'],None,data['binary_sha256'],**overrides)

    def _prefix(self):
        return ['wsl','-d',self.wsl_distribution,'--exec'] if self.wsl_distribution else []

    def probe(self):
        if self.wsl_distribution:
            check=subprocess.run(self._prefix()+['sha256sum',self.executable],capture_output=True,text=True,check=True,timeout=20)
            binary_hash=check.stdout.split()[0]
        else:binary_hash=hashlib.sha256(Path(self.executable).read_bytes()).hexdigest()
        if self.binary_sha256 and binary_hash!=self.binary_sha256:raise ValueError('FastHenry executable hash differs from selected installation')
        run=subprocess.run(self._prefix()+[self.executable,'-h'],capture_output=True,text=True,timeout=20)
        banner=run.stdout+run.stderr
        if '3.0wr' not in banner:raise ValueError('London support requires the Whiteley FastHenry build')
        return dict(binary_sha256=binary_hash,version=next(s.strip() for s in banner.splitlines() if 'Version' in s),
                    variant='Whiteley superconducting FastHenry',command=self.executable)

    def run(self,job,directory):
        positive(self.timeout,'timeout');positive(self.memory_GiB,'memory_GiB')
        if type(self.threads) is not int or not 1<=self.threads<=4:raise ValueError('FastHenry threads must be 1..4')
        info=self.probe();root=Path(directory).resolve();root.mkdir(parents=True,exist_ok=False)
        (root/'input.inp').write_text(job['deck'],encoding='utf-8',newline='\n')
        (root/'job.json').write_text(json.dumps(job,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        # Native GNU timeout and address-space limit remain active if the WSL
        # launcher exits. Arguments are passed without shell interpolation.
        solver=('ludecomp' if len(job['settings']['frequencies_hz'])==1 else 'iterative') if self.solver=='auto' else self.solver
        if solver not in ('ludecomp','iterative'):raise ValueError('Unsupported FastHenry solver')
        args=[self.executable,'input.inp','-a','off','-s',solver,'-m','direct','-p','diag','-t','1e-9','-b','1e-9','-c','1000']
        args=['env',f'OPENBLAS_NUM_THREADS={self.threads}',f'OMP_NUM_THREADS={self.threads}',f'MKL_NUM_THREADS={self.threads}',
              f'SEEQC_FASTHENRY_BLAS_THREADS={self.threads}','stdbuf','-oL','-eL',*args]
        if self.wsl_distribution:
            cmd=['wsl','-d',self.wsl_distribution,'--cd',wsl_path(root),'--exec','timeout','--signal=TERM','--kill-after=5s',str(self.timeout),
                 'prlimit',f'--as={int(self.memory_GiB*2**30)}','--core=0','--',*args]
        elif os.name!='nt':
            cmd=['timeout','--signal=TERM','--kill-after=5s',str(self.timeout),'prlimit',f'--as={int(self.memory_GiB*2**30)}','--core=0','--',*args]
        else:raise ValueError('Use the supported WSL or native Linux FastHenry runtime')
        started=time.monotonic();record=dict(started_utc=datetime.now(timezone.utc).isoformat(),command=cmd,simulator=info,job_sha256=digest(job))
        try:
            r=subprocess.run(cmd,cwd=root,capture_output=True,text=True,timeout=self.timeout+15)
            code=r.returncode;stdout=r.stdout;stderr=r.stderr
        except subprocess.TimeoutExpired as e:
            code=None;stdout=e.stdout or b'';stderr=e.stderr or b''
            if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
            if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
        record.update(returncode=code,elapsed_seconds=time.monotonic()-started,wall_limit_seconds=self.timeout,address_space_limit_GiB=self.memory_GiB)
        for name,value in [('stdout.log',stdout),('stderr.log',stderr)]: (root/name).write_text(value,encoding='utf-8')
        (root/'run.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        if code!=0 or not (root/'Zc.mat').is_file():raise RuntimeError(f'FastHenry failed ({code}); see {root}')
        if any(s in (stdout+stderr).lower() for s in ('did not converge','failed to converge','maximum number of iterations reached')):
            raise RuntimeError(f'FastHenry iteration convergence failed; see {root}')
        samples=parse_impedance((root/'Zc.mat').read_text(),job['ports']);t=np.asarray(job['reduction'])
        output=[]
        for frequency,z in samples:
            reduced=t.T@z@t;l=reduced.imag/(2*np.pi*frequency)
            reciprocity=float(np.max(abs(z-z.T))/max(np.max(abs(z)),1e-30))
            eigen=float(np.min(np.linalg.eigvalsh((l+l.T)/2)))
            resistance=float(np.max(abs(reduced.real)))
            output.append(dict(frequency_Hz=frequency,L_H=l.tolist(),R_ohm=reduced.real.tolist(),
                               native_Z_real_ohm=z.real.tolist(),native_Z_imag_ohm=z.imag.tolist(),
                               reciprocity_relative=reciprocity,minimum_L_eigenvalue_H=eigen))
            if reciprocity>1e-6 or eigen<=0 or resistance>1e-6:raise ValueError('FastHenry reciprocity/positive-energy/lossless gate failed')
        if not np.allclose([f for f,z in samples],job['settings']['frequencies_hz'],rtol=1e-8):raise ValueError('Unexpected frequency sweep')
        result=dict(schema_version=1,kind=job['kind'],source=job['source'],manifest_sha256=job['manifest_sha256'],
                    stack=job['stack'],settings=job['settings'],ports=job['ports'],reduction=job['reduction'],
                    patches=job['patches'],audit=job['audit'],samples=output,run=record,scope=job['scope'])
        result['topology']=job.get('topology','pickup' if len(job['ports'])==3 else 'galvanic')
        result['artifact_sha256']=digest(result)
        (root/'extraction.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        return result
