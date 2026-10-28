"""Isolated openEMS Python runtime backend with immutable inputs and logs."""
from dataclasses import dataclass
from datetime import datetime,timezone
from importlib.resources import files
import hashlib
import json
from pathlib import Path
import subprocess
import time

from ..layout.generate import digest
from ..layout.process import positive


class ExtractionError(RuntimeError):
    pass


@dataclass(frozen=True)
class OpenEMS:
    python_executable: str
    dll_directory: str | None = None
    timeout: float = 600.
    threads: int = 1

    def __post_init__(self):
        if not self.python_executable:
            raise ValueError('An openEMS-capable Python executable is required')
        positive(self.timeout,'timeout')
        if isinstance(self.threads,bool) or not isinstance(self.threads,int) or self.threads<1:
            raise ValueError('threads must be a positive integer')

    def command(self,worker):
        args=[str(self.python_executable),'-u',str(worker),'--threads',str(self.threads)]
        if self.dll_directory:
            args+=['--dll-directory',str(self.dll_directory)]
        return args

    def run(self,job,directory):
        root=Path(directory).resolve(); root.mkdir(parents=True,exist_ok=False)
        text=json.dumps(job,indent=2,allow_nan=False)+'\n'
        worker=files('twpa.em').joinpath('openems_worker.py').read_bytes()
        (root/'job.json').write_text(text,encoding='utf-8')
        (root/'worker.py').write_bytes(worker)
        command=self.command(root/'worker.py')+['--input',str(root/'job.json')]
        manifest=dict(started_utc=datetime.now(timezone.utc).isoformat(),command=command,
                      job_sha256=digest(job),worker_sha256=hashlib.sha256(worker).hexdigest())
        started=time.monotonic()
        try:
            result=subprocess.run(command,cwd=root,capture_output=True,text=True,timeout=self.timeout)
            code=result.returncode; stdout=result.stdout; stderr=result.stderr
        except subprocess.TimeoutExpired as exc:
            code=None
            stdout=exc.stdout.decode(errors='replace') if isinstance(exc.stdout,bytes) else exc.stdout or ''
            stderr=exc.stderr.decode(errors='replace') if isinstance(exc.stderr,bytes) else exc.stderr or ''
        except OSError as exc:
            code=None; stdout=''; stderr=str(exc)
        manifest.update(returncode=code,elapsed_seconds=time.monotonic()-started)
        (root/'run.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        (root/'stdout.log').write_text(stdout,encoding='utf-8')
        (root/'stderr.log').write_text(stderr,encoding='utf-8')
        if code!=0 or not (root/'extraction.json').is_file():
            raise ExtractionError(f'openEMS failed/timed out ({code}); see {root}')
        data=json.loads((root/'extraction.json').read_text(encoding='utf-8'))
        data['run']=manifest
        if not data['fit_accepted']:
            raise ExtractionError(f'Extraction failed numerical gates; see {root}/extraction.json')
        data['artifact_sha256']=digest(data)
        (root/'extraction.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        return data
