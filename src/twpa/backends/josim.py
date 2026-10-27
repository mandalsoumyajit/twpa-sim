"""JoSIM CLI runner, with explicit native or Windows/WSL configuration."""
from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path, PureWindowsPath
import subprocess
import time
from datetime import datetime, timezone
from functools import lru_cache
import numpy as np


class JoSIMError(RuntimeError):
    pass


def wsl_path(path):
    p = PureWindowsPath(str(path))
    if not p.is_absolute() or len(p.drive) != 2 or p.drive[1] != ':':
        raise ValueError('WSL run directories must be absolute Windows drive paths')
    return '/mnt/'+p.drive[0].lower()+'/'+('/'.join(p.parts[1:]))


def load_csv(path):
    """Read upstream-compatible CSV without mangling JoSIM trace names."""
    try:
        with Path(path).open(newline='', encoding='utf-8-sig') as file:
            header = next(csv.reader(file))
    except (OSError, StopIteration) as exc:
        raise JoSIMError(f'Cannot read simulator CSV header: {path}') from exc
    names = [name.strip().strip('"') for name in header]
    if not names or names[0].casefold() != 'time' or len(set(n.casefold() for n in names)) != len(names):
        raise JoSIMError('missing or duplicate output trace names')
    try:
        data = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    except (ValueError, OSError) as exc:
        raise JoSIMError(f'Cannot parse simulator CSV data: {path}') from exc
    if data.shape[1] != len(names) or len(data) < 2 or not np.all(np.isfinite(data)):
        raise JoSIMError('invalid simulator output dimensions or values')
    if np.any(np.diff(data[:, 0]) <= 0):
        raise JoSIMError('simulator output time must strictly increase')
    return dict(zip(names, data.T))


@dataclass
class SimulationResult:
    directory: Path
    traces: dict
    manifest: dict

    @property
    def time(self):
        return next(iter(self.traces.values()))

    def trace(self, name):
        for key, value in self.traces.items():
            if key.casefold() == name.casefold():
                return value
        raise KeyError(f'Trace {name!r} unavailable; traces: {list(self.traces)}')


@dataclass(frozen=True)
class JoSIM:
    executable: str = 'josim-cli'
    wsl_distribution: str | None = None
    timeout: float = 180.
    corrector: int = 0  # Optional enhanced-fork feature; checked before use.
    source_revision: str | None = None
    analysis: int = 1  # JoSIM's default phase formulation preserves loop-flux constraints.

    def __post_init__(self):
        if not self.executable or not np.isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError('executable and positive finite timeout are required')
        if not isinstance(self.corrector, int) or self.corrector < 0:
            raise ValueError('corrector must be a nonnegative integer')
        if self.analysis not in (0, 1):
            raise ValueError('analysis must be 0 (voltage) or 1 (phase)')

    def _command(self, directory=None):
        if self.wsl_distribution:
            command = ['wsl', '-d', self.wsl_distribution]
            if directory is not None:
                command += ['--cd', wsl_path(directory)]
            return command+['--exec', self.executable]
        return [self.executable]

    @lru_cache(maxsize=16)
    def probe(self):
        try:
            version = subprocess.run(self._command()+['-v'], capture_output=True,
                                     text=True, timeout=min(self.timeout, 30), check=True)
            help_run = subprocess.run(self._command()+['-h'], capture_output=True,
                                      text=True, timeout=min(self.timeout, 30), check=True)
        except (OSError, subprocess.SubprocessError) as exc:
            raise JoSIMError(f'Cannot execute JoSIM: {exc}') from exc
        help_text = help_run.stdout+help_run.stderr
        return dict(version=(version.stdout+version.stderr).strip(),
                    corrector_supported='--corrector' in help_text,
                    source_revision=self.source_revision)

    def run(self, netlist, directory, *, metadata=None):
        """Execute a deck in a new directory; preserve inputs, logs and manifest.

        Native and WSL commands use argument arrays, never shell strings. The
        deck must be self-contained; referenced waveform/model files should be
        embedded by the caller or require a separate backend extension.
        """
        info = self.probe()
        if self.corrector and not info['corrector_supported']:
            raise JoSIMError('This JoSIM does not advertise the requested corrector feature')
        root = Path(directory).resolve()
        root.mkdir(parents=True, exist_ok=False)
        (root/'input.cir').write_text(netlist, encoding='utf-8', newline='\n')
        command = self._command(root)+['-a', str(self.analysis), '-m', '-o', 'output.csv']
        # GNU timeout terminates the Linux simulator even if the Windows WSL
        # launcher does not propagate process termination to its child.
        if self.wsl_distribution:
            prefix = self._command(root)
            command = prefix[:-1]+['timeout', '--signal=TERM', '--kill-after=5s',
                                   str(self.timeout), prefix[-1]]+command[len(prefix):]
        if self.corrector:
            command += ['-c', str(self.corrector)]
        command += ['input.cir']
        manifest = dict(simulator=info, command=command,
                        started_utc=datetime.now(timezone.utc).isoformat(),
                        netlist_sha256=hashlib.sha256(netlist.encode()).hexdigest(),
                        metadata={} if metadata is None else metadata)
        started = time.monotonic()
        try:
            result = subprocess.run(command, cwd=root, capture_output=True, text=True,
                                    timeout=self.timeout+(10 if self.wsl_distribution else 0), check=False)
            stdout, stderr, code = result.stdout, result.stderr, result.returncode
        except subprocess.TimeoutExpired as exc:
            def decoded(value):
                return value.decode(errors='replace') if isinstance(value, bytes) else (value or '')
            stdout, stderr, code = decoded(exc.stdout), decoded(exc.stderr), None
        except OSError as exc:
            stdout, stderr, code = '', str(exc), None
        manifest.update(returncode=code, elapsed_seconds=time.monotonic()-started)
        (root/'stdout.log').write_text(stdout, encoding='utf-8')
        (root/'stderr.log').write_text(stderr, encoding='utf-8')
        (root/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        if code != 0:
            raise JoSIMError(f'JoSIM failed or timed out (returncode={code}); see {root}')
        if not (root/'output.csv').is_file():
            raise JoSIMError(f'JoSIM produced no output.csv; see {root}')
        return SimulationResult(root, load_csv(root/'output.csv'), manifest)
