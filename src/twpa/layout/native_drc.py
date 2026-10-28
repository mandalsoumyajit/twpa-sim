"""Run the independent Ruby DRC deck in a full KLayout installation."""
from importlib.resources import as_file, files
import json
import os
from pathlib import Path
import subprocess

from .process import Process


def run_klayout_drc(gds, process, report, executable='klayout', timeout=300):
    """Return violation count from the actual .lyrdb, not the process exit code.

    KLayout normally exits zero even when a deck finds DRC violations.
    ``report`` must be a .lyrdb path. A process snapshot is saved beside it.
    """
    import klayout.rdb as rdb
    if not isinstance(process,Process):
        raise TypeError('A validated Process is required')
    report=Path(report).resolve()
    if report.suffix!='.lyrdb':
        raise ValueError('A .lyrdb output is required')
    report.parent.mkdir(parents=True,exist_ok=True)
    process_path=report.with_suffix('.process.json')
    process_path.write_text(json.dumps(process.to_dict(),indent=2)+'\n',encoding='utf-8')
    # Remove only a stale report at this exact requested output location.
    report.unlink(missing_ok=True)
    env=os.environ.copy(); env.setdefault('QT_QPA_PLATFORM','offscreen')
    with as_file(files('twpa.layout').joinpath('decks','process.drc')) as deck:
        run=subprocess.run([str(executable),'-b','-r',str(deck),'-rd',f'input={Path(gds).resolve()}',
                            '-rd',f'process={process_path}','-rd',f'report={report}'],
                           capture_output=True,text=True,timeout=timeout,env=env)
    if run.returncode or not report.is_file():
        raise RuntimeError(f'KLayout DRC failed ({run.returncode}): {run.stdout}\n{run.stderr}')
    database=rdb.ReportDatabase(); database.load(str(report))
    return dict(passed=database.num_items()==0, violations=database.num_items(),
                report=str(report),stdout=run.stdout,stderr=run.stderr)
