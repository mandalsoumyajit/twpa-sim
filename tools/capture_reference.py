"""Capture regression data from user-supplied ORIGINAL solver files.

Usage: python tools/capture_reference.py PATH/TO/NY_CREATES/scripts
Only numerical definitions and the TLS sweep are executed, excluding all
plotting, environment mutations, GDS work and hardcoded output paths.
"""
import ast
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def execute(nodes, namespace):
    module = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    exec(compile(module, '<original numerical solver>', 'exec'), namespace)


def capture(source):
    tls_file = source/'sim_tls_and_spine.py'
    tree = ast.parse(tls_file.read_text(encoding='utf-8-sig'))
    nodes = []
    for node in tree.body:
        # Stop before the turn-spine analysis. Its output is unrelated to TLS.
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'R_turn' for t in node.targets):
            break
        if isinstance(node, ast.For) and isinstance(node.target, ast.Name) and node.target.id == 'td':
            collect = ast.parse("rows.append(dict(loss=td, gain=float(G_net), noise=float(N_add), "
                                "kp=float(k_p), ks=float(ks), ki=float(ki), "
                                "alpha_s=float(gam_s_exact), alpha_i=float(gam_i_exact)))").body
            node.body.extend(collect)
        nodes.append(node)
    ns = {'rows': []}
    with contextlib.redirect_stdout(io.StringIO()):
        execute(nodes, ns)
    ghz = np.array([1., 3., 4.6, 6., 9.5, 10.5])
    bloch = []
    for td in (0., 1.5e-4, 1e-3):
        k, alpha, matrix, _ = ns['get_bloch_dispersion'](ghz, tan_d_eff=td)
        bloch.append(dict(loss=td, frequency=ghz.tolist(), k=k.tolist(), alpha=alpha.tolist()))
    laplace = source/'generate_reorganized_figures.py'
    funcs = [n for n in ast.parse(laplace.read_text(encoding='utf-8-sig')).body
             if isinstance(n, ast.FunctionDef) and n.name == 'compute_2d_laplace_fields']
    field_ns = dict(np=np, lil_matrix=lil_matrix, spsolve=spsolve,
                    eps0=ns['eps0'], eps_si=ns['eps_si'])
    execute(funcs, field_ns)
    x, y, v, xc, yc, energy = field_ns['compute_2d_laplace_fields']()
    xx, yy = np.meshgrid(xc*1e-9, yc*1e-9, indexing='ij')
    oxide = (((yy >= 0) & (yy <= 60e-9) & (xx >= 250e-9)) |
             ((yy >= 60e-9) & (yy <= 120e-9) & (xx <= 250e-9)) |
             ((yy >= 0) & (yy <= 120e-9) & (xx >= 250e-9) & (xx <= 310e-9)))
    probes = [(0, 0), (100, 80), (120, 100), (149, 159), (60, 40)]
    idc = dict(energy=float(np.nansum(energy)*2.5e-9**2),
               participation=float(np.nansum(energy[oxide])/np.nansum(energy)),
               potential=[float(v[i, j]) for i, j in probes], probes=probes)
    nfile = source/'sim_alvin_twpa.py'
    funcs = [n for n in ast.parse(nfile.read_text(encoding='utf-8-sig')).body
             if isinstance(n, ast.FunctionDef) and n.name == 'cme_rhs']
    ode_ns = {'np': np}; execute(funcs, ode_ns)
    kv = [ns['rows'][0][key] for key in ('kp', 'ks', 'ki')]
    epsilon = [ns['eps1'], ns['eps2']]
    initial = np.array([ns['I_p0'], ns['I_p0']/100, 0.], complex)
    sol = solve_ivp(ode_ns['cme_rhs'], (0, ns['N_c']), initial,
                    args=(epsilon, kv), rtol=1e-10, atol=1e-14)
    nonlinear = dict(k=kv, epsilon=epsilon, cells=ns['N_c'],
                     initial_real=initial.real.tolist(), output_real=sol.y[:, -1].real.tolist(),
                     output_imag=sol.y[:, -1].imag.tolist())
    files = ['sim_tls_and_spine.py', 'generate_reorganized_figures.py',
             'sim_alvin_twpa.py', 'sim_detailed_physics.py', 'sim_gradients_exact.py']
    return dict(source_sha256={f: hashlib.sha256((source/f).read_bytes()).hexdigest() for f in files},
                pump_frequency=ns['f_p_ghz']*1e9, coupled=ns['rows'], bloch=bloch,
                idc=idc, nonlinear=nonlinear)


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    result = capture(Path(sys.argv[1]))
    (root/'tests/fixtures/original_solvers.json').write_text(json.dumps(result, indent=2)+'\n')
    print('Captured original coupled, Bloch, IDC and nonlinear regression data.')
