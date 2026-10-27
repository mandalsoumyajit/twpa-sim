"""The four electrical parameter sets in the original MATLAB model."""
from .models import Cell, Section, Design, Bias


def preset(name='v1-target'):
    if name in ('v2-tls', 'v2-standard', 'v2-folded', 'v2-four-pass'):
        # sim_tls_and_spine.py: Cu/Cl are TOTAL cell capacitances.
        u = Cell(20e-12*(1.5/.58)*1.04, 22.32e-15/2, 20e-12*(75.5/.5)/3)
        l = Cell(u.inductance, 8.38e-15/2, 20e-12*(27.5/.5)/3)
        n = {'v2-tls': 925, 'v2-standard': 925, 'v2-folded': 768, 'v2-four-pass': 726}[name]
        return Design((Section(u, 24), Section(l, 6), Section(u, 24)), n, name=name), Bias(kerr_factor=2.)
    ld, cf, factor, width, lu, ll, repeats, tail = 56e-12, 11.20e-15, 1., .35, 72., 28., 914, 5
    if name in ('v1-cpw', 'v1-idc'):
        ld, cf, factor = 54.08e-12, (10.361e-15 if name == 'v1-cpw' else 8.77e-15), 1.04
    elif name == 'v2':
        ld, cf, width, lu, ll, repeats, tail = 54.45e-12, 11.78e-15, .50, 75.5, 27.5, 1050, 0
    elif name != 'v1-target':
        raise ValueError(f'Unknown preset: {name}')
    u = Cell(ld, cf, 20e-12*(lu/width)/3*factor)
    l = Cell(ld, cf*ll/lu, 20e-12*(ll/width)/3*factor)
    return Design((Section(u, 24), Section(l, 6), Section(u, 24)), repeats,
                  (Section(u, tail),), name), Bias()
