"""Versioned fabrication contracts. All layout lengths are in micrometres.

No KLayout dependency is needed to read or validate a process profile.
Polygons describe retained material/openings, not photoresist mask polarity.
"""
from copy import deepcopy
from importlib.resources import files
import json
import math
from pathlib import Path


def _keys(value, required, optional=()):
    if not isinstance(value, dict) or set(value) - set(required) - set(optional) or set(required) - set(value):
        raise ValueError(f"Expected keys {required}; optional {optional}")


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(text, object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {x}")))


class Process:
    """Validated profile; :meth:`to_dict` returns an independent copy."""

    def __init__(self, data):
        data = deepcopy(data)
        _keys(data, ('schema_version', 'id', 'revision', 'status', 'description', 'units',
                     'dbu_um', 'grid_um', 'layers', 'rules', 'sources', 'parameters', 'limitations'))
        integer(data['schema_version'],'schema_version')
        if data['schema_version'] != 1 or data['units'] != 'um':
            raise ValueError('Only schema_version 1 and units um are supported')
        if data['status'] not in ('research_reference', 'illustrative', 'local_unqualified'):
            raise ValueError('Unsupported qualification status; these decks do not provide signoff')
        for key in ('id', 'revision', 'description'):
            if not isinstance(data[key], str) or not data[key].strip():
                raise ValueError(f'Missing {key}')
        positive(data['dbu_um'], 'dbu_um'); positive(data['grid_um'], 'grid_um')
        ratio = data['grid_um'] / data['dbu_um']
        if not math.isclose(ratio, round(ratio), abs_tol=1e-9, rel_tol=0):
            raise ValueError('Manufacturing grid must be an integer multiple of DBU')
        sources = data['sources']
        if not isinstance(sources, dict) or not sources or any(not isinstance(v, str) or not v for v in sources.values()):
            raise ValueError('Nonempty source descriptions required')
        if not isinstance(data['limitations'], list) or not data['limitations']:
            raise ValueError('Explicit limitations required')
        layers = data['layers']; seen = set()
        if not isinstance(layers, dict) or 'DIE' not in layers:
            raise ValueError('DIE layer is required')
        for name, layer in layers.items():
            _keys(layer, ('gds', 'purpose', 'semantics', 'source'))
            if not name or not isinstance(layer['purpose'],str) or not layer['purpose']:
                raise ValueError('Layer names and purposes must be nonempty strings')
            pair = layer['gds']
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError(f'Explicit layer/datatype pair required: {name}')
            for v in pair:
                integer(v, 'GDS layer/datatype', 0)
                if v > 32767:
                    raise ValueError('GDS layer/datatype out of range')
            if tuple(pair) in seen or layer['source'] not in sources:
                raise ValueError('Duplicate stream mapping or unknown layer source')
            seen.add(tuple(pair))
            if layer['semantics'] not in ('conductor', 'opening', 'junction', 'capacitor', 'inductor', 'boundary'):
                raise ValueError('Unknown layer semantics')
        if layers['DIE']['semantics']!='boundary':
            raise ValueError('DIE must have boundary semantics')
        seen = set()
        if not isinstance(data['rules'], list) or not data['rules']:
            raise ValueError('Rules must be a nonempty list')
        for rule in data['rules']:
            _keys(rule, ('id', 'kind', 'layers', 'value_um', 'source', 'status'), ('ignore_angle_deg',))
            arity = {'width': 1, 'space': 1, 'enclosure': 2, 'inside': 2, 'grid': 1}
            if not isinstance(rule['id'],str) or not rule['id'] or rule['id'] in seen or rule['kind'] not in arity:
                raise ValueError('Duplicate rule ID or unsupported rule kind')
            seen.add(rule['id'])
            if not isinstance(rule['layers'],list) or len(rule['layers']) != arity[rule['kind']] or any(l not in layers for l in rule['layers']):
                raise ValueError('Invalid rule layer references')
            if rule['source'] not in sources or rule['status'] not in ('reference_geometry', 'engineering_assumption'):
                raise ValueError('Unknown rule provenance')
            if 'ignore_angle_deg' in rule and (rule['kind'] not in ('width','space') or
                    not isinstance(rule['ignore_angle_deg'], (int,float)) or not 0 < rule['ignore_angle_deg'] <= 90):
                raise ValueError('ignore_angle_deg is restricted to width/space and (0,90]')
            if rule['kind'] == 'inside':
                if isinstance(rule['value_um'],bool) or rule['value_um'] != 0:
                    raise ValueError('inside uses value_um=0; use enclosure for a margin')
            else:
                positive(rule['value_um'], rule['id'])
            ticks = rule['value_um'] / data['dbu_um']
            if not math.isclose(ticks, round(ticks), abs_tol=1e-8, rel_tol=0):
                raise ValueError('Rule distance must be exactly representable in DBU')
        if not isinstance(data['parameters'], dict):
            raise ValueError('parameters must be a mapping')
        for parameter in data['parameters'].values():
            _keys(parameter, ('value', 'unit', 'source', 'status'))
            positive(parameter['value'], 'parameter')
            if not isinstance(parameter['unit'], str) or not parameter['unit'] or parameter['source'] not in sources:
                raise ValueError('Invalid parameter provenance or unit')
            if parameter['status'] not in ('reference_geometry', 'engineering_assumption'):
                raise ValueError('Invalid parameter status')
        self._data = data

    def to_dict(self):
        return deepcopy(self._data)

    @property
    def dbu(self):
        return self._data['dbu_um']

    @property
    def grid(self):
        return self._data['grid_um']

    @property
    def id(self):
        return self._data['id']

    def layer(self, name):
        return tuple(self._data['layers'][name]['gds'])

    def require_roles(self, roles):
        for name,semantics in roles.items():
            if name not in self._data['layers'] or self._data['layers'][name]['semantics']!=semantics:
                raise ValueError(f'This generator requires {name} with {semantics} semantics')

    def parameter(self, name, unit):
        p = self._data['parameters'][name]
        if p['unit'] != unit:
            raise ValueError(f'{name} requires {unit}, got {p["unit"]}')
        return p['value']

    def coordinate(self, value):
        """Reject off-grid user dimensions; curved polygon vertices are snapped explicitly."""
        if not math.isfinite(value) or not math.isclose(value / self.grid, round(value / self.grid), abs_tol=1e-7, rel_tol=0):
            raise ValueError(f'Off-grid coordinate/dimension: {value} um')
        return round(value / self.dbu)


def load_process(name_or_path='kit_nycreates_research'):
    """Load a bundled profile name or an explicit JSON path."""
    path = Path(name_or_path)
    if path.is_file():
        text = path.read_text(encoding='utf-8-sig')
    else:
        if path.name != str(name_or_path) or path.suffix:
            raise FileNotFoundError(path)
        text = files('twpa.layout').joinpath('processes', f'{name_or_path}.json').read_text(encoding='utf-8')
    return Process(_json(text))
