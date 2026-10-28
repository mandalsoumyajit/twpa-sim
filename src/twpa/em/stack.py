"""Explicit physical stack for the illustrative MIM extraction coupon."""
from dataclasses import asdict, dataclass
from importlib.resources import files
from pathlib import Path
from ..layout.process import _json, positive


@dataclass(frozen=True)
class EMStack:
    id: str
    status: str
    description: str
    dielectric_thickness_um: float
    dielectric_epsilon_r: float
    substrate_thickness_um: float
    substrate_epsilon_r: float
    electrode_model: str = 'PEC'

    def __post_init__(self):
        if self.status not in ('illustrative','local_unqualified'):
            raise ValueError('Only explicitly unqualified EM stacks are supported')
        if not self.id or not self.description or self.electrode_model!='PEC':
            raise ValueError('A named PEC stack with provenance description is required')
        for name in ('dielectric_thickness_um','dielectric_epsilon_r','substrate_thickness_um','substrate_epsilon_r'):
            positive(getattr(self,name),name)
        if self.dielectric_epsilon_r<1 or self.substrate_epsilon_r<1:
            raise ValueError('Relative permittivities must be >=1')

    def to_dict(self):
        return asdict(self)


def load_stack(path=None):
    text=(Path(path).read_text(encoding='utf-8-sig') if path else
          files('twpa.em').joinpath('stacks','mim_illustrative.json').read_text(encoding='utf-8'))
    data=_json(text)
    if data.pop('schema_version',None)!=1:
        raise ValueError('Unsupported EM stack schema')
    try:
        return EMStack(**data)
    except TypeError as e:
        raise ValueError(f'Invalid EM stack fields: {e}') from e
