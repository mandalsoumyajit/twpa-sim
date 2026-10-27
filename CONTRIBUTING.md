# Contributing

Install a source checkout with `python -m pip install -e ".[plot]"`.
The core models need NumPy and SciPy; JoSIM remains optional and external.

Run the relevant checks before proposing a change:

```powershell
python -m unittest discover -s tests -v
python tools/check_documentation.py --run-examples
```

For a small documentation correction, the local-link check alone is usually
sufficient. Run changed examples when their inputs or code change. The
workflow also checks the full numerical suite on Python 3.10 and 3.12.

## Numerical and model changes

Explain the physical problem, units, amplitude normalization and supported
domain. Preserve existing public conventions or document a migration.
Use analytic limits and independent comparisons where possible; generating
a fixture from the same calculation is not independent validation.

For quantum noise, include loss baths as well as damping. Check commutators,
thermal equilibrium, covariance uncertainty and the appropriate finite-gain
or quadrature limit. Do not clip or renormalize a result to hide a defect.
Distinguish Gaussian/classical approximations from broader quantum dynamics.

For circuit validation, record physical parameters, ports, source phases,
initial/bias history, simulator build, controls and measurement windows.
Pin the executable for a convergence study. Retain successful netlists,
CSV hashes and logs locally under ignored `validation/runs/`; commit compact
measurements with provenance. Do not replace measured simulator outputs
with predictions labeled as measurements. Check timestep/settling and
spectral truncation separately when possible.

## Documentation and references

Keep the [HOWTO](HOWTO.md), [documentation map](docs/README.md) and device
guide aligned with an API change. New example scripts should be commented,
import safely via a `main()` guard, and run with small synthetic parameters.
External simulator execution should be explicitly requested.

Add relevant primary literature to [references.bib](references.bib) and
explain its role in [References](docs/references.md). State whether a paper
is background, a model reference or an actually reproduced benchmark.
Contributions are covered by the repository's [MIT License](LICENSE).
