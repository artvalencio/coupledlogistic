# Coupled Logistic

Generates time-series from networks of coupled logistic maps. These chaotic
systems with a known coupling structure are a workbench for nonlinear and
causality analysis tools such as [CaMI](https://github.com/artvalencio/cami-python).

## Installation

```
pip install git+https://github.com/artvalencio/coupledlogistic
pip install "coupledlogistic[fast] @ git+https://github.com/artvalencio/coupledlogistic"   # with the compiled backend
```

Requirements: Python 3.9 or later and numpy (>= 1.20). Optional: numba
(>= 0.55) for the compiled backend, about 50x faster. See `requirements.txt`.

## Usage

```python
import numpy as np
from coupledlogistic import coupledlogistic

# node 0 drives node 1 (A[j, i] is the link j -> i)
A = [[0, 1],
     [0, 0]]
x = coupledlogistic(100000, 4, A, 0.2, 'diffusive', seed=42)
x.shape                    # (100000, 2): one column per node, rescaled to [0, 1]

# heterogeneous nodes, weighted links, Kaneko coupling, raw values
x = coupledlogistic(10000, [3.9, 3.8, 4.0], [[0, 0.5, 1], [0, 0, 2], [0, 0, 0]], 0.1,
                    'kaneko', normalize=False)

# 100 independent realizations at once, stored as float32
x = coupledlogistic(10000, 4, A, 0.2, 'kaneko', n_realizations=100, dtype=np.float32, seed=0)
x.shape                    # (100, 10000, 2)
```

### Model

Each node follows `f_i(x) = r_i x (1 - x)` and is coupled to the nodes that
link to it, with in-strength `k_i = sum_j A[j, i]`:

- diffusive: `x_i' = (1 - sigma) f_i(x_i) + sigma/k_i sum_j A[j,i] (x_j - x_i)`
- Kaneko: `x_i' = (1 - sigma) f_i(x_i) + sigma/k_i sum_j A[j,i] f_j(x_j)`

Nodes without incoming links follow their own map.

### Parameters

| Parameter | Meaning |
|---|---|
| `tslength` | points per series (after the transient) |
| `r` | logistic parameter, one value or one per node |
| `A` | adjacency matrix, `A[j, i]` = weight of the link j → i (binary or weighted) |
| `sigma` | coupling strength |
| `couplingtype` | `'diffusive'` or `'kaneko'` |
| `seed` | int or numpy Generator, for reproducible series |
| `normalize` | rescale each series to [0, 1] (default `True`) |
| `transient` | iterations discarded first (default 10000) |
| `max_retries` | diverging or stuck orbits are restarted this many times, then `RuntimeError` (default 10) |
| `n_realizations` | number of independent orbits returned together (default: one) |
| `dtype` | type of the output, e.g. `numpy.float32` |
| `backend` | `'auto'` (numba when installed), `'numba'` or `'numpy'` |
| `verbose` | print progress |

Orbits that leave [0, 1] diverge; they are detected as soon as they do and
restarted from new initial conditions.

**Backends and reproducibility.** Both backends follow the equations exactly,
but their floating-point rounding differs; because the maps are chaotic,
their orbits separate after some steps. A seed reproduces a series exactly
within the same backend (pass `backend=` to pin one).

## Performance

| Case | 0.01 | 0.1.0 |
|---|---|---|
| 4 nodes, 10^6 steps | 5.8 s / 96 MB | 0.09 s / 32 MB (numba); 16 MB with `dtype=float32` |
| 50 nodes, 10^5 steps | 0.72 s | 0.05 s (numba) |
| 100 realizations × 10^4 steps | 6.2 s (loop) | 0.28 s (numpy backend, `n_realizations=100`) |

Benchmarks: `python benchmarks/run.py LABEL`.

## Changes in 0.1.0

Output has exactly `tslength` rows; seeds make series reproducible; Kaneko
coupling uses the source node's `r`; weighted links; diverging orbits are
detected instead of being returned as NaN; single-node networks and input
validation; compiled backend, batched realizations and float32 output.

## Tests

```
pip install -e .[test]
pytest              # fast tests (both backends when numba is installed)
pytest -m slow      # time and memory budgets
```

Every orbit is checked step by step against the model equations.

## Licence and acknowledgements

Open source licence: MIT

Developers: Arthur Valencio [1], Murilo S. Baptista [2]

[1] Institute of Computing/University of Campinas and Research, Innovation and Dissemination Center for Neuromathematics
[2] Institute for Complex Systems and Mathematical Biology/University of Aberdeen

Funders: FAPESP (S. Paulo Research Foundation), Brazil (2018-onwards). CNPq (National Council for Scientific and Technologic Development), Brazil (2015-2018).

Acknowledgements: Arthur Valencio is supported by the FAPESP grant 2018/09900-8 and acknowledges this work as part of the activities of the FAPESP Research Innovation and Dissemination Center for Neuromathematics (grant 2013/07699-0, S. Paulo Research Foundation). A.V. contributions are partly derived from previous research funded by CNPq grant 206246/2014-5.

```
Copyright (c) 2018 The Python Packaging Authority

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
