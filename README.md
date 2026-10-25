# Coupled Logistic

Generates a time-series from a coupled logistic map system.
This dynamical system is useful as a workbench for nonlinear
and causality analytical tools.

-----------------------------------------------------------------------------
-----------------------------------------------------------------------------

Requirements:

Python>=3.8 and numpy (installed automatically).

-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
Installation:

    pip install git+https://github.com/artvalencio/coupledlogistic

Now you can import the package in Python3 as: "import coupledlogistic"

-----------------------------------------------------------------------------
Usage:

    from coupledlogistic import coupledlogistic

    # node 0 drives node 1 (A[j, i] is the link j -> i)
    A = [[0, 1],
         [0, 0]]
    x = coupledlogistic(100000, 4, A, 0.2, 'diffusive', seed=42)
    x.shape          # (100000, 2): one column per node, rescaled to [0, 1]

Parameters: tslength, r (one value or one per node), A (binary or weighted
adjacency matrix), sigma (coupling strength), couplingtype ('diffusive' or
'kaneko'), and optionally seed, normalize (default True), transient
(default 10000), max_retries (default 10) and verbose.

Orbits that leave [0, 1] (and therefore diverge) are restarted from new
initial conditions; if every attempt fails a RuntimeError is raised.

Changes in 0.1.0: output has exactly tslength rows; reproducible output
with seed; Kaneko coupling uses r of the source node; weighted links;
diverging orbits are detected instead of returned as NaN; about 5x faster.

-----------------------------------------------------------------------------
Running the tests:

    pip install -e .[test]
    pytest

-----------------------------------------------------------------------------

-----------------------------------------------------------------------------

Licence and acknoledgements

Open source licence: MIT

Developers: Arthur Valencio [1], Murilo S. Baptista [2]

Institutions: 

[1] Institute of Computing/University of Campinas and Research, Innovation and Dissemination Center for Neuromathematics

[2] Institute for Complex Systems and Mathematical Biology/University of Aberdeen

Funders: FAPESP (S. Paulo Research Foundation), Brazil (2018-onwards). CNPq (National Council for Scientific and Technologic Development), Brazil (2015-2018).

Acknowledgements: Arthur Valencio is supported by the FAPESP grant 2018/09900-8 and 	     acknowledges this work as part of the activities of the FAPESP Research 	     Innovation and Dissemination Center for Neuromathematics (grant 2013/07699-0, S. Paulo Research Foundation). A.V. contributions are partly derived from previous research funded by CNPq grant 206246/2014-5.

Licence text:

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

-----------------------------------------------------------------------------