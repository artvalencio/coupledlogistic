'''Compiled orbit kernel (optional: needs numba).

Same update as the numpy backend, written as plain loops that numba compiles
to machine code (cache=True stores the compiled code on disk, so only the very
first call pays the compilation time).
'''
import numba
import numpy as np


@numba.njit(cache=True)
def simulate_numba(x0, r, a, rowptr, cols, vals, transient, diffusive, out):
    '''Iterates x' = a*f(x) + M@v (v = x if diffusive else f(x)) from x0, with M in
    compressed-row form (row i holds vals[rowptr[i]:rowptr[i+1]] at columns cols[...]).

    Rows transient, transient+1, ... of the orbit (row 0 is x0) are written to out.
    Returns False as soon as a value leaves [0, 1] or a node is 0 twice in a row.
    '''
    n = x0.shape[0]
    total = transient + out.shape[0]
    x = x0.copy()
    f = np.empty(n)
    new = np.empty(n)
    if transient == 0:
        for i in range(n):
            out[0, i] = x[i]
    for t in range(1, total):
        # logistic map of every node, then the coupling of each node with its sources
        for i in range(n):
            f[i] = r[i] * x[i] * (1.0 - x[i])
        for i in range(n):
            acc = 0.0
            for k in range(rowptr[i], rowptr[i + 1]):
                j = cols[k]
                acc += vals[k] * (x[j] if diffusive else f[j])
            v = a[i] * f[i] + acc
            # leaving [0, 1] means the orbit diverges; stop at once (NaN fails too)
            if not (v >= 0.0 and v <= 1.0):
                return False
            if v == 0.0 and x[i] == 0.0:
                return False
            new[i] = v
        for i in range(n):
            x[i] = new[i]  # synchronous update: all nodes use the previous step's values
        if t >= transient:  # the transient is iterated but not stored
            for i in range(n):
                out[t - transient, i] = x[i]
    return True
