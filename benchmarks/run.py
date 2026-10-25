"""Benchmarks for coupledlogistic: best-of-k wall time and tracemalloc peak.

    python benchmarks/run.py LABEL [--only SUBSTRING]
"""
import json
import os
import sys
import time
import tracemalloc

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from coupledlogistic import coupledlogistic  # noqa: E402

RESULTS = os.path.join(HERE, 'results.json')
WHEATSTONE = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]]
CHAIN50 = np.eye(50, k=1)


def cases():
    return {
        'kaneko 4 nodes 1e5': lambda: coupledlogistic(10**5, 4, WHEATSTONE, 0.15, 'kaneko', seed=1),
        'diffusive 4 nodes 1e6': lambda: coupledlogistic(10**6, 4, WHEATSTONE, 0.15, 'diffusive', seed=1),
        'kaneko 50 nodes 1e5': lambda: coupledlogistic(10**5, 4, CHAIN50, 0.15, 'kaneko', seed=1),
        '100 realizations 4 nodes 1e4 (loop)': lambda: [coupledlogistic(10**4, 4, WHEATSTONE, 0.15, 'kaneko', seed=s, transient=1000)
                                                        for s in range(100)],
    }


def main():
    label = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else 'latest'
    only = sys.argv[sys.argv.index('--only') + 1] if '--only' in sys.argv else ''
    stored = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    baseline = stored.get('baseline', {})
    results = dict(stored.get(label, {}))
    for name, func in cases().items():
        if only not in name:
            continue
        func()
        times = []
        for __ in range(3):
            start = time.perf_counter()
            func()
            times.append(time.perf_counter() - start)
        tracemalloc.start()
        func()
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        results[name] = r = {'time_ms': round(min(times) * 1000, 2), 'peak_mb': round(peak / 1e6, 2)}
        base = baseline.get(name)
        compare = ''
        if base and label != 'baseline':
            compare = '  (baseline %9.1f ms %7.1f MB, %5.1fx faster)' % (
                base['time_ms'], base['peak_mb'], base['time_ms'] / max(r['time_ms'], 1e-3))
        print('%-38s %9.1f ms %7.1f MB%s' % (name, r['time_ms'], r['peak_mb'], compare), flush=True)
    stored[label] = results
    with open(RESULTS, 'w') as fh:
        json.dump(stored, fh, indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
