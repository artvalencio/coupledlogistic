def coupledlogistic(tslength,r,A,sigma,couplingtype,verbose=False,seed=None,normalize=True,transient=10000,max_retries=10,
                    n_realizations=None,dtype=None,backend='auto'):
    r'''Generates the time-series of a network of coupled logistic maps.

    Each node i follows the logistic map f_i(x) = r_i x (1 - x), coupled to the
    nodes j that link to it (A[j, i] != 0), with in-strength k_i = sum_j A[j, i]:

        couplingtype='diffusive':
            x_{n+1}^i = (1 - sigma) f_i(x_n^i) + sigma/k_i sum_j A[j,i] (x_n^j - x_n^i)
        couplingtype='kaneko':
            x_{n+1}^i = (1 - sigma) f_i(x_n^i) + sigma/k_i sum_j A[j,i] f_j(x_n^j)

    Nodes without incoming links (k_i = 0) follow their own map, x_{n+1}^i = f_i(x_n^i).

    Parameters
    ----------
    tslength: int
        Number of points of each series (after the transient).
    r: float or sequence of floats
        Logistic map parameter: one value for every node, or one per node.
    A: 2-D array-like (nodes x nodes)
        Adjacency matrix: A[j, i] is the weight of the link j -> i (binary or weighted).
    sigma: float
        Coupling strength.
    couplingtype: str
        'diffusive' or 'kaneko' (see above).
    verbose: bool, optional
        Prints progress and retries. Default: False.
    seed: int, numpy Generator or None, optional
        Seed of the random initial conditions, for reproducible series. Default: None.
    normalize: bool, optional
        Rescales each node's series to [0, 1] (constant series become 0). Default: True.
    transient: int, optional
        Number of initial iterations discarded (not stored). Default: 10000.
    max_retries: int, optional
        Orbits that leave [0, 1] (and therefore diverge) or get stuck at 0 are
        restarted from new initial conditions up to this many times; then a
        RuntimeError is raised. Default: 10.
    n_realizations: int or None, optional
        Number of independent orbits (each from its own initial conditions,
        each with its own seed derived from `seed`) to generate at once.
        Default: None (one orbit).
    dtype: numpy float type or None, optional
        Type of the returned series, e.g. numpy.float32 to halve memory (the
        simulation always runs in float64). Default: None (float64).
    backend: str, optional
        'numba' (compiled kernel, ~50x faster; needs the numba package), 'numpy',
        or 'auto' (default: numba when installed). Both follow the same
        equations, but their rounding differs; since the maps are chaotic,
        orbits from the two backends separate after some steps, so a seed
        reproduces a series exactly within the same backend.

    Returns
    -------
    out: numpy array
        Shape (tslength, nodes): column i is the series of node i; with
        n_realizations, shape (n_realizations, tslength, nodes).

    Raises
    ------
    ValueError
        For a non-square A, an unknown couplingtype or backend, or an r vector
        of the wrong length.
    RuntimeError
        When every attempt diverged or got stuck (check r and sigma).

    Examples
    --------
    Simple X -> Y system (node 0 drives node 1):

        A = [[0, 1], [0, 0]]
        x = coupledlogistic(1e5, 4, A, 0.2, 'diffusive', seed=1)

    Serial chain [1]->[2]->[3]->[4], r=4, diffusive coupling 0.2:

        A = [[0,1,0,0], [0,0,1,0], [0,0,0,1], [0,0,0,0]]
        x = coupledlogistic(1e7, 4, A, 0.2, 'diffusive')

    Parallel paths [1]->[2]->[4] and [1]->[3]->[4], coupling 0.1:

        A = [[0,1,1,0], [0,0,0,1], [0,0,0,1], [0,0,0,0]]
        x = coupledlogistic(1e7, 4, A, 0.1, 'diffusive')

    Wheatstone-bridge-like ([1]->[2], [1]->[3], [2]->[3], [2]->[4], [3]->[4]), coupling 0.15:

        A = [[0,1,1,0], [0,0,1,1], [0,0,0,1], [0,0,0,0]]
        x = coupledlogistic(1e7, 4, A, 0.15, 'diffusive')

    100 independent realizations of a Kaneko-coupled pair, as float32:

        x = coupledlogistic(1e4, 4, [[0, 1], [0, 0]], 0.2, 'kaneko', n_realizations=100,
                            dtype=numpy.float32, seed=0)

    (C) Arthur Valencio(1)* and Murilo Baptista(1), 15 December 2017 (Python translation: 16 Sep 2020)
       (1) ICSMB, University of Aberdeen, UK
       *Support: CNPq (Brazil)
    '''

    import numpy as np
    A=np.array(A,dtype=float)
    tslength=int(tslength)
    transient=int(transient)
    if A.ndim!=2 or A.shape[0]!=A.shape[1]:
        raise ValueError('adjacency matrix A must be square')
    if couplingtype not in ('diffusive','kaneko'):
        raise ValueError("couplingtype must be 'diffusive' or 'kaneko'")
    if backend not in ('auto','numba','numpy'):
        raise ValueError("backend must be 'auto', 'numba' or 'numpy'")
    nonodes=A.shape[0]
    r=np.atleast_1d(np.asarray(r,dtype=float))
    if len(r)==1:
        r=np.ones(nonodes)*r[0]
    elif len(r)!=nonodes:
        raise ValueError('r must be a single number or have one value per node')
    dtype=np.dtype(np.float64 if dtype is None else dtype)
    #the compiled kernel is optional: without numba, 'auto' falls back to numpy
    kernel=None
    if backend!='numpy':
        try:
            from ._kernel import simulate_numba as kernel
        except ImportError:
            if backend=='numba':
                raise ImportError("backend='numba' needs numba (pip install numba)") from None

    #Matrix form of the update: x' = a*f(x) + M@v, with v=x (diffusive) or v=f(x) (kaneko).
    #W holds the normalized incoming weights, W[i,j] = A[j,i]/k_i, for coupled nodes;
    #a is (1-sigma) for coupled nodes and 1 for input nodes (pure logistic dynamics).
    #Diffusive coupling sigma/k_i sum_j A[j,i] (x_j - x_i) = sigma (W x - x_i) puts
    #-sigma on the diagonal of M.
    instrength=A.sum(axis=0)
    coupled=instrength!=0
    W=np.zeros_like(A)
    W[coupled,:]=A.T[coupled,:]/instrength[coupled,None]
    a=np.where(coupled,1-sigma,1.0)
    diffusive=couplingtype=='diffusive'
    if diffusive:
        M=sigma*(W-np.diag(coupled.astype(float)))
    else:
        M=sigma*W

    #nonzero entries of M by row (compressed rows, increasing column order): the
    #compiled kernel only visits existing links, so large sparse networks stay fast
    rows,cols=np.nonzero(M)
    rowptr=np.searchsorted(rows,np.arange(nonodes+1))
    vals=M[rows,cols]

    def simulate(x0, out):
        '''Orbits from initial conditions x0 (n,) or (R, n) into out (T, n) or (R, T, n);
        returns which ones stayed in [0,1] and never got stuck at 0.'''
        if kernel is not None:
            if x0.ndim==1:
                return np.array([kernel(x0,r,a,rowptr,cols,vals,transient,diffusive,out)])
            return np.array([kernel(x0[k],r,a,rowptr,cols,vals,transient,diffusive,out[k]) for k in range(len(x0))])
        return _simulate_numpy(x0,r,a,M,transient,diffusive,out,verbose)

    if verbose:
        print('Generating time-series')
    if n_realizations is None:
        #one orbit: draw new initial conditions until an attempt succeeds
        rng=np.random.default_rng(seed)
        out=np.empty((tslength,nonodes),dtype=dtype)
        for attempt in range(max_retries+1):
            if simulate(rng.random(nonodes),out)[0]:
                break
            if verbose:
                print('recalculating...')
        else:
            raise RuntimeError('orbit diverged or degenerated in %d attempts; check r and sigma' % (max_retries+1))
    else:
        #each realization has its own generator, so results do not depend on the batching
        if isinstance(seed,np.random.Generator):
            rngs=[np.random.default_rng(seed.integers(2**63)) for __ in range(n_realizations)]
        else:
            rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(n_realizations)]
        out=np.empty((n_realizations,tslength,nonodes),dtype=dtype)
        #simulate all realizations together; only the failed ones are retried
        pending=np.arange(n_realizations)
        for attempt in range(max_retries+1):
            x0=np.array([rngs[k].random(nonodes) for k in pending])
            block=np.empty((len(pending),tslength,nonodes),dtype=dtype)
            ok=simulate(x0,block)
            out[pending[ok]]=block[ok]
            pending=pending[~ok]
            if len(pending)==0:
                break
            if verbose:
                print('recalculating %d realizations...' % len(pending))
        else:
            raise RuntimeError('orbit diverged or degenerated in %d attempts; check r and sigma' % (max_retries+1))
    if normalize:
        #quick 0 to 1 normalization of each node's series, in place (constant series become 0)
        low=out.min(axis=-2,keepdims=True)
        span=out.max(axis=-2,keepdims=True)-low
        span[span==0]=1
        out-=low
        out/=span
    return out


def _simulate_numpy(x0,r,a,M,transient,diffusive,out,verbose,block=4096):
    '''numpy backend of the simulation.

    Steps are computed in blocks of `block` iterations into a small buffer. After
    each block the new rows are checked (values outside [0,1] mean divergence;
    a node at 0 twice in a row is stuck), so a failing orbit stops within one
    block, and only rows after the transient are copied to `out`.
    x0 is (n,) for one orbit or (R, n) for R orbits computed together.
    Returns a boolean array: which orbits succeeded.
    '''
    import numpy as np
    batch=x0.ndim==2
    R=x0.shape[0] if batch else 1
    n=x0.shape[-1]
    T=out.shape[-2]
    total=transient+T          #number of orbit rows, row 0 being the initial condition
    Mt=M.T
    buf=np.empty((block+1,)+x0.shape)  #buf[0] is the last row of the previous block
    buf[0]=x0
    ok=np.ones(R,dtype=bool)
    if transient==0:
        out[...,0,:]=x0
    report=max(1,total//10)
    step=0                     #orbit index of buf[0]
    with np.errstate(all='ignore'):  #diverging orbits overflow; they are rejected below
        while step<total-1:
            m=min(block,total-1-step)
            x=buf[0]
            for k in range(1,m+1):
                f=r*x*(1-x)
                v=x if diffusive else f
                x=buf[k]=a*f+(v@Mt if batch else M@v)
            #check the block (including the step from the previous block's last row)
            segment=buf[:m+1].reshape(m+1,R,n)
            new,old=segment[1:],segment[:-1]
            failed=~((new>=0)&(new<=1)).all(axis=(0,2))|((new==0)&(old==0)).any(axis=(0,2))
            ok&=~failed
            if not ok.any():
                return ok
            #store the rows of this block that come after the transient
            first,last=max(step+1,transient),step+m
            if first<=last:
                rows=segment[first-step:last-step+1]
                out_rows=out.reshape(R,T,n)[:,first-transient:last-transient+1]
                out_rows[...]=rows.transpose(1,0,2)
            if verbose and (step//report)!=((step+m)//report):
                print('.',end='')
            buf[0]=buf[m]
            step+=m
    return ok
