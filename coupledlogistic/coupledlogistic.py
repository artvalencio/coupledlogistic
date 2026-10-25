def coupledlogistic(tslength,r,A,sigma,couplingtype,verbose=False,seed=None,normalize=True,transient=10000,max_retries=10,
                    n_realizations=None,dtype=None,backend='auto'):
    '''Generates time-series dynamics for coupled logistic networks of parameter r
    --------------------------------
    Inputs:
           tslength: length of the time-series (number of points)
           r: logistic map free parameter. Can be a single number (all nodes
              are equal) or a vector (specifying different r for each node). 
           A: adjacency matrix, A[j,i] is the weight of the link j->i
              (binary or weighted)
           sigma: coupling strength
           couplingtype: one of the options: 'diffusive' or 'kaneko'
           verbose: whether to display progress bar or not
           seed: seed (int) or numpy Generator for the random initial
                 conditions, for reproducible series. Default: None
           normalize: whether to rescale each node's series to [0,1].
                 Default: True
           transient: number of initial iterations discarded. Default: 10000
           max_retries: orbits that leave [0,1] (and so diverge) or get stuck
                 at 0 are restarted from new initial conditions up to this
                 many times before a RuntimeError is raised. Default: 10
           n_realizations: number of independent orbits (each from its own
                 initial conditions) to generate at once. Default: None (one orbit)
           dtype: type of the returned series, e.g. numpy.float32 to halve
                 memory (the simulation itself runs in float64). Default: float64
           backend: 'numba' (compiled, much faster; needs numba), 'numpy', or
                 'auto' (default: numba when installed). Both follow the same
                 equations, but for chaotic dynamics their rounding differs
                 and orbits separate after some steps: a seed reproduces a
                 series exactly within the same backend.
    --------------------------------
    Output:
           out: each column is time-series for a node 
                described in the adjacency matrix
                (shape: tslength x nodes, or n_realizations x tslength x nodes)
    --------------------------------
    Usage examples:
       - Simple X->Y system:
           A=[[0,1],[0,0]]
           x=coupledlogistic(1e5,4,A,0.2,'diffusive')
    
    
     - Serial:
           A=[[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]]
           x=coupledlogistic(1e7,4,A,0.2,'diffusive')
            
           (This adjacency matrix A defines a system with nodes [i] connected
            as: 
                               [1]->[2]->[3]->[4]. 
            The dynamics of each node is a r=4 logistic map.
            Each link is a linear diffusive coupling with strength 0.2.
            The time-series from each node has 1*10^7 points)
    
     - Parallel:
           A=[[0,1,1,0],[0,0,0,1],[0,0,0,1],[0,0,0,0]]
           x=coupledlogistic(1e7,4,A,0.1,'diffusive')
            
           (This adjacency matrix A defines a system with nodes [i] connected
            as: 
                                   [1]->[2]
                                    |    |
                                    V    V
                                   [3]->[4]. 
            The dynamics of each node is a r=4 logistic map.
            Each link is a linear diffusive coupling with strength 0.1.
            The time-series from each node has 1*10^7 points)
    
     - Wheatstone-bridge-like:
           A=[[0,1,1,0],[0,0,1,1],[0,0,0,1],[0,0,0,0]]
           x=coupledlogistic(1e7,4,A,0.15,'diffusive')
            
           (This adjacency matrix A defines a system with nodes [i] connected
            as: 
                                   [1]->[2]
                                    |  / |
                                    V V  V
                                   [3]->[4]. 
            The dynamics of each node is a r=4 logistic map.
            Each link is a linear diffusive coupling with strength 0.15.
            The time-series from each node has 1*10^7 points)
    --------------------------------
    LaTeX expression:
    
           If couplingtype='diffusive':
           $x_{n+1}^i=(1-\sigma)f(x_n^i)+\frac{\sigma}{k_i}\sum_j{A_{ji}(x_n^j-x_n^i)}$
           Or if couplingtype='kaneko':
           $x_{n+1}^i=(1-\sigma)f(x_n^i)+\frac{\sigma}{k_i}\sum_j{A_{ji}f(x_n^j)}$
           where $f(x)=r*x*(1-x)$
           
    --------------------------------
    (C) Arthur Valencio(1)* and Murilo Baptista(1), 15 December 2017 (Python translation: 16 Sep 2020)
       (1)ICSMB, University of Aberdeen,UK
      *Support: CNPq (Brazil)'''

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
    kernel=None
    if backend!='numpy':
        try:
            from ._kernel import simulate_numba as kernel
        except ImportError:
            if backend=='numba':
                raise ImportError("backend='numba' needs numba (pip install numba)") from None

    #matrix form of the update: x' = a*f(x) + M@v, with v=x (diffusive) or v=f(x) (kaneko)
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

    #nonzero entries of M by row (increasing column order), for the compiled kernel
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
    '''numpy backend of the simulation: steps are computed in blocks through a small
    buffer, so a failing orbit stops within one block and the transient is not stored.'''
    import numpy as np
    batch=x0.ndim==2
    R=x0.shape[0] if batch else 1
    n=x0.shape[-1]
    T=out.shape[-2]
    total=transient+T
    Mt=M.T
    buf=np.empty((block+1,)+x0.shape)
    buf[0]=x0
    ok=np.ones(R,dtype=bool)
    if transient==0:
        out[...,0,:]=x0
    report=max(1,total//10)
    step=0
    with np.errstate(all='ignore'):
        while step<total-1:
            m=min(block,total-1-step)
            x=buf[0]
            for k in range(1,m+1):
                f=r*x*(1-x)
                v=x if diffusive else f
                x=buf[k]=a*f+(v@Mt if batch else M@v)
            segment=buf[:m+1].reshape(m+1,R,n)
            new,old=segment[1:],segment[:-1]
            failed=~((new>=0)&(new<=1)).all(axis=(0,2))|((new==0)&(old==0)).any(axis=(0,2))
            ok&=~failed
            if not ok.any():
                return ok
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
