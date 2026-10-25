def coupledlogistic(tslength,r,A,sigma,couplingtype,verbose=False,seed=None,normalize=True,transient=10000,max_retries=10):
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
    --------------------------------
    Output:
           out: each column is time-series for a node 
                described in the adjacency matrix
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
    nonodes=A.shape[0]
    r=np.atleast_1d(np.asarray(r,dtype=float))
    if len(r)==1:
        r=np.ones(nonodes)*r[0]
    elif len(r)!=nonodes:
        raise ValueError('r must be a single number or have one value per node')
    rng=np.random.default_rng(seed)
    total=transient+tslength

    def step(x):
        new=np.empty(nonodes)
        for k in range(nonodes):
            sumterm=0
            deg=0
            for l in range(nonodes):
                if A[l,k]!=0:
                    if couplingtype=='diffusive':
                        sumterm=sumterm+A[l,k]*(x[l]-x[k])
                    else:
                        sumterm=sumterm+A[l,k]*r[l]*x[l]*(1-x[l])
                    deg=deg+A[l,k]
            if deg>0:
                new[k]=(1-sigma)*r[k]*x[k]*(1-x[k])+sigma*sumterm/deg
            else: #deg=0 means it's an input node, so calc only logistic dynamics
                new[k]=r[k]*x[k]*(1-x[k])
        return new

    def simulate():
        out=np.empty([total,nonodes])
        out[0,:]=rng.random(nonodes)
        for n in range(total-1):
            if verbose and n % max(1,total//10)==0:
                print('.',end='')
            out[n+1,:]=step(out[n,:])
            #diverging or degenerate orbit: repeat for new initial cond
            escaped=not np.all((out[n+1,:]>=0) & (out[n+1,:]<=1))
            if escaped or np.any((out[n+1,:]==0) & (out[n,:]==0)):
                return None
        return out

    if verbose:
        print('Generating time-series')
    for attempt in range(max_retries+1):
        out=simulate()
        if out is not None:
            break
        if verbose:
            print('recalculating...')
    else:
        raise RuntimeError('orbit diverged or degenerated in %d attempts; check r and sigma' % (max_retries+1))
    out=out[transient:,:]
    if normalize:
        #quick 0 to 1 normalization (constant series become 0)
        span=out.max(axis=0)-out.min(axis=0)
        span[span==0]=1
        out=(out-out.min(axis=0))/span
    return out
