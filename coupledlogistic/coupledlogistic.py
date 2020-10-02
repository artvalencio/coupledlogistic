def coupledlogistic(tslength,r,A,sigma,couplingtype,verbose=False):
    '''Generates time-series dynamics for coupled logistic networks of parameter r
    --------------------------------
    Inputs:
           tslength: length of the time-series (number of points)
           r: logistic map free parameter. Can be a single number (all nodes
              are equal) or a vector (specifying different r for each node). 
           A: adjacency matrix
           sigma: coupling strength
           couplingtype: one of the options: 'diffusive' or 'kaneko'
           verbose: whether to display progress bar or not
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

    import numpy as  np
    import pandas as pd
    A=np.array(A)
    tslength=int(tslength)
    if verbose==True:
        print('Generating time-series')
    nonodes=len(A[1,:])
    
    if type(r)==int or type(r)==float:
        r=np.ones(nonodes)*r
    elif len(r)==1:
        r=np.ones(nonodes)*r[0]


    def diffusivecalc(tslength,r,A,sigma,nonodes):
    #calculation when diffusive
        cond=1
        count=0
        while cond:    
            temp=0
            count=count+1
            #initial cond
            out=np.full([tslength+10001,nonodes],np.nan)
            out[0,:]=np.random.random(nonodes)
            #calculate
            for n in range(tslength+10000):
                #progress bar
                if verbose==True:
                    if n % np.floor(tslength/10)==0:
                        print('.',end='')
                #actual calc
                deg=[]
                for k in range(nonodes):
                        #calc coupling    
                        sumterm=0
                        deg.append(0)
                        for l in range(nonodes):
                            if A[l,k]==1:
                                sumterm=sumterm+out[n,l]-out[n,k]
                                deg[k]=deg[k]+1
                        #calc next step
                        if deg[k]>0:
                            sumterm=sumterm/deg[k]
                            out[n+1,k]=(1-sigma)*r[k]*out[n,k]*(1-out[n,k])+sigma*sumterm
                        else: #deg=0 means it's an input node, so calc only logistic dynamics
                            out[n+1,k]=r[k]*out[n,k]*(1-out[n,k])
                        
                        #error: repeat for new initial cond
                        if out[n+1,k]==np.nan:
                            temp=1   
                            break
                        elif out[n+1,k]==0 and out[n,k]==0:
                            temp=1
                            break
            #error handling
            if temp==0:
                break
            else:
                print('recalculating...', end=' ')
            print(count)
            if count>10:
                break
        return out

    def kanekocalc(tslength,r,A,sigma,nonodes):
    #calculation when kaneko
        cond=1
        count=0
        while cond:
            count=count+1
            temp=0
            #initial cond
            out=np.full([tslength+10001,nonodes],np.nan)
            out[0,:]=np.random.random(nonodes)
            #calculate
            for n in range(tslength+10000):
                #progress bar
                if verbose==True:
                    if n % np.floor(tslength/10)==0:
                        print('.',end='')
                #actual calc
                deg=[]
                for k in range(nonodes):
                        #calc coupling    
                        sumterm=0
                        deg.append(0)
                        for l in range(nonodes):
                            if A[l,k]==1:
                                sumterm=sumterm+r[k]*out[n,l]*(1-out[n,l])
                                deg[k]=deg[k]+1
                        #calc next step
                        if deg[k]>0:
                            sumterm=sumterm/deg[k]
                            out[n+1,k]=(1-sigma)*r[k]*out[n,k]*(1-out[n,k])+sigma*sumterm
                        else: #deg=0 means it's an input node, so calc only logistic dynamics
                            out[n+1,k]=r[k]*out[n,k]*(1-out[n,k])
                        
                        #error: repeat for new initial cond
                        if out[n+1,k]==np.nan:
                            temp=1   
                            break
                        elif out[n+1,k]==0 and out[n,k]==0:
                            temp=1
                            break           
            #error handling
            if temp==0:
                break
            else:
                print('recalculating...',end=' ')
            print(count)
            if count>10:
                break
        return out

    def normal(x):
        #NORMAL Quick 0 to 1 nomalization
        return (x-min(x))/(max(x)-min(x))
        
    if couplingtype=='diffusive':
        out=diffusivecalc(tslength,r,A,sigma,nonodes)
    elif couplingtype=='kaneko':
        out=kanekocalc(tslength,r,A,sigma,nonodes)
    
    #cut transient
    out=out[10002:,:]
    #normalize
    for i in range(nonodes):
        out[:,i]=normal(out[:,i])
    return out
