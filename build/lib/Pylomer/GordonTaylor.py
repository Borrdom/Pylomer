import numpy as np

def Tggt(wi,Tg0i,q,Ki=None,rho0i=None):
    """
    Compute the glass transition temperature of a mixture  
    
    Args:
        wi (array_like): 2D Array of weight fractions [ number of components,number of Points]
        Tg0i (array_like): pure component glass transition temperature /K
        q (array_like): Kwei parameter /-
        rho0i (optional,array_like) : pure component densities /kg/m^3
        Ki (optional,array_like): Gordon-Taylor parameters         /-
    Returns:
        ndarray:   
        glass transition temperature of a mixture  /K 
    """
    nc=Tg0i.shape[0]
    qmat=np.zeros((nc,nc))
    qmat[np.triu_indices(nc, k=1)]=q
    qmat[np.tril_indices(nc, k=-1)]=q
    Excess=np.asarray([np.sum(np.prod(np.meshgrid(wi[:,i],wi[:,i]),axis=0)*qmat) for i,val in enumerate(wi[0,:])])
    Ideal=np.sum(wi*1/rho0i[:,None],axis=0)/np.sum(wi*1/rho0i[:,None]/Tg0i[:,None],axis=0)
    return Ideal+Excess
