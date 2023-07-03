import numpy as np

def Tggt(wi,Tgi,q,Ki=None,rhoi=None):
    """
    Compute the glass transition temperature of a mixture  
    
    Args:
        wi (array_like): 2D Array of weight fractions [number of Points, number of components]
        Tgi (array_like): pure component glass transition temperature /K
        q (array_like): Kwei parameter /-
        rhoi (optional,array_like) : pure component densities /kg/m^3
        Ki (optional,array_like): Gordon-Taylor parameters         /-
    Returns:
        ndarray:   
        glass transition temperature of a mixture  /K 
    """
    nc=Tgi.shape[0]
    qmat=np.zeros((nc,nc))
    qmat[np.triu_indices(nc, k=1)]=q
    #qmat=qmat[0]
    Excess=np.asarray([np.sum(np.prod(np.meshgrid(wi[i,:],wi[i,:]),axis=0)*qmat) for i,val in enumerate(wi[:,0])])
    Ideal=np.sum(wi*1/rhoi,1)/np.sum(wi*1/rhoi/Tgi,axis=1)
    #Ideal=np.sum(wi*Tgi*ki,1)/np.sum(wi*ki,1)
    return Ideal+Excess
    #Triangular numbers for binary interaction parameters
