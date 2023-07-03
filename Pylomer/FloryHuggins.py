
import numpy as np

def lngi(wi,Mi,rho0i,chi=None):
    """
    Compute the glass transition temperature of a mixture  
    
    Args:
        wi (array_like): 2D Array of weight fractions [number of components,number of Points]
        Mi (array_like): pure component glass transition temperature /K
        rho0i (optional,array_like) : pure component densities /kg/m^3
        Ki (optional,array_like): Gordon-Taylor parameters         /-
    Returns:
        ndarray:   
        glass transition temperature of a mixture  /K 
    """
    xi=wi/Mi/np.sum(wi/Mi,axis=0)
    Voi=Mi/rho0i/1000.
    nc=len(xi)
    chi_mat=np.zeros((nc,nc))
    chi_mat[np.triu_indices(nc, k=1)]=chi
    chi_mat[np.tril_indices(nc,k=-1)]=chi
    def Vi(xi): return xi*Voi/np.sum(xi*Voi,axis=0)
    def delGERT(xi): return np.sum(xi*np.log(Vi(xi)),axis=0)+np.sum((chi_mat@Vi(xi))*xi,axis=0)
    h=1E-26
    idx=-1
    lnai=np.zeros(nc)
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'c16')
        dx[i] = h * 1j
        dx[idx] = - h * 1j 
        lnai[i]=np.imag((delGERT(xi+dx))/h)
    return lnai-np.log(xi)

