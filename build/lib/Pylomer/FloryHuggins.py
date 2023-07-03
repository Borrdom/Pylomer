
import numpy as np

def lnai(wi,Mi,chi=None):
    """
    Compute the glass transition temperature of a mixture  
    
    Args:
        wi (array_like): weight fractions number of components
        Mi (array_like): Molar mass / g/mol
        chi (optional,array_like): Flory Huggins Chi Prameter        /-
    Returns:
        ndarray:   
        logarithmic activity of component i  /K 
    """
    xi=wi/Mi/np.sum(wi/Mi,axis=0)
    ri=Mi/np.min(Mi)
    nc=len(xi)
    chi_mat=np.zeros((nc,nc))
    chi_mat[np.triu_indices(nc, k=1)]=chi
    def Vi(xi): return xi*ri/np.sum(xi*ri,axis=0)
    def delGERT(xi): return np.sum(Vi(xi)/ri*np.log(Vi(xi)),axis=0)+np.sum(chi_mat*np.outer(Vi(xi),Vi(xi)))
    h=1E-26
    idx=-1
    lngi=np.zeros(nc)
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'c16')
        dx[i] = h * 1j
        dx[idx] = - h * 1j 
        lngi[i]=np.imag((delGERT(xi+dx))/h)
    return lngi+np.log(xi)

# Mi=np.asarray([18.015,721])
# rho0i=np.asarray([997.,1150.])
# chi=np.asarray([2.272])
# w1=np.linspace(0.01,0.07,100)
# w2=1-w1

# wi=np.stack((w1,w2))

# lnaivec=np.asarray([lnai(wi[:,i],Mi,chi=chi) for i,val in enumerate(wi[0,:])]).T
# RH=np.exp(lnaivec[0,:])

