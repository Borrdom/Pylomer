
import numpy as np

def lnai(wi,rho0i,Mi,chi=None):
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
    v0i=Mi/rho0i
    ri=v0i/v0i[np.argmin(Mi)]
    # ri=Mi/np.min(Mi)
    nc=len(xi)
    chi_mat=np.zeros((nc,nc))
    chi_mat[np.triu_indices(nc, k=1)]=chi
    def Vi(xi): return xi*ri/np.sum(xi*ri,axis=0)
    def delGERT(xi): return np.sum(xi*np.log(Vi(xi)/xi),axis=0)+np.sum(chi_mat*np.outer(Vi(xi),Vi(xi))*np.sum(xi*ri,axis=0))
    h=1E-26
    idx=-1
    lngi=np.zeros(nc)
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'c16')
        dx[i] = h * 1j
        dx[idx] = - h * 1j 
        lngi[i]=np.imag((delGERT(xi+dx))/h)
    return lngi+np.log(xi)

# Mi=np.asarray([18.015,357.787])
# rho0i=np.asarray([997.,1320.])
# chi=np.asarray([2.7067])
# w1=np.linspace(0.001,0.05,100)
# w2=1-w1

# wi=np.stack((w1,w2))

# lnaivec=np.asarray([lnai(wi[:,i],rho0i,Mi,chi=chi) for i,val in enumerate(wi[0,:])]).T
# RH=np.exp(lnaivec[0,:])

# import matplotlib.pyplot as plt
# plt.plot(RH,w1)
# plt.show()
# plt.show()
