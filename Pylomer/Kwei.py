import numpy as np

def Tggt(wi,rhoi,Tgi,q):
    nc=rhoi.shape[0]
    #ki=rhoi[idx]*Tgi[idx]/(rhoi*Tgi)
    qmat=np.zeros((nc,nc))
    qmat[np.triu_indices(nc, k=1)]=q
    #qmat=qmat[0]
    Excess=np.asarray([np.sum(np.prod(np.meshgrid(wi[i,:],wi[i,:]),axis=0)*qmat) for i,val in enumerate(wi[:,0])])
    Ideal=np.sum(wi*1/rhoi,1)/np.sum(wi*1/rhoi/Tgi,1)
    #Ideal=np.sum(wi*Tgi*ki,1)/np.sum(wi*ki,1)
    return Ideal+Excess
    #Triangular numbers for binary interaction parameters
