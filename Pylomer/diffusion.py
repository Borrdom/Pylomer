import numpy as np
from scipy.integrate import solve_ivp
from numba import njit

@njit(['f8[:,:](f8, f8[:,::1], b1[::1], f8[:,:], f8[::1],f8[:,:],f8[::1]'],cache=True)
def dcvdt(t,cv,mob,Lvv,ci0,ciB):
    def roll(a): return (a[1:,:]+a[:-1,:])/2.
    nv,nz_1=ck.shape
    nc=ci0.shape
    ci=np.ones((nz_1,nc))*ci0
    ci[:,mob]=cv
    ci[-1,:]=ciB
    wi=ci/np.sum(ci,axis=1)
    wi_ = roll(wi)
    ci_ = roll(ci)
    cv_=ci_[:,mob]
    cv=ci[:,mob]
    # dlnwi=np.diff(wi)/wi_
    dlnci=np.diff(ci)/ci_
    dlncv=np.diff(cv)/cv_
    # Gij=np.asarray([np.eye(nc)]*nz_1)
    # dlnai=rho*wi_*dlnwi@Gij
    # ji=(rho*wi_*dlnwi)@Lij
    jv=dlncv@Lvv
    djv=np.diff(np.hstack((np.zeros((1,nv)),jv)))
    dcvdt=np.hstack((dji,np.zeros((1,nv))))
    return dcvdt
    

def diffusion(t,length,Lvec,wi0,wi8,mob,**kwargs):
    """
    Args:
        t (array_like): time
        length (float) : diffusion distance /m
        Lij (array_like): matrix of linear coefficients /m^2/s
        wi0 (array_like): Mass fractions at t=0               /-
        wi8 (array_like): Mass fraction at t=infinity         /-
        mob(array_like): boolean vector indicating the mobility of a component
    Returns:
        Matrix of mass fractions at t,z     /- \n
    """
    nc=len(wi0)
    nt=len(t)
    nz=20
    dz=length/nz
    Lij=L_Matrix(Lvec/dz**2,nc)
    space=np.linspace(0,length,nz+1)
    nv=np.sum(mob)
    ct=1200.
    wi0_nomob=wi0[~mob]/np.sum(wi0[~mob])
    ci0=ct*wi0*np.ones((nz+1,nc))
    cv0=ci0[:,mob]
    ciB=wi8*rho
    def ode(t,x,ciB):
        cv0=np.ascontiguousarray(np.reshape(x,(nz+1,nk)))
        return dckdt(t,cv,mob,Lij,ci0,ciB)     
    sol=solve_ivp(ode,(t[0],t[-1]),cv0,args=(ciB),method="BDF",t_eval=t)
    if not sol["success"]: raise Exception(sol["message"])
    x_sol=sol["y"]
    
    wit=np.zeros((nt,nz+1,nc))
    cit=np.zeros((nt,nz+1,nc))
    for k in range(nt):
        cvt=np.reshape(x_sol[:(nz+1)*nv,k],(nz+1,nv))
        cit[k,:,:]=ci0
        cit[k,:,:]=cv0
        cit[k,-1,:]=ciB
        for i in range(nc):
            wik[k,:,i]=cit[k,:,i]/np.sum(cit[k,:,:],axis=1)
    return wik


def L_Matrix(Lvec,nv):
    Lij=np.zeros((nv,nv))
    Lij[np.triu_indices_from(D,k=0)]=Lvec
    Lij=Lij.T
    Lij[np.triu_indices_from(D,k=0)]=Lvec
    return Lij