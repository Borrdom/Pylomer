import numpy as np
from scipy.integrate import solve_ivp
from numba import njit,config

# config.DISABLE_JIT = True

@njit(['f8[:,:](f8, f8[:,::1], b1[::1], f8[:,::1], f8[::1],f8[::1])'],cache=True)
def dcvdt(t,cv,mob,Lvv,ci0,ciB):
    def roll(a): return (a[1:,:]+a[:-1,:])/2.
    def diff(a): return (a[1:,:]-a[:-1,:])
    nz_1,nv=cv.shape
    nc=ci0.shape[0]
    ci=np.ones((nz_1,nc))*ci0
    ci[:,mob]=cv
    ci[-1,:]=ciB
    wi=ci/np.atleast_2d(np.sum(ci,axis=1))
    wi_ = roll(wi)
    ci_ = roll(ci)
    cv_=ci_[:,mob]
    cv=ci[:,mob]
    # dlnwi=np.diff(wi)/wi_
    # dlnci=np.diff(ci)/ci_
    dlncv=diff(cv)/cv_
    # Gij=np.asarray([np.eye(nc)]*nz_1)
    # dlnai=rho*wi_*dlnwi@Gij
    # ji=(rho*wi_*dlnwi)@Lij
    jv=dlncv@Lvv
    djv=diff(np.vstack((np.zeros((1,nv)),jv)))
    dcvdt=np.vstack((djv,np.zeros((1,nv))))
    return dcvdt
    

def diffuse(t,length,Lvec,wi0,wi8,mob,**kwargs):
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
    nv=np.sum(mob)
    Lij=L_Matrix(Lvec/dz**2,nv)
    space=np.linspace(0,length,nz+1)
    
    ct=1200.
    ci0=ct*wi0
    ciB=ct*wi8
    cv0=ci0[mob]*np.ones((nz+1,nv))
    def ode(t,x,ciB):
        cv=np.ascontiguousarray(np.reshape(x,(nz+1,nv)))
        return dcvdt(t,cv,mob,Lij,ci0,ciB).flatten()     
    sol=solve_ivp(ode,(t[0],t[-1]),cv0.flatten(),args=(ciB,),method="BDF",t_eval=t)
    if not sol["success"]: raise Exception(sol["message"])
    x_sol=sol["y"]
    witz=np.zeros((nt,nz+1,nc))
    citz=np.zeros((nt,nz+1,nc))
    for k in range(nt):
        cvtz=np.reshape(x_sol[:(nz+1)*nv,k],(nz+1,nv))
        citz[k,:,:]=ci0
        citz[k,:,mob]=cvtz.T
        citz[k,-1,:]=ciB
        for i in range(nc):
            witz[k,:,i]=citz[k,:,i]/np.sum(citz[k,:,:],axis=1)
    return witz


def L_Matrix(Lvec,nv):
    Lij=np.zeros((nv,nv))
    Lij[np.triu_indices_from(Lij,k=0)]=Lvec
    Lij=Lij.T
    Lij[np.triu_indices_from(Lij,k=0)]=Lvec
    return Lij