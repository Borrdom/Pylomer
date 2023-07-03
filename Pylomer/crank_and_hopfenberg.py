import numpy as np

def crank(t,D,L,w0=0,winfty=1):
    """    
    Model solvent diffusion in a 1D slab according to cranks equation p.48 (equation 4.18)
    Crank, The mathematics of diffusion, 2nd Ed., Clarendon, Oxford, United Kingdom, 1976, 85 (ISBN 0 19 853344 6) 
    Args:
        t (array_like): time /s
        D (float):  Solvent diffusion coefficient  / m^2/s
        L (float):  Thickness /m  
        w0 (float): solvent mass fraction at t=0
        winfty (float): solvent mass fraction in equilibrium
    Returns:
        ndarray:   
        solvent mass fraction at t /
    """
    ninfty=100
    kf=D/L**2
    X0=w0/(1-w0)
    Xinfty=winfty/(1-winfty) if winfty!=1 else 1
    M_M=(1-sum([8/np.pi**2*1/(2*n+1)**2*np.exp(-(2*n+1)**2*kf*t) for n in range(ninfty)]))
    Xt=(1-M_M)*X0+M_M*Xinfty
    return Xt/(1+Xt)



def BH_Model(t,D,L,kr,difffrac,w0=0,winfty=1):
    """
    Model solvent diffusion in a 1D slab with Hophenberg modification for relaxation https://doi.org/10.1016/0032-3861(78)90269-0 
    Args:
        t (array_like): time /s
        D (float):  Solvent diffusion coefficient  / m^2/s
        L (float):  Thickness /m  
        kr (float): time constant for relaxation
        difffrac (float): specifies the fractional solvent amount that is sorbed soley through diffusion 
        w0 (float): solvent mass fraction at t=0
        winfty (float): solvent mass fraction in equilibrium
        
    Returns:
        ndarray:   
        solvent mass fraction at t /
    """
    def Relaxation(t,kr):
        """"used by the Hophenberg modification """
        return (1-np.exp(-kr*t))

    X0=w0/(1-w0)
    Xinfty=winfty/(1-winfty) if winfty!=1 else 1
    M_M=crank(t,D,L)*difffrac+Relaxation(t,kr)*(1-difffrac)
    Xt=(1-M_M)*X0+M_M*Xinfty
    return Xt/(1+Xt)