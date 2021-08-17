import casadi as cs
import numpy as np
R,Nav=8.31448,6.022E23
class FloryFilm:
    T=298.15
    p=1.013E5
    def __init__(self,comp,chi):
        self.comp=comp
        self.nc=len(comp)
        self.chi=chi

    def FH(self):
        nc=self.nc
        comp=self.comp
        ni=cs.SX.sym("ni",nc)
        Mi=cs.SX.ones(nc)
        rho0i=cs.SX.ones(nc)
        for i in range(nc):
            Mi[i]=comp[i]["Mi"]
            rho0i[i]=comp[i]["rho0"]
        Voi=Mi/rho0i
        wi=ni*Mi/cs.sum1(ni*Mi)
        xhi=np.zeros((nc,nc))
        xhi[np.triu_indices(nc, k=1)]=self.chi
        xhi=cs.SX(xhi)
        #phii=xi*Mi/rho0i/cs.sum1(xi*Mi/rho0i)
        Vi=ni*Voi/cs.sum1(ni*Voi)
        #wi=xi*Mi/cs.sum1(xi*Mi)
        delGERT=cs.sum1(ni*cs.log(Vi))+cs.sum1((xhi@Vi)*ni)
        lnai=cs.gradient(delGERT,ni)
        ai=cs.exp(lnai)

        THFaktor=2*ni*cs.jacobian(lnai,ni)-cs.sum2(ni*cs.jacobian(lnai,ni))
        THFaktor_fun=cs.Function("TH_fun",[ni],[THFaktor])
        #print(THFaktor.shape)
        ai_fun=cs.Function("ai_fun",[ni],[ai])
        wi_fun=cs.Function("wi_fun",[ni],[wi])
        self.wi_fun=wi_fun
        self.ai_fun=ai_fun
        self.THFaktor_fun=THFaktor_fun
        #wi_fun=cs.Function("gammai_fun",[ni],[wi])
    def Isotherm(self):
        self.FH()
        ai_fun=self.ai_fun
        wi_fun=self.wi_fun
        THFaktor_fun=self.THFaktor_fun
        n=1000
        x2vec=cs.linspace(0,1,n)

        x1vec=1-x2vec
        xvec=cs.horzcat(x1vec,x2vec)
        aivec=np.asarray([ai_fun(xvec[i,:]).full().T[0] for i in range(n)])
        wivec=np.asarray([wi_fun(xvec[i,:]).full().T[0] for i in range(n)])
        THivec=np.asarray([THFaktor_fun(xvec[i,:]).full().T for i in range(n)])
        return wivec,aivec,THivec
