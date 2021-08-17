import matplotlib.pyplot as plt
import numpy as np
from Pylomer.Kwei import Tggt
w1=np.linspace(0,1,100)
w2=1-w1
w3=np.zeros_like(w1)
rhoi=np.asarray([1320,997,1180])
Tgi=np.asarray([317.6,136,420])
q=np.asarray([-600,0,0])
#w1=np.linspace(0,1,100)
#w2f=np.linspace(0,1,100)
#w2=w2f*(1-w1)
#w3=(1-w2f)*(1-w1)
wi=np.asarray([w1,w2,w3])

rhoi[0]=(Tgi[1]*rhoi[1])/(0.11*Tgi[0])
Tgw1=Tggt(wi.T,rhoi,Tgi,q=q)
Tgw2=Tggt(wi.T,rhoi,Tgi,q=0)
plt.plot(w1,Tgw1)
plt.plot(w1,Tgw2)
plt.show()
