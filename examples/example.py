import numpy as np
from Pylomer.FloryHuggins import FloryFilm
import matplotlib.pyplot as plt
T=298.15
sol={"Mi":0.01815,"rho0":997}
pol={"Mi":0.721,"rho0":1150}
chi=np.asarray([2.272])
pure=[sol,pol]
FilmFlory=FloryFilm(pure,chi)
wivec,aivec,THivec=FilmFlory.Isotherm()
RH=aivec[:,0]
w1vec=wivec[:,0]
fig,ax=plt.subplots()
fig1,ax1=plt.subplots()
ax.plot(RH,w1vec)
ax1.plot(w1vec,THivec[:,0,0])

plt.show()
