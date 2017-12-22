import matplotlib.pyplot as plt
import numpy as np

Ec, Ep,q,x =np.loadtxt("data.txt",delimiter=";",unpack=True)
interaccion, dist1,dist2 = np.loadtxt("interaccion.txt",delimiter = ';',unpack = True)

plt.figure(1)
plt.plot(interaccion)
plt.title("Interaccion")
plt.figure(2)
plt.plot(dist1**0.5,label="Espacio posicion")
plt.plot(dist2**0.5,'r.-',label="Espacio impulsos")
plt.legend(loc='best')
plt.title("Distancia")
plt.figure(3)
plt.plot(Ec,"r-")
plt.plot(Ep+Ec,"b-")
plt.figure(4)
plt.plot(Ep,"k-")
plt.title("Energia")
"""plt.figure(4)
plt.plot(q,"b-")
plt.plot(x,"r-")
plt.xlabel("Tiempo")
plt.ylabel("Posicion")
plt.legend(["q","x"])
plt.title("Espacio fases")"""
plt.show()
