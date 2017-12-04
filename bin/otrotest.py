import matplotlib.pyplot as plt
import numpy as np

Ec, Ep,q,x =np.loadtxt("data.txt",delimiter=";",unpack=True)
"""
plt.figure(1)
plt.plot(H)
plt.title("-Entropia")"""
plt.figure(2)
plt.plot(Ec,"r-")
plt.plot(Ep+Ec,"b-")
plt.figure(3)
plt.plot(Ep,"k-")
plt.title("Energia")
plt.figure(4)
plt.plot(q,"b-")
plt.plot(x,"r-")
plt.xlabel("Tiempo")
plt.ylabel("Posicion")
plt.legend(["q","x"])
plt.title("Espacio fases")
plt.show()
