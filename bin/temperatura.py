import numpy as np
import matplotlib.pyplot as plt

name = "temperaturas_64.txt"

T, Tb, Tv = np.loadtxt(name, unpack = True, delimiter=";")

m,b = np.polyfit(Tb,Tv,1)
R = np.mean((Tv-np.mean(Tv))*(Tb-np.mean(Tb)))/(np.std(Tb)*np.std(Tv))

print np.mean(T), "+-", np.std(T)
print T,Tb, Tv

plt.plot(Tb,Tv,"bo")
plt.plot(Tb,m*Tb+b,"r-")
plt.legend(["Puntos", "m="+str(m)+"\nb="+str(b)+"\nR="+str(R)])
plt.xlabel("Temperatura Boltzmann")
plt.ylabel("Temperatura Virial")
plt.show()
