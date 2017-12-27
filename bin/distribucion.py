import numpy as np
import matplotlib.pyplot as plt

name = "distribucion.txt"

P, proba = np.loadtxt(name, unpack = True, delimiter = ';')

plt.plot(P,proba,"o")
plt.xlabel("p^2")
plt.ylabel("Probabilidad")
plt.title("Distribucion de impulsos")
plt.show()
