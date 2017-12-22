import matplotlib.pyplot as plt
import numpy as np

dist = np.loadtxt("distribucion.txt")

plt.figure(0)
plt.plot(dist)
plt.xlabel("p^2")
plt.ylabel("Cantidad de particulas")
plt.show()
