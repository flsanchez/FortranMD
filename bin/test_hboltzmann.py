import numpy as np
import matplotlib.pyplot as plt

name = 'data.txt'

h,Ec,Ep = np.loadtxt(name,delimiter=';',unpack=True)

print Ec[0:10]

plt.figure()
plt.plot(h)
plt.figure()
plt.plot(Ec+Ep,'r')
plt.show()
