import numpy as np
import matplotlib.pyplot as plt

name = "choque.txt"

q, p, x,y,energia = np.loadtxt(name, unpack = True, delimiter = ';')
print p[0:10]
print q[0:10]
print energia[0:10]

plt.figure()
plt.plot(range(0,len(q)),q,'b')
plt.plot(range(0,len(q)),x,'g')
plt.xlim((0,100))
plt.figure()
plt.plot(range(0,len(q)),p,'r')
plt.plot(range(0,len(q)),y,'g')
plt.xlim((0,100))
plt.figure()
plt.plot(p,q)
plt.figure()
plt.plot(energia,'g')
plt.xlim((0,100))
plt.show()
