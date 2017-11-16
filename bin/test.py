import numpy as np
import matplotlib.pyplot as plt

name = "choque.txt"

q, p, x,y,energia,energiaC = np.loadtxt(name, unpack = True, delimiter = ';')
print p[0:10]
print y[0:10]
print q[0:10]
print x[0:10]
print energia[0:10]

n=50000

plt.figure()
plt.plot(range(0,n),q[:n],'b-')
plt.plot(range(0,n),x[:n],'g-')
plt.ylabel("Posicion")
plt.figure()
plt.plot(range(0,n),p[:n],'r-')
plt.plot(range(0,n),y[:n],'g-')
plt.ylabel("Momento")
plt.figure()
plt.plot(q[0:n],p[0:n],"-")
plt.ylabel("Momento")
plt.xlabel("Posicion")
plt.figure()
plt.plot(energia[0:n],'g')
"""plt.plot(energiaC[0:n],'r')"""
plt.ylabel("Energia")
plt.ylim((0.4,.6))
"""plt.legend(["Total", "Cinetica"])"""
plt.show()
