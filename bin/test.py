import numpy as np
import matplotlib.pyplot as plt

name = "choque.txt"

q, p, x,y,energia,energiaC,H = np.loadtxt(name, unpack = True, delimiter = ';')
print p[0:10]
print y[0:10]
print q[0:10]
print x[0:10]
print energia[0:10]

n=50000
"""
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
plt.plot(energiaC[0:n],'r')
plt.ylabel("Energia")
#plt.ylim((0.4,.6))
plt.legend(["Total", "Cinetica"])
plt.figure()
plt.plot(H,'b')
plt.show()"""

plt.figure(1)
plt.title("Energia Total")
plt.plot(energia,'g')
plt.grid()
plt.figure(2)
plt.title("Energia Cinetica")
plt.plot(energiaC,'r')
plt.grid()
plt.figure(3)
plt.title("Energia Potencial")
plt.plot(energia-energiaC,'b')
plt.grid()
plt.figure(4)
plt.show()
