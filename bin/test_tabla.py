import matplotlib.pyplot as plt
import numpy as np

name = "tabla.txt"

f = open(name, "r")
f.readline()
data = f.readline()
print data[0:1000]
tabla = data.split(" ")
tabla = [float(tabla[i]) for i in range(0,len(tabla)-1)]

print len(tabla)
plt.plot(tabla,"r")
plt.show()
