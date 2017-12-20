import matplotlib.pyplot as plt
import numpy as np
import sys

name = sys.argv[1]

p1,p2,p3,y1,y2,y3=np.loadtxt(name,delimiter=';',unpack=True)

plt.plot(p1)
plt.plot(y1,'r')
plt.xlim(0,5000)
plt.show()
