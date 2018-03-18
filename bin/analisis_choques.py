import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])

for i in range(0,N):
  if(i<10):
      name = "choque_ %d.txt" %(i)
  else:
      name = "choque_%d.txt" %(i)

  q, p = np.loadtxt(name, unpack = True, delimiter = ';')

  plt.figure(1)
  plt.plot(q,p,"b-")

plt.figure(1)
plt.ylabel("Momento")
plt.xlabel("Posicion")
plt.axis([-3,3,-4,4])
plt.show()
