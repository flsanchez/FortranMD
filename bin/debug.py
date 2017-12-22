import numpy as np
import matplotlib.pyplot as plt
import sys

name = sys.argv[1]
n = 8
ini = int(sys.argv[2])
fin = int(sys.argv[3])

q = np.zeros((n,fin-ini,3))
p = np.zeros((n,fin-ini,3))
x = np.zeros((n,fin-ini,3))
y = np.zeros((n,fin-ini,3))

with open(name,'r') as f:
	data = f.readlines()
	for step in range(fin-ini): # para cada paso temporal
		line = data[step+ini].rstrip(';\n').split(';')
		for part in range(n): # para cada particula
			aux = np.array([float(coord) for coord in line[part].rstrip(',').split(',')]) #las 12 coordenadas
			q[part,step,:] = aux[0:3]
			y[part,step,:] = aux[3:6]
			x[part,step,:] = aux[6:9]
			p[part,step,:] = aux[9:12]
	
etiq = ['x','y','z']
for part in range(n):
	plt.figure(part)
	plt.suptitle("Particula {0}".format(part))
	for i in range(3):
		plt.subplot(3,2,2*i+1)
		plt.plot(q[part,:,i],'b',label='q')
		plt.plot(x[part,:,i],'r',label='x')
		plt.title("Posicion {0}".format(etiq[i]))
		plt.legend(loc='best')
		plt.subplot(3,2,2*i+2)
		plt.plot(p[part,:,i],'b',label='p')
		plt.plot(y[part,:,i],'r',label='y')
		plt.title("Momento {0}".format(etiq[i]))
		plt.legend(loc='best')
plt.show()
	
		
	
