import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.optimize as sco

N = int(sys.argv[1])
temp = sys.argv[2]
Nbins=10
if len(sys.argv) == 4:
    Nbins = int(sys.argv[3])
name = "distribucion_" + str(N)+ "_" + "%1.4f" %float(temp) + ".txt"


data = np.loadtxt(name, unpack = True, delimiter=";")
n = len(data)/(3*N+2)       # Cantidad de muestras

T = []
T_corr = []
p = np.array([])
for i in range(0,n):
  T.append(data[(3*N+2)*i])
  T_corr.append(data[(3*N+2)*i+1])
  p = np.concatenate((p,data[(3*N+2)*i+2:(3*N+2)*(i+1)]))

p = np.abs(p)
T = np.array(T)
T_corr = np.array(T_corr)
fermi = lambda p,mu,To: 2.0/(np.exp((0.5*p**2-mu)/To)+1.0)
#boltzmann = lambda p,To: np.sqrt(2*np.pi/(1.88**2*To))*np.exp(-0.5*p**2/To)*N/(2*512)
boltzmann = lambda p,To,C: C*np.exp(-0.5*p**2/To)

print len(p)

plt.figure(0)
H,bins = np.histogram(p,Nbins)
bins = (bins[0:len(H)]+bins[1:len(H)+1])/2
H = H/float(3*n)            # Hay 3n impulsos muestreados, me interesa saber cuantas particulas tenian impulso en cada bin
plt.plot(bins,H, "b-")
paramsf,coso = sco.curve_fit(fermi,bins,H)
paramsb,coso = sco.curve_fit(boltzmann,bins,H,bounds=(0,np.inf))
F = fermi(bins,paramsf[0],paramsf[1])
B = boltzmann(bins,paramsb[0],paramsb[1])
#z = np.exp(params[0]/params[1])
To = np.mean(T/(1.5*N)+T_corr)
LOT = np.sqrt(2*np.pi/(1.88**2*To))
plt.plot(bins,F, "r-.")
plt.plot(bins,B, "k--")
plt.legend(("Datos", "Ajuste Fermi", "Ajuste Boltzmann"))
print paramsf, paramsb, sum(F),sum(B), LOT#, z

# plt.figure(1)
# plt.hist(p,bins=Nbins)
# plt.xlabel("Impulso")
# plt.ylabel("Ocurrencias")
plt.figure(2)
plt.plot(T/(1.5*N)+T_corr, "r-.")
plt.plot(T/(1.5*N), "b--")
plt.legend(["Temperatura virial", "Temperatura ideal"])
print "T_ideal=", np.mean(T/(1.5*N)), "+-", np.std(T)/(1.5*N)
print "T_virial=", np.mean(T/(1.5*N)+T_corr), "+-", np.std(T/(1.5*N)+T_corr)
plt.show()
