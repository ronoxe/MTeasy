from matplotlib import pyplot as plt
import time
from pylab import *
import math
import numpy as np
from scipy import interpolate
import MTeasyPy
import subprocess
import pylab

# This example is intended to be run for the double-quadratic potential defined in 
# the potentialSetup.py file. Ensure the potential is the one compiled into the MTeasyPy module
# It calculates and plots a slice through the bispectrum of this model
# It is recommended you look through example.py first.

# initial field values
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(np.array([24.3,0.48])) # calculate potential from some initial conditions
dV = np.zeros(nF) # calculat1e derivatives of potential
MTeasyPy.dV(np.array([24.3,0.48]), dV )

print dV
# calculate derivatives of potential (changes dV to derivatives)
initial = np.array([24.3,0.48,-dV[0]/sqrt(3.*V),-dV[1]/sqrt(3.*V)])
Nstart = 0.0
Nend = 75.0
t=np.linspace(Nstart, Nend, 1000)
MTeasyPy.backEvolve(t, initial)
back = np.genfromtxt('../../runData/back.dat', dtype=None)

fig1 = plt.figure()
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'r')

# set a pivot scale which exits after certain time 
# in this example we treat this scale as k_t
Nexit = 10.0
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
ki = np.exp(Nexit) *  MTeasyPy.H(backExit); 
# other scales can then be defined wrt to k



start_time = time.time()

# example bispectrum run
# runs through many vaues of alpa and beta (throwing away those not within range, and calculates bispectrum at end of run)
# side is number of alphas and betas on each side it tries
biA=np.array([])
spA = np.array([])
kA= np.array([])
k=ki
for l in range(0,20):
        print l
        alpha=0
        beta = 1/3.        
        k1 = k/2 - beta*k/2.
        k2 = k/4*(1+alpha+beta)
        k3 = k/4*(1-alpha+beta)
        
        Nstart = Nexit - max(log(ki/k1),log(ki/k2),log(ki/k3))-4.0
        backExitMinus = np.zeros(2*nF)
        for i in range (1,2*nF+1):
            backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

# run solver for this triangle
        t = np.linspace(Nstart,Nend,5) 
        MTeasyPy.alphaEvolve(t ,k1,k2,k3, backExitMinus,0)
        zzz = np.genfromtxt('../../runData/zzz.dat', dtype=None)         
         
        bi=5.0/6.0*np.divide(zzz[-1,4], (np.multiply(zzz[-1,2],zzz[-1,3])+np.multiply(zzz[-1,1],zzz[-1,2]) + np.multiply(zzz[-1,1],zzz[-1,3])))
        bi2=5.0/6.0*np.divide(zzz[-1,4], (1/k1**2 * 1/k2**2 + 1/k1**2 * 1/k3**2 + 1/k2**2 * 1/k3**2))
        biA=np.append(biA,bi)
        spA=np.append(spA,zzz[-1,1])        
        kA=np.append(kA,k1)
        k=k+0.7*k;
             
        
#fig3 = plt.figure(3)
#plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
print  time.time() - start_time

plt.show(fig1)

fig2 = plt.figure(2)
plt.plot(log(kA/kA[0]), biA, 'r',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
pylab.ylabel(r'$fNL$', fontsize=20)
pylab.xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("oscFnl.png")



fig3 = plt.figure(3)
plt.plot(log(kA/kA[0]),log(spA/spA[0]*kA**3/kA[0]**3), 'g',linewidth = 2)
title('Power spectrum',fontsize=15)
grid(True)
plt.legend(fontsize=15)
pylab.ylabel(r'${\cal P}$', fontsize=20)
pylab.xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("oscP.png")


fig4 = plt.figure(4)
plt.plot(log(kA[1:]/kA[0]),diff(log(spA/spA[0]*kA**3/kA[0]**3))/diff(log(kA/kA[0])), 'g',linewidth = 2)
title('Spectral index',fontsize=15)
grid(True)
plt.legend(fontsize=15)
pylab.xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("oscNS.png")

plt.show(fig2)
plt.show(fig3)
plt.show(fig4)