from matplotlib import pyplot as plt
import time
from pylab import *
import math
import numpy as np
import subprocess
from scipy import interpolate
import MTeasyPy
reload(MTeasyPy)
# This example is intended to be run for the axion-quadratic potential defined in 
# the potentialSetup.py file. Ensure the potential is the one compiled into the MTeasyPy module.
# Its purpose is to show the use of the key functions provided by MTeasyPy

nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(np.array([24.3,0.46])) # calculate potential from some initial conditions
dV = np.zeros(nF) # calculat1e derivatives of potential
MTeasyPy.dV(np.array([24.3,0.46]), dV )

print dV
# calculate derivatives of potential (changes dV to derivatives)
initial = np.array([24.3,0.46,-dV[0]/sqrt(3.*V),-dV[1]/sqrt(3.*V)])
Nstart = 0.0
Nend = 75.0
t=np.linspace(Nstart, Nend, 1000)
# solve a fiducial background run
MTeasyPy.backEvolve(t, initial)
back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'r')
plt.plot(back[:,0], back[:,4 ], 'r')


ep1 = np.zeros(size(back[:,1]))
for ii in range (1, size(back[:,1]) ):
    MTeasyPy.dV(np.array([back[ii,1],back[ii,2]]), dV ) 
    MTeasyPy.V (np.array([back[ii,1],back[ii,2]]))   
    ep1[ii] = 1/2.*(dV[0]*dV[0])/V**2
eta1 = np.diff(ep1)/np.diff(back[:,0])/sqrt(ep1[1:])
    
fig5 = plt.figure(5)
plt.plot(back[:,0], ep1, 'r')
plt.plot(back[1:,0], eta1, 'r')


# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
Nexit = 10.0
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
k = np.exp(Nexit) *  MTeasyPy.H(backExit); 
# other scales can then be defined wrt to k


#example 2pt run 

#find conditions for 4 e-folds sub-horizon evolution using another spline
Nstart = Nexit - 4.0
backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

t=np.linspace(Nstart,Nend, 1000)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum
MTeasyPy.sigEvolve(t, k, backExitMinus, 1)
zz = np.genfromtxt('../../runData/zz.dat', dtype=None)
zz1=zz[-1,1]
MTeasyPy.sigEvolve(t, k+.1*k, backExitMinus, 1)
zz = np.genfromtxt('../../runData/zz.dat', dtype=None)
zz2=zz[-1,1]
n_s = (log(zz2)-log(zz1))/(log(k+.1*k)-log(k))+4;
print n_s
fig2 = plt.figure(2)
plt.plot(zz[:,0], zz[:,1],'r')


# example bispectrum run
# set three scales in FLS manner
#alpha = 1./6 
alpha=0.
#beta = 1./2
beta = 1/3.
#beta = 0.95

k1 = k/2 - beta*k/2.
k2 = k/4*(1+alpha+beta)
k3 = k/4*(1-alpha+beta)

k1 = 10 *k
k2 = 10* k
k3= 10* k

# find initial conditions that allows longest mode 4 e-folds inside horizon
Nstart = Nexit - max(log(k/k1),log(k/k2),log(k/k3))-4.0
backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

start_time = time.time()

# run solver for this triangle
t=np.linspace(Nstart,Nend, 1000)
MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus,1)
zzz = np.genfromtxt('../../runData/zzz.dat', dtype=None)

fig3 = plt.figure(3)
P=4.*10**(-10)
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
#plt.plot(zzz[500:,0], k1**2*k2**2*k3**2/P/(2*pi)**7*zzz[500:,4])

print  time.time() - start_time

plt.show(fig1)

plt.show(fig2)

plt.show(fig3)

plt.show(fig5)