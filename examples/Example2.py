from matplotlib import pyplot as plt
import time
from pylab import *
import math
import numpy as np
from scipy import interpolate
import MTeasyPy
import subprocess


# This example is intended to be run for the axion-quadratic potential defined in 
# the potentialSetup.py file. Ensure the potential is the one compiled into the MTeasyPy module
# It calculates and plots a slice through the bispectrum of this model
# It is recommended you look through example.py first.

# initial field values
nF=MTeasyPy.nF()
initialFields = np.array([24.2,0.5-0.001])
V = MTeasyPy.V(initialFields)
dV = np.zeros(nF) # calculate derivatives of potential
MTeasyPy.dV(initialFields,dV) # calculate derivatives of potential (changes dV to derivatives)
initial = np.array([24.2,0.5-0.001,-dV[0]/sqrt(3.*V),-dV[1]/sqrt(3.*V)]) # set initial conditions to be in slow roll

Nend = 75.0
# solve a fiducial background run
t=np.linspace(0.0, Nend, 1000)
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
k = np.exp(Nexit) *  MTeasyPy.H(backExit); 
# other scales can then be defined wrt to k



start_time = time.time()

# example bispectrum run
# runs through many vaues of alpa and beta (throwing away those not within range, and calculates bispectrum at end of run)
# side is number of alphas and betas on each side it tries
side = 30
step =1.0/side
biA= np.array([])
#np.zeros([side,side])
alphaA = np.array([])
betaA = np.array([]) 
for l in range(0,side):
      print l
      alpha =  -1+step + l*2*step    
      for j in range(0,side):
         beta =j*step
         if alpha>-(1-beta) and alpha < 1-beta :            
             k1 = k/2 - beta*k/2.
             k2 = k/4*(1+alpha+beta)
             k3 = k/4*(1-alpha+beta)

# find conditions that allows longest mode 4 e-folds inside horizon
             Nstart = Nexit - max(log(k/k1),log(k/k2),log(k/k3))-4.0
             backExitMinus = np.zeros(2*nF)
             for i in range (1,2*nF+1):
                 backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

# run solver for this triangle
             t = np.linspace(Nstart,Nend,5) 
             MTeasyPy.alphaEvolve(t ,k1,k2,k3, backExitMinus,0)
             zzz = np.genfromtxt('../../runData/zzz.dat', dtype=None)         
         
             bi=5.0/6.0*np.divide(zzz[-1,4], (np.multiply(zzz[-1,2],zzz[-1,3])+np.multiply(zzz[-1,1],zzz[-1,2]) + np.multiply(zzz[-1,1],zzz[-1,3])))
             bi2=5.0/6.0*np.divide(zzz[-1,4], (1/k1**3 * 1/k2**3 + 1/k1**3 * 1/k3**3 + 1/k2**3 * 1/k3**3))
             betaA=np.append(betaA,beta)
             alphaA=np.append(alphaA,alpha)
             biA=np.append(biA,bi)
#fig3 = plt.figure(3)
#plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
print  time.time() - start_time

plt.show(fig1)


from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
# note this: you can skip rows!
X = alphaA
Y = betaA
Z = biA

xi = np.linspace(X.min(),X.max(),100)
yi = np.linspace(Y.min(),Y.max(),100)
# VERY IMPORTANT, to tell matplotlib how is your data organized
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

CS = plt.contour(xi,yi,zi,15,linewidths=1.5,color='k')
ax = fig.add_subplot(1, 2, 2, projection='3d')

xig, yig = np.meshgrid(xi, yi)

surf = ax.plot_surface(xig, yig, zi,cmap=cm.coolwarm,
        linewidth=0.5)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show(fig)

#call biPlot.py which makes a better 3D plot of bispectrum as function of alpha beta(file is in useful directory, ensure that current directory is in path)
import biPlot

