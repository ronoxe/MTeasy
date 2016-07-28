from matplotlib import pyplot as plt
from mayavi.mlab import *

import time
from pylab import *
import math
import numpy as np
import subprocess
from scipy import interpolate
import pylab
import MTeasyPy
reload(MTeasyPy)
# This example is intended to be run for the heavy model




nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(np.array([sqrt(500./2),0.00000000000000001])) # calculate potential from some initial conditions
dV = np.zeros(nF) # calculat1e derivatives of potential
MTeasyPy.dV(np.array([sqrt(500./2),0.0000000000000001]),dV) # calculate derivatives of potential (changes dV to derivatives)

initial = np.array([sqrt(500./2),0.00000000000000001,-dV[0]/sqrt(3.*V),0.0]) # set initial conditions to be in slow roll
Nstart = 0.0
Nend = 30.0
t=np.linspace(Nstart, Nend, 1000)
# solve a fiducial background run
MTeasyPy.backEvolve(t, initial)
back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background

fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'r')



# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t

kOut=np.array([])
side = 3
nk=4
biP= np.zeros((side,side,nk))


PhiExit = 14.2
backExit = np.zeros(2*nF)
Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[::-1,1], back[::-1,0], s=1e-15),der=0)

for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
    ks = np.exp(Nexit) *  MTeasyPy.H(backExit); 

PhiExit = 15.0
backExit = np.zeros(2*nF)
Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[::-1,1], back[::-1,0], s=1e-15),der=0)

for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
    ke = np.exp(Nexit) *  MTeasyPy.H(backExit); 

lkOut = np.linspace(0.0,log(ks/ke),num=nk)

for ii in range(0,nk):
    print ii
    start_time = time.time()

    k=exp(lkOut[ii]) *ke
# other scales can then be defined wrt to k
    kOut=np.append(kOut,k)
    
#example 2pt run 
    step =1.0/side
#np.zeros([side,side])
    alphaA = np.array([])
    num=0    
    for l in range(0,side):
        alpha =  -1+step + l*2*step  
        alphaA=np.append(alphaA,alpha)    
        betaA = np.array([]) 
        for j in range(0,side):
            beta =j*step
            betaA=np.append(betaA,beta)               
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
                bi2=5.0/6.0*np.divide(zzz[-1,4], (1/k1**2 * 1/k2**2 + 1/k1**2 * 1/k3**2 + 1/k2**2 * 1/k3**2))
                biP[l,j,ii]=bi
                num=num+1
            else:  
                biP[l,j,ii]=0.   
    print  time.time() - start_time
    print num
plt.show(fig1)

betaP, alphaP, kP = np.mgrid[betaA[0]:betaA[-1]:size(betaA)*1j, alphaA[0]:alphaA[-1]:size(alphaA)*1j, lkOut[0]:lkOut[-1]:size(lkOut)*1j]
obj = contour3d(betaP,alphaP,kP,biP, contours=4, transparent=True)
print  time.time() - start_time
