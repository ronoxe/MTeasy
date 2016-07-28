####################################### Axion quadratic simple example of basic MTeasy functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
#from pylab import *
import math 
import numpy as np
from mayavi.mlab import *
from scipy import interpolate 
import sympy as sym
import subprocess
############################################################################################################################################

#### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################### NB restart the 
nF=2
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + p[2] * (1-sym.cos(2*math.pi * f[1] / p[1]))
import potentialSetup as p 
pvalue = np.zeros(3)
p.potential(V,nF,nP) # writes this potential into c file when run 
subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy with new potential
############################################################################################################################################

import MTeasyPy; imp.reload(MTeasyPy)   # import recompilled module

########################### initial field values ###########################################################################################
fields = np.array([17.,.5-0.0001])
pvalue[0]=1.*pow(10.,-5); pvalue[1]=1.; pvalue[2]=25*pvalue[0]**2/4.0/math.pi**2;
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPy.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 50.0
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPy.backEvolve(t, initial,pvalue)
############################################################################################################################################

NexitS = 13.0
Npoints = 30
Nstep = .1/3
side = 30

k1Out = np.array([])
k2Out = np.array([])
k3Out = np.array([])
biOut = np.array([])
biROut = np.array([])


############################################################################################################################################
for ll in range(0,Npoints) :
    Nexit = NexitS+Nstep*ll
    backExit = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        k = np.exp(Nexit) *  MTeasyPy.H(backExit,pvalue); 
# other scales can then be defined wrt to k
############################################################################################################################################


# example bispectrum run
# runs through many vaues of alpa and beta (throwing away those not within range, and calculates bispectrum at end of run)
# side is number of alphas and betas on each side it tries
    step =1.0/side

#np.zeros([side,side])
    
    for l in range(0,side+1):
        alpha =  -1.0 + l*2*step          
        for j in range(0,side+1):
            print "ll", ll, "l " ,  l,   " j " , j  
            beta =j*step
            if alpha>-(1-beta) and alpha < 1-beta :            
                k1 = k/2 - beta*k/2.
                k2 = k/4*(1+alpha+beta)
                k3 = k/4*(1-alpha+beta)
                #print k1, k2, k3
                # find conditions that allows longest mode 4 e-folds inside horizon
                NstartN = Nexit - max(math.log(k/k1),math.log(k/k2),math.log(k/k3))-4.0
                backExitMinus = np.zeros(2*nF)
                for i in range (1,2*nF+1):
                    backExitMinus[i-1] = interpolate.splev(NstartN,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

# run solver for this triangle
                t=np.linspace(NstartN,Nend, 5)
                threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus, pvalue,0)
                zzz= threePt[:,:5]
                #dim= np.abs(zzz[:,4]*(k1*k2*k3)**2)
                biR=5.0/6.0*np.divide(zzz[-1,4], (np.multiply(zzz[-1,2],zzz[-1,3])+np.multiply(zzz[-1,1],zzz[-1,2]) + np.multiply(zzz[-1,1],zzz[-1,3])))
                biROut=np.append(biROut,biR)
                biOut=np.append(biOut,zzz[-1,4])
                k1Out = np.append(k1Out, k1)
                k2Out = np.append(k2Out, k2)
                k3Out = np.append(k3Out, k3)
#fig3 = plt.                 
#plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
pts = 200j
ks = np.hstack((k1Out.reshape((len(k1Out),1)),k2Out.reshape((len(k2Out),1)),k3Out.reshape((len(k3Out),1))))
X,Y,Z = np.mgrid[min(k1Out):max(k2Out):pts, min(k2Out):max(k2Out):pts, min(k3Out):max(k3Out):pts]

#F = griddata(ks, biROut, (X,Y,Z))
S=1.0e22*biOut*(k1Out**2*k2Out**2*k3Out**2)
F=griddata(ks,S, (X,Y,Z))
axes()
contour3d(F,contours=10,opacity=.5 )

