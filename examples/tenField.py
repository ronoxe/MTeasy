####################################### ten field quadratic example of basic MTeasy functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
from pylab import *
import math 
import numpy as np
from scipy import interpolate 
import sympy as sym
import subprocess
############################################################################################################################################


#### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################### NB restart the kernel before running if previously using a differnet potential. ############################
nF=10
nP=10
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V=1./2.*(p[0]**2*f[0]**2 + p[1]**2*f[1]**2 + p[2]**2*f[2]**2 + p[3]**2*f[3]**2 + p[4]**2*f[4]**2 + p[5]**2*f[5]**2
  + p[6]**2*f[6]**2 + p[7]**2*f[7]**2 + p[8]**2*f[8]**2 + p[9]**2*f[9]**2)
import potentialSetup as p 
pvalue = np.zeros(3)
pvalue = np.array([2.45391E-10, 1.72583E-10, 1.44638E-10, 1.06326E-10, 8.53891E-11, 7.74336E-11, 5.1783E-11, 4.34577E-11, 2.37254E-11, 1.77517E-11])
p.potential(V,nF,nP,pvalue) # writes this potential into c file when run 
subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy
############################################################################################################################################


import MTeasyPy; imp.reload(MTeasyPy)   # import recompilled module


########################### initial field values ###########################################################################################
fields = np.array([sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), 
                   sqrt(142.0/5.0), sqrt(142.0/5.0), sqrt(142.0/5.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields) # calculate potential from some initial conditions
dV=MTeasyPy.dV(fields) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/sqrt(3.*V))) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 65.0
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPy.backEvolve(t, initial)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'r')
############################################################################################################################################


############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
Nexit = 13.0
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-12),der=0)
k = np.exp(Nexit) *  MTeasyPy.H(backExit); 
# other scales can then be defined wrt to k
############################################################################################################################################


################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds sub-horizon evolution using another spline
Nstart = Nexit - 4.0

backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-12),der=0)

t=np.linspace(Nstart,Nend, 1000)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an estimate for ns
twoPt = MTeasyPy.sigEvolve(t, k, backExitMinus, 1)
zz=twoPt[:,1]
zz1=zz[-1]
twoPt=MTeasyPy.sigEvolve(t, k+.1*k, backExitMinus, 1)
zz=twoPt[:,1]
zz2=zz[-1]
n_s = (log(zz2)-log(zz1))/(log(k+.1*k)-log(k))+4;
print n_s
fig2 = plt.figure(2)
plt.plot(t, zz,'r')
############################################################################################################################################


###################################### example bispectrum run ##############################################################################
# set three scales in FLS manner (using alpha beta notation)
alpha=0
beta =1/3.

k1 = k/2 - beta*k/2.
k2 = k/4*(1+alpha+beta)
k3 = k/4*(1-alpha+beta)

# find initial conditions that allows longest mode 4 e-folds inside horizon
Nstart = Nexit - max(log(k/k1),log(k/k2),log(k/k3))-4.0
backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-12),der=0)

start_time = time.time()
# run solver for this triangle
t=np.linspace(Nstart,Nend, 1000)
threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus,1)
zzz= threePt[:,:5]

fig3 = plt.figure(3)
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
print  time.time() - start_time
############################################################################################################################################

plt.show(fig1)
plt.show(fig2)
plt.show(fig3)

