##################################### Axion quadratic simple example of basic MTeasy functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
from pylab import *
import math 
import numpy as np
from scipy import interpolate as interp
import sympy as sym
import subprocess
from gravipy import *
location = "/home/jwr/Code/MTeasy"
sys.path.append(location)
## Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################# NB restart the 
nF=2
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + p[2] * (1-sym.cos(2*math.pi * f[1] / p[1]))
G=Matrix([[1,0],[0,1]])
import MTeasySetup as MTS
MTS.tol(1e-12,1e-18)
MTS.potential(G,V,nF,nP)
MTS.fieldmetric(G,nF)
MTS.compileName("AQ")
MTS.pathSet()  # his add sets the other paths that MTeasy uses
import MTeasyPyAQ

pvalue = np.zeros(3)
pvalue[0]=1.*pow(10.,-5); pvalue[1]=1.; pvalue[2]=25*pvalue[0]**2/4.0/math.pi**2;

########################## initial field values ###########################################################################################

fields = np.array([17.,.5-0.0001])
nF=MTeasyPyAQ.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
nP=MTeasyPyAQ.nP() # gets number of parameters needed (useful check)

V = MTeasyPyAQ.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPyAQ.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll


################################ run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 50.0
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPyAQ.backEvolve(t, initial,pvalue)

####plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
###############################################################################

############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
Nexit = 13.0
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interp.splev(Nexit,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)
k = np.exp(Nexit) *  MTeasyPyAQ.H(backExit,pvalue);
################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds sub-horizon evolution using another spline
Nstart = Nexit - 4.0

backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interp.splev(Nstart,interp.splrep(back[:,0], back[:,i], s=1e-12),der=0)
print backExitMinus, 'BKexit'
t=np.linspace(Nstart,Nend, 1000)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an estimate for ns
twoPt = MTeasyPyAQ.sigEvolve(t, k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz1=zz[-1]
twoPt=MTeasyPyAQ.sigEvolve(t, k+.1*k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz2=zz[-1]
n_s = (log(zz2)-log(zz1))/(log(k+.1*k)-log(k))+4;
print n_s
fig2 = plt.figure(2)
plt.plot(t, zz,'r')

plt.show()
