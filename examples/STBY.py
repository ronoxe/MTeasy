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
location = "/home/jwr/Code/MTeasy"
sys.path.append(location)

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyST
pvalue = np.zeros(2)
M=1.*pow(10.,-5)
pvalue[0]=3.0/4.0 * M**2; pvalue[1]=math.sqrt(2.0/3.0);

########################## initial field values ###########################################################################################
 
fields = np.array([4.80])
nF=MTeasyPyST.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
nP=MTeasyPyST.nP() # gets number of parameters needed (useful check)

V = MTeasyPyST.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPyST.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll


################################ run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 70.0
t=np.linspace(Nstart, Nend, 100)
back = MTeasyPyST.backEvolve(t, initial,pvalue)

####plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,1 ], 'g')
###############################################################################
plt.show()
############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t

Nexit = 13.3
backExit = np.zeros(2*nF)

for i in range (1,2*nF+1):
    backExit[i-1] = interp.splev(Nexit,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)
ke = np.exp(Nexit) *  MTeasyPyST.H(backExit,pvalue);
print ke, 'end'

Nexit = 13.1
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interp.splev(Nexit,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)
k = np.exp(Nexit) *  MTeasyPyST.H(backExit,pvalue);
print k, 'start'
################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds sub-horizon evolution using another spline
Nstart = Nexit - 4.0

backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interp.splev(Nstart,interp.splrep(back[:,0], back[:,i], s=1e-12),der=0)

t=np.linspace(Nstart,Nend, 2)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an estimate for ns

twoPt = MTeasyPyST.sigEvolve(t, k, backExitMinus,pvalue, 1)

zz=twoPt[:,1]
zz1=zz[-1]
twoPt=MTeasyPyST.sigEvolve(t, ke, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz2=zz[-1]

n_s = (math.log(zz2)-math.log(zz1))/(math.log(k+.1*k)-math.log(k))+4;
print 'n_s = ', n_s

twoPt=MTeasyPyST.sigEvolve(t,1.0*k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
pp= [1.0]
kmode=100
kkk=np.linspace(k,ke,kmode)

#Calculates the power Spectrum
for kk in range(kmode-1):
	tp = MTeasyPyST.sigEvolve(t, kkk[kk], backExitMinus,pvalue, 1)
	zz3=tp[:,1]
	pp = np.append(pp,[zz3[-1]/zz1])

	print 'k mode # ',kk, ' '
fig2 = plt.figure(2)
plt.plot(t, zz,'r')
print pp
fig3= plt.figure(3)
kkk=kkk/k
kp=np.log(kkk)
pp=np.log(abs(pp))
np.savetxt('PSout.out',(kp,pp))
plt.plot(kp[1:],pp[1:])
plt.show()








