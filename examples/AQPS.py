from matplotlib import pyplot as plt
import time
import math 
import numpy as np
import sys 
from scipy import interpolate as interp

location = "/home/jwr/Code/MTeasy"
sys.path.append(location)

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyAQ;  # import module  
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

t=np.linspace(Nstart,Nend, 1000)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an estimate for ns
print backExitMinus, 'BKexit'

print k, 'k'

twoPt = MTeasyPyAQ.sigEvolve(t, k, backExitMinus,pvalue, 1)

zz=twoPt[:,1]
zz1=zz[-1]
twoPt=MTeasyPyAQ.sigEvolve(t, k+.1*k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz2=zz[-1]
n_s = (math.log(zz2)-math.log(zz1))/(math.log(k+.1*k)-math.log(k))+4;
print n_s
fig2 = plt.figure(2)
plt.plot(t, zz,'r')

plt.show()