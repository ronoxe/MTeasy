##################################### Axion quadratic simple example of basic MTeasy functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
from pylab import *
import math 
import numpy as np
from scipy import interpolate 
import sympy as sym
import subprocess
from gravipy import *
location = "~/Code/MTeasy/"
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
MTS.compileName("AQ")
MTS.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyAQ

pvalue = np.zeros(3)
pvalue[0]=1.*pow(10.,-5); pvalue[1]=1.; pvalue[2]=25*pvalue[0]**2/4.0/math.pi**2;

########################### initial field values ###########################################################################################

fields = np.array([17.,.5-0.0001])
nF=MTeasyPyAQ.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
nP=MTeasyPyLH.nP() # gets number of parameters needed (useful check)

V = MTeasyPyAQ.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPyAQ.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 50.0
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPyAQ.backEvolve(t, initial,pvalue)

# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
plt.show()
#################################################################################