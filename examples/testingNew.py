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
############################### NB restart the 
nF=2
nP=5
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + 1/2.*p[1]**2 * sym.cos(p[2]/2.)**2*(f[1] - (f[0]-p[3])*sym.tan(p[2]/math.pi*sym.atan(p[4]*(f[0]-p[3]))))**2
import potentialSetup as pset 
pvalue = np.zeros(nP)
shift=231.
pvalue[0]=1.*pow(10.,-7); pvalue[1]=1.*pow(10.,-4); pvalue[2]=math.pi/10.; pvalue[3] = -100.*math.sqrt(6.) +shift; pvalue[4]=1000.*math.sqrt(3.);
#pset.potential(V,nF,nP) # writes this potential into c file when run 

subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy with new potential

############################################################################################################################################
import MTeasyPy; imp.reload(MTeasyPy)   # import recompilled module
fields = np.array([-2.-100.*math.sqrt(6.) +shift, 2.*math.tan(math.pi/20.)])
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields,pvalue) #dV=MTeasyPy.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions to be in slow roll

Nstart = 0.0
Nend = 25
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPy.backEvolve(t, initial,pvalue)