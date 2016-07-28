####################################### Setup file for the heavy field example of Langlois ###########################################
import sympy as sym
import numpy as np
import math
import sys
############################################################################################################################################

location = "/Users/David/Dropbox/MTeasyDist/MTeasy/" # this should be the location of the MTeasy folder
sys.path.append(location)  # we add this location to the python path

import MTeasySetup

### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ###
### Restart the python kernel after running this file

nF=2
nP=5
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + 1/2.*p[1]**2 * sym.cos(p[2]/2.)**2*(f[1] - (f[0]-p[3])*sym.tan(p[2]/math.pi*sym.atan(p[4]*(f[0]-p[3]))))**2

MTeasySetup.tol(1e-12,1e-18)
MTeasySetup.potential(V,nF,nP) # writes this potential into c file when run

MTeasySetup.compileName("LH") # this compiles the module with the new potential and places it in the location folder, and adds this folder to the path ready for use
############################################################################################################################################

