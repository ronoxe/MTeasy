####################################### Axion quadratic simple example of basic MTeasy functions ###########################################
import sympy as sym
import subprocess
import math
import numpy as np
############################################################################################################################################


#### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################### NB restart the 
nF=2
nP=5
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + 1/2.*p[1]**2 * sym.cos(p[2]/2.)**2*(f[1] - (f[0]-p[3])*sym.tan(p[2]/math.pi*sym.atan(p[4]*(f[0]-p[3]))))**2
import potentialSetup as pset 

pset.potential(V,nF,nP) # writes this potential into c file when run

subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy with this new potential

############################################################################################################################################
