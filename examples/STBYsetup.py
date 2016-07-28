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
## Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################# NB restart the 
nF=1
nP=2
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= p[0]*(1-sym.exp(-1.0 * p[1]*f[0]))
import MTeasySetup as MTS
MTS.tol(1e-12,1e-18)
MTS.potential(V,nF,nP)
MTS.compileName("ST")
MTS.pathSet()  # his add sets the other paths that MTeasy uses
import MTeasyPyST


