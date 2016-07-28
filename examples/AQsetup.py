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