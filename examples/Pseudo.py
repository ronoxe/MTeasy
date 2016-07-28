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
nF=3
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + 1./2. * p[1]**2 * f[1]**2 + 1./2. * p[2]**2 * f[2]**2
R=(0.9)/((sym.cosh(2.0*((1.0*f[0])-7.0)/0.12))**(2.0))
G=Matrix([[1,R,0],[R,1,0],[0,0,1]])
import MTeasySetup as MTS
MTS.tol(1e-12,1e-18)
MTS.potential(G,V,nF,nP)
MTS.fieldmetric(G,nF)
MTS.compileName("PS")
MTS.pathSet()  # his add sets the other paths that MTeasy uses
import MTeasyPyPS