####################################### Axion quadratic simple example of basic MTeasy functions ###########################################
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
#pset.potential(V,nF,nP,pvalue) # writes this potential into c file when run 

subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy with this new potential

############################################################################################################################################
import MTeasyPy; imp.reload(MTeasyPy)   # import recompilled module

########################### initial field values ###########################################################################################
fields = np.array([-2.-100.*math.sqrt(6.) +shift, 2.*math.tan(math.pi/20.)])
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPy.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 33
t=np.linspace(Nstart, Nend, 10)
back = MTeasyPy.backEvolve(t, initial, pvalue)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################


zzzOut=np.array([])
zzOut=np.array([])
fnlOut=np.array([])
kOut=np.array([])
PhiOut=np.array([])


for ii in range(0,10):
    PhiExit = -100.*math.sqrt(6.) +shift + ii*0.1
    print ii
    backExit = np.zeros(2*nF)
    Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)

    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        k = np.exp(Nexit) *  MTeasyPy.H(backExit,pvalue); 
# other scales can then be defined wrt to k

#example 2pt run 

#find conditions for 4 e-folds sub-horizon evolution using another spline
    Nstart = Nexit - 10.0
    backExitMinus = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

    t=np.linspace(Nstart,Nend, 10)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum
#MTeasyPy.sigEvolve(t, k, backE xitMinus, 1)
#zz = np.genfromtxt('../../runData/zz.dat', dtype=None)

# example bispectrum run
# set three scales in FLS manner
#alpha = 1./6 
    alpha=0.
#beta = 1./2
    beta = 1/3.
#beta = 0.99

    k1 = k/2 - beta*k/2.
    k2 = k/4*(1+alpha+beta)
    k3 = k/4*(1-alpha+beta)
    k1=k;k2=k;k3=k;
# find initial conditions that allows longest for this model we find 4 e-folds before massless condition
    Nstart = Nexit - max(log(k/k1),log(k/k2),log(k/k3))-3.0
    backExitMinus = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

    start_time = time.time()

# run solver for this triangle
    t=np.linspace(Nstart,Nend, 10)
    threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus, pvalue, 1)
    zzz= threePt[:,:5]
    zzzOut = np.append(zzzOut, zzz[-1,4])
    zzOut = np.append(zzOut, zzz[-1,1])
    kOut= np.append(kOut, k)
    fnlOut = np.append(fnlOut, 5.0/6.0*np.divide(zzz[-1,4], (np.multiply(zzz[-1,2],zzz[-1,3])+np.multiply(zzz[-1,1],zzz[-1,2]) + np.multiply(zzz[-1,1],zzz[-1,3]))))
    PhiOut = np.append(PhiOut, PhiExit)

    print  time.time() - start_time

    
fig2 = plt.figure(2)
plt.plot(PhiOut, fnlOut, 'g',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel('$fNL$', fontsize=20)
xlabel('$\phi^*$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.xlim(min(PhiOut),max(PhiOut))
plt.savefig("LH1.png")


fig3 = plt.figure(3)
plt.plot(log(kOut/kOut[-1]), fnlOut, 'g',linewidth = 2)
title('Reduced bispectrum in equilateral configuration',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$fNL$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH2.png")




fig4 = plt.figure(4)
G = 9./10*5./6.*1/3.*zzzOut*kOut**6/zzOut[-1]**2/kOut[-1]**6
plt.plot(log(kOut/kOut[-1]), G, 'g',linewidth = 2)
title(r'G quantity used in previous work',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'$G/k^3$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH3.png")


fig5 = plt.figure(5)
plt.plot(log(kOut/kOut[-1]), log(zzOut/zzOut[-1]*kOut**3/kOut[-1]**3), 'g',linewidth = 2)
title('Power spectrum',fontsize=15)
grid(True)
plt.legend(fontsize=15)
ylabel(r'${\cal P}/{\cal P}_{\rm pivot}$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15)
grid(True)
plt.legend(fontsize=15)
#plt.xlim((0,5.65))
plt.savefig("LH4.png")


plt.show(fig1)
plt.show(fig2)
plt.show(fig3)
plt.show(fig4)
plt.show(fig5)