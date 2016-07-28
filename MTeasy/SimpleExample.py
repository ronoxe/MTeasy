from matplotlib import pyplot as plt
import time
import math 
import numpy as np
import sys 
from scipy import interpolate 
############################################################################################################################################

#This file contains simple examples of using the MTeasy package for the heavy field example of Langlois.
#It assumes the LangHeavySeptup file has been run to install a LH version of MTeasyPy
#It is recommended you restart the kernel to insure any updates to MTeasyPyLH are imported 

location = "/Users/David/Dropbox/MTeasyDist/MTeasy/" # this should be the location of the MTeasy folder 
sys.path.append(location) # sets up python path to give access to MTeasySetup

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyLH;  # import module  


# Example 

########################### set initial field values and parameters for a simple example run ###################################################
shift=231.
fields = np.array([-2.-100.*math.sqrt(6.) +shift, 2.*math.tan(math.pi/20.)])

nF=MTeasyPyLH.nF() # gets number of fields (useful check)
nP=MTeasyPyLH.nP() # gets number of parameters needed (useful check)

params = np.zeros(nP)
params[0]=1.*pow(10.,-7); params[1]=1.*pow(10.,-4); params[2]=math.pi/10.; params[3] = -100.*math.sqrt(6.) +shift; params[4]=1000.*math.sqrt(3.);

V = MTeasyPyLH.V(fields,params) # calculate potential from some initial conditions
dV=MTeasyPyLH.dV(fields,params) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions using slow roll expression  
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 30
t=np.linspace(Nstart, Nend, 1000) # array at which output is returned
back = MTeasyPyLH.backEvolve(t, initial, params) # The output is read into the back numpy array 

# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################


############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
PhiExit = -100.*math.sqrt(6.) +shift +.9
import tests
k = tests.kexitPhi(PhiExit, 1, back, params, MTeasyPyLH) 
k= tests.kexitN(22.0, back, params, MTeasyPyLH) 
# other scales can then be defined wrt to k
############################################################################################################################################


################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds of massless sub-horizon evolution using a spline routine
NBMassless = 4.0
Nstart, backExitMinus = tests.ICsBMassless(NBMassless, k, back, params, MTeasyPyLH)  


t=np.linspace(Nstart,Nend, 1000)  # array at which output is returned -- itial valye should correspond to initial field values
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an crude estimate for ns

twoPt = MTeasyPyLH.sigEvolve(t, k, backExitMinus,params, 1) # puts information about the two point fuction in twoPt array
zz=twoPt[:,1] # the second column is the 2pt of zeta
zz1=zz[-1] 
twoPt=MTeasyPyLH.sigEvolve(t, k+.1*k, backExitMinus,params, 1) 
zz=twoPt[:,1]
zz2=zz[-1]
n_s = (np.log(zz2)-np.log(zz1))/(np.log(k+.1*k)-np.log(k))+4.0;
print n_s
fig2 = plt.figure(2)
plt.plot(t, zz,'r')
############################################################################################################################################


###################################### example bispectrum run ##############################################################################
# set three scales in FLS manner (using alpha beta notation)
alpha=0.
beta =1/3.

k1 = k/2 - beta*k/2.
k2 = k/4*(1+alpha+beta)
k3 = k/4*(1-alpha+beta)

# find initial conditions for 4 e-folds of massless evolution for the largest k (which stays inside the longest)

kM = np.max(np.array([k1,k2,k3]))
Nstart, backExitMinus = tests.ICsBMassless(NBMassless,kM, back, params, MTeasyPyLH)

print Nstart
print backExitMinus

start_time = time.time()
# run solver for this triangle
t=np.linspace(Nstart,Nend, 100000)
threePt = MTeasyPyLH.alphaEvolve(t,k1,k2,k3, backExitMinus,params,1) # all data from three point run goes into threePt array
ppp= threePt[:,4+2*nF+6*2*nF*2*nF+1:]        

zzz= threePt[:,:5] # this the 3pt of zeta

fig3 = plt.figure(3)
# plots the reduced bispectrum
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')

print  time.time() - start_time
############################################################################################################################################

plt.show(fig1)
plt.show(fig2)
plt.show(fig3)


