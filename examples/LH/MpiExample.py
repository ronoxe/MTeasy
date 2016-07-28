from matplotlib import pyplot as plt
import time
import imp  
from pylab import *
import math 
import numpy as np
from scipy import interpolate 
import sympy as sym
import subprocess
import MTeasyPy
from mpi4py import MPI

comm = MPI.COMM_WORLD

print "rank", comm.Get_rank()

########################### initial field values ###########################################################################################
nP=5
shift=231.
fields = np.array([-2.-100.*math.sqrt(6.) +shift, 2.*math.tan(math.pi/20.)])
pvalue = np.zeros(nP)
pvalue[0]=1.*pow(10.,-7); pvalue[1]=1.*pow(10.,-4); pvalue[2]=math.pi/10.; pvalue[3] = -100.*math.sqrt(6.) +shift; pvalue[4]=1000.*math.sqrt(3.);
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPy.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,np.array([0.,0.]))) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 33
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPy.backEvolve(t, initial, pvalue)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
# plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,2 ], 'r')
plt.plot(back[:,0], back[:,1 ], 'g')
############################################################################################################################################


points = 10
kOut=np.array([])
phiOut = np.array([])


for ii in range(0,points):
    PhiExit = P-100.*math.sqrt(6.) +shift + ii*0.1
    backExit = np.zeros(2*nF)
    Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)

    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        k = np.exp(Nexit) *  MTeasyPy.H(backExit,pvalue); 
    kOut= np.append(kOut, k)
    phiOut = np.append(phiOut, PhiExit)

rank=comm.Get_rank()
size=comm.Get_size()
num = points/size;

kOutL=np.array([])
zzzOutL=np.array([])
zzOutL=np.array([])
fnlOutL=np.array([])
kOutL=np.array([])
PhiOutL=np.array([])

start_time = time.time()
for ii in range(0,num):
    PhiExit = phiOut[rank*num + ii]
    print "element", rank*num + ii
    backExit = np.zeros(2*nF)
    Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)
    k=kOut[rank*num + ii]    
    k1 = k
    k2= k
    k3= k

    Nstart = Nexit - 10.0
    backExitMinus = np.zeros(2*nF)

    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

    Nstart = Nexit - max(log(k/k1),log(k/k2),log(k/k3))-3.0
    backExitMinus = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

    start_timeIn = time.time()

# run solver for this triangle
 
    t=np.linspace(Nstart,Nend, 10)
    threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus, pvalue, 1)
    zzz= threePt[:,:5]
    zzzOutL = np.append(zzzOutL, zzz[-1,4])
    zzOutL = np.append(zzOutL, zzz[-1,1])
    fnlOutL = np.append(fnlOutL, 5.0/6.0*np.divide(zzz[-1,4], (np.multiply(zzz[-1,2],zzz[-1,3])+np.multiply(zzz[-1,1],zzz[-1,2]) + np.multiply(zzz[-1,1],zzz[-1,3]))))
    PhiOutL = np.append(PhiOutL, PhiExit)
    print "element", rank*num+ii, "time taken", time.time() - start_timeIn

comm.Send(fnlOutL,dest=0)
comm.Send(PhiOutL,dest=0)
print "rank", rank, "time taken", time.time() - start_timeIn


if rank == 0:
    fnltot = np.array([])    
    phitot = np.array([])    
    for jj in range(0,size):    
        comm.Recv(fnlOutL,source = jj)
        comm.Recv(PhiOutL,source = jj)
        fnltot=np.append(fnltot,fnlOutL)
        phitot=np.append(phitot,PhiOutL)
    print fnltot
    print phitot   
    print phiOut
    print "total time taken",  time.time()-start_time