####################################### Axion quadratic simple example of basic MTeasy functions ###########################################
from matplotlib import pyplot as plt
import time
import imp  
#from pylab import *
import math 
import numpy as np
from scipy import interpolate 
import sympy as sym
import subprocess
############################################################################################################################################

#### Sets potential and compiles MTeasy, users may prefer to do this only once in a separate file (or comment after running below once) ####
############################### NB restart the 
nF=2
nP=3
f=sym.symarray('f',nF)
p=sym.symarray('p',nP)
V= 1./2. * p[0]**2 * f[0]**2 + p[2] * (1-sym.cos(2*math.pi * f[1] / p[1]))
import potentialSetup as p 
pvalue = np.zeros(3)
pvalue[0]=1.*pow(10.,-5); pvalue[1]=1.; pvalue[2]=25*pvalue[0]**2/4.0/math.pi**2;
p.potential(V,nF,nP) # writes this potential into c file when run 
subprocess.call(["python", "../moduleSetup.py", "install"]) # reintstalls MTeasy with new potential
############################################################################################################################################

import MTeasyPy; imp.reload(MTeasyPy)   # import recompilled module

########################### initial field values ###########################################################################################
fields = np.array([17.,.5-0.0001])
nF=MTeasyPy.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
V = MTeasyPy.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPy.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll
############################################################################################################################################


################################## run the background fiducial run #########################################################################
Nstart = 0.0
Nend = 50.0
t=np.linspace(Nstart, Nend, 1000)
back = MTeasyPy.backEvolve(t, initial,pvalue)
#back = np.genfromtxt('../../runData/back.dat', dtype=None) #load data into python
## plot background
#fig1 = plt.figure(1)
#plt.plot(back[:,0], back[:,2 ], 'r')
#plt.plot(back[:,0], back[:,1 ], 'r')
############################################################################################################################################


############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
Nexit = 13.0
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
k = np.exp(Nexit) *  MTeasyPy.H(backExit,pvalue); 
# other scales can then be defined wrt to k
############################################################################################################################################

start_time = time.time()

# example bispectrum run
# runs through many vaues of alpa and beta (throwing away those not within range, and calculates bispectrum at end of run)
# side is number of alphas and betas on each side it tries
side = 4
step =1.0/side

#np.zeros([side,side])
alphaA = np.array([])
betaA = np.array([]) 
nsnaps=5
snaps =  np.linspace(8.,Nend,nsnaps) 

biAOut = np.zeros(nsnaps) 
biROut = np.zeros(nsnaps)

print Nstart

for l in range(0,side+1):
      alpha =  -1.0 + l*2*step          
      for j in range(0,side+1):
         print "l " ,  l,   " j " , j  
         beta =j*step
         if alpha>-(1-beta) and alpha < 1-beta :            
             k1 = k/2 - beta*k/2.
             k2 = k/4*(1+alpha+beta)
             k3 = k/4*(1-alpha+beta)
             #print k1, k2, k3
# find conditions that allows longest mode 4 e-folds inside horizon
             Nstart = Nexit - max(math.log(k/k1),math.log(k/k2),math.log(k/k3))-4.0
             backExitMinus = np.zeros(2*nF)
             for i in range (1,2*nF+1):
                 backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

# run solver for this triangle
             t= np.concatenate((np.array([Nstart]),snaps))
             threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus, pvalue,0)
             zzz= threePt[:,:5]
             dim= np.abs(zzz[:,4]*(k1*k2*k3)**2)
             biA= np.array([])
             biR= np.array([])
             for ii in range(1,nsnaps+1):
                 bi=5.0/6.0*np.divide(zzz[ii,4], (np.multiply(zzz[ii,2],zzz[ii,3])+np.multiply(zzz[ii,1],zzz[ii,2]) + np.multiply(zzz[ii,1],zzz[ii,3])))
                 biR=np.append(biR,bi)
                 biA=np.append(biA,dim[ii])
#fig3 = plt.                 
             biAOut=np.vstack((biAOut,biA))
             biROut=np.vstack((biROut,biR))
             betaA=np.append(betaA,beta)
             alphaA=np.append(alphaA,alpha)
#plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
print  time.time() - start_time


for ii in range(0,nsnaps):
    print ii
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.interpolate import griddata
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    #ax = fig.gca(projection='3d')
    # note this: you can skip rows!
    X = alphaA
    Y = betaA
    Z = biROut[1:,ii]
    #Z = biAOut[1:,ii]
    ax.set_xticks([])                               
    ax.set_yticks([])                               
    ax.set_zticks([np.min(Z), np.max(Z)])
    ax.grid(False)   
    ax.set_ylabel(r'$\alpha$')
    ax.set_xlabel(r'$\beta$')
    ax.set_zlabel(r'$fNL$')    
    Nlabel =     'N = ' + format(t[ii+1], '.2f')
    #ax.title(Nlabel,fontsize=15)
    #plt.legend(fontsize=13)       
    xi = np.linspace(X.min(),X.max(),300)
    yi = np.linspace(Y.min(),Y.max(),300)
    # VERY IMPORTANT, to tell matplotlib how is your data organized
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')
    CS = plt.contourf(xi,yi,zi,50, extend3d=True)
    ax.clabel(CS, fontsize=9, inline=1)
    #ax.view_init(elev = 55)
#    ,vmin=-50, vmax=50)
    frame = 'movie/frame' + str(ii).zfill(3)+".png"
    plt.savefig(frame)    
    plt.close(fig)


#call biPlot.py which makes a better 3D plot of bispectrum as function of alpha beta(file is in useful directory, ensure that current directory is in path)
#import biPlot
