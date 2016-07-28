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

def equilData(krange,back, Nbefore):
    zzzOutL=np.array([])
    zzOutL=np.array([])
    fnlOutL=np.array([])    
    times=np.array([])    
    Hin = np.array([])    
    
    for ii in range(0,np.size(krange)):
        k=krange[ ii]    
        k1 = k
        k2=k
        k3=k
        for jj in range(0,np.size(back[:,0])):
            Hin[jj]=MTeasyPy.H(back[jj,1:])
            
        Nexit = interpolate.splev(k, np.exp(back[:,0])*Hin,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)

        Nstart = Nexit - Nbefore
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
        times = np.append(times, time.time()-start_timeIn) 
    return (zzOutL, zzzOutL, times)



def alpBetData(k,back, Nstart, alphaA, betaA, nsnaps, Nbefore):
    step =1.0/side


    biAOut = np.array([]) 
    biRout = np.array([])
    times = np.array([])
    
    Hin = np.array([])    

    Nend = back[0,-1]
    step =1.0/sides

    snaps =  np.linspace(Nstart,Nend,nsnaps) 
    if nsnaps ==1:
        snaps= np.array([Nend]) 


    biAOut = np.zeros(nsnaps) 
    biROut = np.zeros(nsnaps)
    biROut = np.zeros(nsnaps)


    for jj in range(0,np.size(back[:,0])):
            Hin[jj]=MTeasyPy.H(back[jj,1:])
            
    Nexit = interpolate.splev(k, np.exp(back[:,0])*Hin,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)


    for l in range(0,np.size(alphaA)):
      alpha =  alphaA[l]          
      for j in range(0,side+1):
         print "l " ,  l,   " j " , j  
         beta =j*step
         if alpha>-(1-beta) and alpha < 1-beta :            
             k1 = k/2 - beta*k/2.
             k2 = k/4*(1+alpha+beta)
             k3 = k/4*(1-alpha+beta)
             #print k1, k2, k3
# find conditions that allows longest mode 4 e-folds inside horizon
             Nstart = Nexit - max(math.log(k/k1),math.log(k/k2),math.log(k/k3))-Nbefore
             backExitMinus = np.zeros(2*nF)
             for i in range (1,2*nF+1):
                 backExitMinus[i-1] = interpolate.splev(Nstart,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)

# run solver for this triangle
             t= np.concatenate((np.array([Nstart]),snaps))
             timebefore =     time.tim()            
             threePt = MTeasyPy.alphaEvolve(t,k1,k2,k3, backExitMinus,0)
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
             times = np.append(times,timebefore - time.time())
    return (alphaA, betaA, biAOut, biROut, times)
#plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))),'r')
