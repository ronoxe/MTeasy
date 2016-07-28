import numpy as np
from scipy import interpolate 
import time


def ICsBMassless(NBMassless,k, back, params, MTE):
    eigen = np.zeros(np.size(back[:,0]))
    nF = np.size(back[0,:])/2
    for jj in range(0, np.size(back[:,0])):
        w, v = np.linalg.eig(MTE.ddV(back[jj,1:1+nF],params))
        eigen[jj] = np.max(w)
    massEff = np.add(-k**2*np.exp(-2.0*back[:,0]),eigen)
    NMassless = interpolate.splev(0.0,interpolate.splrep(np.sort(massEff), back[:,0], s=1e-15),der=0)
    backExitMinus = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(NMassless-NBMassless,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        
    return NMassless-NBMassless, backExitMinus
    
def ICsBExit(NBexit, k, back, params, MTE):
    kvaH = np.zeros(np.size(back[:,0]))    
    nF = np.size(back[0,:])/2    
    for jj in range(0, np.size(back[:,0])):
        H=MTE.H(back[jj,1:1+2*nF],params)
        kvaH[jj] = -k+np.exp(back[jj,0])*H
    Nexit = interpolate.splev(0.0,interpolate.splrep(np.sort(kvaH), back[:,0], s=1e-15),der=0)
    
    backExitMinus = np.zeros(2*nF)
    for i in range (1,2*nF+1):
        backExitMinus[i-1] = interpolate.splev(Nexit-NBexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        
    return Nexit-NBexit, backExitMinus

def eqSpectraBM(kA, back, params, NBM, MTE):
    zzzOut=np.array([])
    zzOut=np.array([])
    times = np.array([])    
    for ii in range(0,np.size(kA)):
        k=kA[ii]        
        Nstart, backExitMinus = ICsBMassless(NBM, k, back, params, MTE)  
        t=np.linspace(Nstart,back[-1,0], 10)

        k1=k;k2=k;k3=k;
        # run solver for this triangle
        start_time = time.time()        
        threePt = MTE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,1) # all data from three point run goes into threePt array
        zzOut=np.append(zzOut, threePt[-1,1])
        zzzOut=np.append(zzzOut, threePt[-1,4]) 
        times = np.append(times, time.time()-start_time)
        
    return zzOut, zzzOut, times



def alpBetSpectraBM(kt,alpha, beta, back, params, NBM, nsnaps, MTE):

    biAOut = np.array([]) 
    k1 = np.array([])
    k2 = np.array([])
    k3 = np.array([])
  
    times = np.array([])
    
    Hin = np.array([])    

    Nend = back[0,-1]
    
    for jj in range(0,np.size(back[:,0])):
            Hin[jj]=MTE.H(back[jj,1:])
            
    Nexit = interpolate.splev(kt/3., np.exp(back[:,0])*Hin,interpolate.splrep(back[:,1], back[:,0], s=1e-15),der=0)

    if (nsnaps ==1 or nsnaps == 0):
        snaps= np.array([Nend]) 

    snaps =  np.linspace(Nexit-2.0,Nend,nsnaps) 


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

        
def eqSpectraBE(kA, back, params, NBM, MTE):
    zzzOut=np.array([])
    zzOut=np.array([])
    times = np.array([])    
    for ii in range(0,np.size(kA)):
        k=kA[ii]        
        Nstart, backExitMinus = ICsBexit(NBM, k, back, params, MTE)  
        t=np.linspace(Nstart,back[-1,0], 10)

        k1=k;k2=k;k3=k;
        # run solver for this triangle
        start_time = time.time()        
        threePt = MTE.alphaEvolve(t,k1,k2,k3, backExitMinus,params,1) # all data from three point run goes into threePt array
        zzOut=np.append(zzOut, threePt[-1,1])
        zzzOut=np.append(zzzOut, threePt[-1,4]) 
        times = np.append(times, time.time()-start_time)
        
    return zzOut, zzzOut, times
                
        
def kexitN(Nexit, back, params, MTE):
    nF = np.size(back[0,:])/2    
    backExit = np.zeros(2*nF)

    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        k = np.exp(Nexit) *  MTE.H(backExit,params);
    return k


def kexitPhi(PhiExit, n, back, params, MTE):
    nF = np.size(back[0,:])/2    
    backExit = np.zeros(2*nF)
    Nexit = interpolate.splev(PhiExit,interpolate.splrep(back[:,n], back[:,0], s=1e-15),der=0)

    for i in range (1,2*nF+1):
        backExit[i-1] = interpolate.splev(Nexit,interpolate.splrep(back[:,0], back[:,i], s=1e-15),der=0)
        k = np.exp(Nexit) *  MTE.H(backExit,params);
    return k
    

