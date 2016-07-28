from matplotlib import pyplot as plt
import time
import math 
import numpy as np
import sys 
from scipy import interpolate as interp

location = "/home/jwr/Code/MTeasy"
sys.path.append(location)

import MTeasySetup
MTeasySetup.pathSet()  # his add sets the other paths that MTeasy uses

import MTeasyPyPS;  # import module  
pvalue = np.zeros(3)
M=1.*pow(10.,-5)
pvalue[0]=1./2.*30.0 * M**2; pvalue[1]=1./2.*300.0 * M**2; pvalue[2]=1./2. * 30./81. *M**2;

########################## initial field values ###########################################################################################

fields = np.array([10.0,0.01,13.0])
nF=MTeasyPyPS.nF() # use number of fields routine to get number of fields (can be used as a useful cross check)
nP=MTeasyPyPS.nP() # gets number of parameters needed (useful check)

V = MTeasyPyPS.V(fields,pvalue) # calculate potential from some initial conditions
dV=MTeasyPyPS.dV(fields,pvalue) # calculate derivatives of potential (changes dV to derivatives)
initial = np.concatenate((fields,dV/np.sqrt(3.*V))) # set initial conditions to be in slow roll


################################ run the background fiducial run #########################################################################
Nstart = 0
Nend = 50.0
t=np.linspace(Nstart, Nend, 500)
back = MTeasyPyPS.backEvolve(t, initial,pvalue)

####plot background
fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,3 ], 'b',label=r'$\phi_3$')
plt.plot(back[:,0], back[:,2 ], 'r',label=r'$\phi_2$')
plt.plot(back[:,0], back[:,1 ], 'g',label=r'$\phi_1$')
plt.title(r'Evolution of the Fields',fontsize=15)
plt.grid(True)
plt.legend(fontsize=15)
plt.ylabel(r'Field Value', fontsize=20)
plt.xlabel(r'e-folds', fontsize=15)
plt.grid(True)
plt.legend(fontsize=15)
plt.savefig("FieldsD3.png")
plt.show()
###############################################################################
############################################################################################################################################
# set a pivot scale which exits after certain time using the background run -- a spline
# is used to find field and field velocity values after Nexit number of e-folds, this gives H, and 
# then k=aH gives the k pivot scale
# in this example we treat this scale as k_t
Nexit = 15.7
backExit = np.zeros(2*nF)

for i in range (1,2*nF+1):
    backExit[i-1] = interp.splev(Nexit,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)
ke = np.exp(Nexit) *  MTeasyPyPS.H(backExit,pvalue);
print ke, 'end'


Nexit = 13
backExit = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExit[i-1] = interp.splev(Nexit,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)
k = np.exp(Nexit) *  MTeasyPyPS.H(backExit,pvalue);
print k, 'start'
################################# example 2pt run ##########################################################################################
#find conditions for 4 e-folds sub-horizon evolution using another spline
Nstart = Nexit - 4.0

backExitMinus = np.zeros(2*nF)
for i in range (1,2*nF+1):
    backExitMinus[i-1] = interp.splev(Nstart,interp.splrep(back[:,0], back[:,i], s=1e-15),der=0)

t=np.linspace(Nstart,Nend, 500)
# run the sigma routine to calc and plot the evolution of power spectrum value for this k -- can be 
# repeated to build up the spectrum, here we run twice to get an estimate for ns
print 'quack'

twoPt = MTeasyPyPS.sigEvolve(t, k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz1=zz[-1]
#twoPt=MTeasyPyPS.sigEvolve(t, k+.1*k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
zz2=zz[-1]
fig2 = plt.figure(2)
plt.plot(t, zz,'r')
plt.title(r'Evolution of the two-point function at the pivot scale',fontsize=15)
plt.grid(True)
plt.ylabel(r'$\left<\zeta\zeta\right>$', fontsize=20)
plt.xlabel(r'e-folds', fontsize=15)
plt.grid(True)
plt.savefig("twopointD3.png")
plt.show()
twoPt=MTeasyPyPS.sigEvolve(t, ke, backExitMinus,pvalue, 1)
print 'QQQQQQQQQQQQ'

#n_s = (math.log(zz2)-math.log(zz1))/(math.log(k+.1*k)-math.log(k))+4;
#print 'n_s = ', n_s
#twoPt=MTeasyPyPS.sigEvolve(t,k, backExitMinus,pvalue, 1)
twoPt=MTeasyPyPS.sigEvolve(t,1.0*k, backExitMinus,pvalue, 1)
#twoPt=MTeasyPyPS.sigEvolve(t,k, backExitMinus,pvalue, 1)
zz=twoPt[:,1]
pp= [1.0]
kmode=50
kkk=np.linspace(k,ke,kmode)

#print np.log(ke/k), np.log(abs(1.0)), np.log(abs(zz2/zz1))
#plt.axis([0, np.log(ke/k), np.log(abs(1.0)), np.log(abs(zz2/zz1))])
#plt.ion()

#for kk in range(kmode-1):
#	tp = MTeasyPyPS.sigEvolve(t, kkk[kk], backExitMinus,pvalue, 1)
#	zz3=tp[:,1]
#	pp = np.append(pp,[zz3[-1]/zz1])
	#print np.log(abs(pp[kk]))
	# plt.scatter(np.log(kkk[kk]/k),np.log(abs(pp[kk]/zz1)))
	# plt.show()
	# plt.pause(0.05)
#	print 'k mode # ',kk, ' '
fig2 = plt.figure(2)
plt.plot(t, zz,'r')
plt.title(r'Evolution of the two-point function at the pivot scale',fontsize=15)
plt.grid(True)
plt.ylabel(r'$\left<\zeta\zeta\right>$', fontsize=20)
plt.xlabel(r'e-folds', fontsize=15)
plt.grid(True)
plt.savefig("twopointD3.png")
plt.show()
print pp
#plt.loglog(t, zz, basex=2)
#plt.grid(True)
fig3= plt.figure(3)
kkk=kkk/k
kp=np.log(kkk)
pp=np.log(abs(pp))
np.savetxt('PSout.out',(kp,pp))
plt.plot(kp[1:],pp[1:])
plt.show()