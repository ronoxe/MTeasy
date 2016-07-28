# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:34:08 2015

@author: David
"""
import numpy as np
import MTeasyPy
import pylab
from matplotlib import pyplot as plt
nF=2

back = np.genfromtxt('../../runData/back.dat')
zz = np.genfromtxt('../../runData/zz.dat')
zzz = np.genfromtxt('../../runData/zzz.dat')
alpha = np.genfromtxt('../../runData/alpha.dat')
sig = np.genfromtxt('../../runData/sigma.dat')


fig1 = plt.figure(1)
plt.plot(back[:,0], back[:,1], label = '$\phi$',linewidth = 2)
plt.plot(back[:,0], back[:,2], label ='$\chi$',linewidth = 2)
pylab.ylabel('fields', fontsize=20)
pylab.xlabel('N', fontsize=15)
title('Backgroud field evolution',fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("fields.png")


#plt.savefig("back.png")
 
fig2 = plt.figure(2)
plt.semilogy(alpha[:,0], 5.8e-11* np.absolute(alpha[:,nF*2 + 1 + (1-1)* 2*nF]), label = '$ < \phi \phi> $',linewidth = 2)
plt.semilogy(alpha[:,0], 5.8e-11*np.absolute(alpha[:,nF*2 + 1 + (2-1)* 2*nF]), label ='$ <\phi \chi>  $',linewidth = 2)
plt.semilogy(alpha[:,0], 5.8e-11*np.absolute(alpha[:,nF*2 + 2 + (2-1)* 2*nF]), label ='$ <\chi \chi> $',linewidth = 2)
pylab.ylabel('correlations', fontsize=20)
pylab.xlabel('N', fontsize=15)
title('2pt correlations evolution',fontsize=15)
grid(True)
plt.legend(fontsize=15,loc=7)
plt.savefig("fields2pt.png")


fig3 = plt.figure(3)
plt.semilogy(alpha[:,0], (5.8e-11)**2*np.absolute(alpha[:,nF*2 + 6*2*nF*2*nF + 1 + (1-1)* 2*nF +(1-1)*2*nF*2*nF]), label = '$ < \phi \phi \phi> $',linewidth = 2)
plt.semilogy(alpha[:,0], (5.8e-11)**2*np.absolute(alpha[:,nF*2 +6*2*nF*2*nF + 1 + (1-1)* 2*nF +(2-1)*2*nF*2*nF]), label = '$ < \phi \phi \chi >$',linewidth = 2)
plt.semilogy(alpha[:,0], (5.8e-11)**2*np.absolute(alpha[:,nF*2 +6*2*nF*2*nF + 1 + (2-1)* 2*nF +(2-1)*2*nF*2*nF]), label = '$ < \phi \chi \chi > $',linewidth = 2)
plt.semilogy(alpha[:,0], (5.8e-11)**2*np.absolute(alpha[:,nF*2 +6*2*nF*2*nF + 2 + (2-1)* 2*nF +(2-1)*2*nF*2*nF]), label = '$ < \chi \chi \chi > $',linewidth = 2)


pylab.ylabel('correlations', fontsize=20)
pylab.xlabel('N', fontsize=15)
title('3pt correlations evolution',fontsize=15)
grid(True)
plt.legend(fontsize=15,loc=7)
plt.savefig("fields3pt.png")


fig4 = plt.figure(4)
plt.semilogy(zzz[:,0], (5.8e-11)*np.absolute(zzz[:,1]), label = '$  <\zeta \zeta >$',linewidth = 2)
plt.semilogy(zzz[:,0], (5.8e-11)**2*np.absolute(zzz[:,4]), label = '$ < \zeta \zeta \zeta>$',linewidth = 2)
pylab.ylabel('correlations', fontsize=20)
pylab.xlabel('N', fontsize=15)
title('2pt and 3pt correlations of $\zeta$',fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("zeta.png")



zzz65 = np.genfromtxt('../../savedResults/paper/DQ6d5/zzz.dat')
zzz6 = np.genfromtxt('../../savedResults/paper/DQ6/zzz.dat')
zzz55 = np.genfromtxt('../../savedResults/paper/DQ5d5/zzz.dat')
zzz5 = np.genfromtxt('../../savedResults/paper/DQ5/zzz.dat')
zzz4 = np.genfromtxt('../../savedResults/paper/DQ4/zzz.dat')
zzz3 = np.genfromtxt('../../savedResults/paper/DQ3/zzz.dat')
zzz35 = np.genfromtxt('../../savedResults/paper/DQ3d5/zzz.dat')
zzz45 = np.genfromtxt('../../savedResults/paper/DQ4d5/zzz.dat')
zzz425 = np.genfromtxt('../../savedResults/paper/DQ4d25/zzz.dat')
zzz475 = np.genfromtxt('../../savedResults/paper/DQ4d75/zzz.dat')
zzz575 = np.genfromtxt('../../savedResults/paper/DQ5d75/zzz.dat')
zzz325 = np.genfromtxt('../../savedResults/paper/DQ3d25/zzz.dat')
zzz625 = np.genfromtxt('../../savedResults/paper/DQ6d25/zzz.dat')
zzz375 = np.genfromtxt('../../savedResults/paper/DQ3d75/zzz.dat')
zzz525 = np.genfromtxt('../../savedResults/paper/DQ5d25/zzz.dat')
zzz675 = np.genfromtxt('../../savedResults/paper/DQ6d75/zzz.dat')



fig5 = plt.figure(5)
plt.semilogy(zzz65[:,0], np.absolute(zzz65[:,4]), label = '$  <\zeta \zeta \zeta>_{6.5}$',linewidth = 2)
plt.semilogy(zzz6[:,0], np.absolute(zzz6[:,4]), label = '$  <\zeta \zeta \zeta_2>_6$', linewidth = 2)
plt.semilogy(zzz5[:,0], np.absolute(zzz5[:,4]), label = '$  <\zeta \zeta \zeta_3>_5$',linewidth = 2)
plt.semilogy(zzz475[:,0], np.absolute(zzz475[:,4]), label = '$  <\zeta \zeta \zeta>_{4.25}$',linewidth = 2)
plt.semilogy(zzz45[:,0], np.absolute(zzz45[:,4]), label = '$  <\zeta \zeta \zeta>_{4.5}$',linewidth = 2)
plt.semilogy(zzz425[:,0], np.absolute(zzz425[:,4]), label = '$  <\zeta \zeta \zeta>_{4.25}$',linewidth = 2)
plt.semilogy(zzz4[:,0], np.absolute(zzz4[:,4]), label = '$  <\zeta \zeta \zeta>_4$',linewidth = 2)
plt.semilogy(zzz35[:,0], np.absolute(zzz35[:,4]), label = '$  <\zeta \zeta \zeta>_{3.5}$',linewidth = 2)
plt.semilogy(zzz3[:,0], np.absolute(zzz3[:,4]), label = '$  <\zeta \zeta \zeta>_3$',linewidth = 2)
plt.semilogy(zzz675[:,0], np.absolute(zzz675[:,4]), label = '$  <\zeta \zeta \zeta>_{6.75}$',linewidth = 2)




B=np.array([3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75])
zzzD = np.array([(5.8e-11)**2*zzz3[-1,4],(5.8e-11)**2*zzz325[-1,4],(5.8e-11)**2*zzz35[-1,4],(5.8e-11)**2*zzz375[-1,4],(5.8e-11)**2*zzz4[-1,4],(5.8e-11)**2*zzz425[-1,4],(5.8e-11)**2*zzz45[-1,4],(5.8e-11)**2*zzz475[-1,4],(5.8e-11)**2*zzz5[-1,4],(5.8e-11)**2*zzz525[-1,4],(5.8e-11)**2*zzz55[-1,4],(5.8e-11)**2*zzz575[-1,4],(5.8e-11)**2*zzz6[-1,4],(5.8e-11)**2*zzz625[-1,4],(5.8e-11)**2*zzz65[-1,4],(5.8e-11)**2*zzz675[-1,4]])

fig6 = plt.figure(6)
scatter(B,zzzD)
plt.ylim((4.6e-19,4.9e-19))
title('3pt depedance on e-folds pre horizon crossing',fontsize=15)
grid(True)
plt.legend(fontsize=15)
pylab.ylabel('$<\zeta \zeta \zeta>_{end}$', fontsize=20)
pylab.xlabel('N before horizon crossing', fontsize=15)
grid(True)
plt.legend(fontsize=15)
plt.savefig("stability.png")

zzz4 = np.genfromtxt('../../savedResults/paper/AX4/zzz.dat')
fig7 = plt.figure(7)
plt.plot(zzz4[:,0], 5.0/6.0*np.divide(zzz4[:,4], (np.multiply(zzz4[:,2],zzz4[:,3])+np.multiply(zzz4[:,1],zzz4[:,2]) + np.multiply(zzz4[:,1],zzz4[:,3]))), label = 'equilateral',linewidth = 2)
plt.plot(zzz[:,0], 5.0/6.0*np.divide(zzz[:,4], (np.multiply(zzz[:,2],zzz[:,3])+np.multiply(zzz[:,1],zzz[:,2]) + np.multiply(zzz[:,1],zzz[:,3]))), label = 'squeezed',linewidth = 2)

title('Evolution of the reduced bispectrum',fontsize=15)
grid(True)
plt.ylim((-70,30))
pylab.ylabel(r'$f_{\rm NL}$', fontsize=20)
pylab.xlabel('N ', fontsize=15)
grid(True)
plt.legend(fontsize=15,loc=4)
plt.savefig("fnlDQ.png")


plt.show(fig1)
plt.show(fig2)
plt.show(fig3)
plt.show(fig4)
plt.show(fig5)
plt.show(fig6)
plt.show(fig7)
