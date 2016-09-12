#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from corelib import JAMLIB
from tools import save,load


x=1.0
save(x,'test.dat')
print load('test.dat')
sys.exit()






###################################
# pions
###################################
FFpion=JAMLIB('JAM16/FFpion')




# Quick test
x=0.05
Q2=1.0
print 'alphaS     = ',FFpion.get_alphaS(Q2)
print 'num pos    = ',FFpion.npos
print 'xF(ipos=0) = ',FFpion.get_XF(0,'up',x,Q2)   #<--- up == u + ub

# Plot FF
X=np.linspace(0.05,1)
XF=[[FFpion.get_XF(i,'up',x,Q2) for x in X] for i in range(FFpion.npos)]
XF0=np.mean(XF,axis=0)
dXF=np.std(XF,axis=0)

py.figure(figsize=(7,7))
ax=py.subplot(111)
ax.fill_between(X,XF0-dXF,XF0+dXF,edgecolor='r',facecolor='r',alpha=0.5)
ax.plot(X,XF0,'r-')
ax.set_ylim(0,1.5)
ax.set_xlim(np.amin(X),np.amax(X));
ax.set_xlabel(r'$z$',size=50)
ax.set_ylabel(r'$zD(z)$',size=50)
ax.tick_params(axis='both', which='major', labelsize=20)
#py.savefig('FFplot.pdf')    #<--- uncomment to save the figure in pdf format





