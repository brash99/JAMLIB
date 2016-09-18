#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from corelib import JAMLIB

###################################
# select distribution
###################################
#jamlib=JAMLIB('JAM15/PPDF')
#jamlib=JAMLIB('JAM15/T3PPDF')
#jamlib=JAMLIB('JAM15/g1T4')
jamlib=JAMLIB('JAM16/FFpion')
#jamlib=JAMLIB('JAM16/FFkaon')

###################################
# Quick test
###################################
x=0.05
Q2=1.0
flav='up'
print 'alphaS     = ',jamlib.get_alphaS(Q2)
print 'num pos    = ',jamlib.npos
print 'xF(ipos=0) = ',jamlib.get_XF(0,flav,x,Q2)


###################################
# Plot xF
###################################
X=np.linspace(0.01,1)
XF=[[jamlib.get_XF(i,flav,x,Q2) for x in X] for i in range(jamlib.npos)]
XF0=np.mean(XF,axis=0)
dXF=np.std(XF,axis=0)

py.figure(figsize=(7,7))
ax=py.subplot(111)
ax.fill_between(X,XF0-dXF,XF0+dXF,edgecolor='r',facecolor='r',alpha=0.5)
ax.plot(X,XF0,'r-')
ax.set_xlim(np.amin(X),np.amax(X));
ax.set_xlabel(r'$x$',size=50)
ax.set_ylabel(r'$xF(x)$',size=50)
ax.tick_params(axis='both', which='major', labelsize=20)
py.tight_layout()
py.savefig('xF.pdf')    





