#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from corelib import JAMLIB

jamlib=JAMLIB('FFpions')

x,Q2=0.5,10.0

print 'alphaS     = ',jamlib.get_alphaS(Q2)
print 'num pos    = ',jamlib.npos
print 'xF(ipos=0) = ',jamlib.get_XF(0,'FFpion','up',x,Q2)





