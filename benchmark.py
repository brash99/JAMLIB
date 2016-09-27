#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from python.corelib import JAMLIB
import lhapdf
from fortran import flib

P1=JAMLIB('python/JAM16/FFpion')
P2 = lhapdf.mkPDF("JAM16_FF_pi_Ceven",1)
#P3 = flib.grid_init('fortran/JAM16','FFpion',0)
flib.grid_init('JAM16','FFpion',0)

x=0.01
Q2=10.0
print 
print 'x=',x
print 'Q2=',Q2
print 'python  = ',P1.get_XF(0,'up',x,Q2)
print 'lhadpf  = ',P2.xfxQ2(901,x,Q2)
print 'fortran = ',flib.get_xf(x,Q2,'up')


x=0.9
Q2=10.0
print 
print 'x=',x
print 'Q2=',Q2
print 'python  = ',P1.get_XF(0,'up',x,Q2)
print 'lhadpf  = ',P2.xfxQ2(901,x,Q2)
print 'fortran = ',flib.get_xf(x,Q2,'up')


x=1.0
Q2=10.0
print 
print 'x=',x
print 'Q2=',Q2
print 'python  = ',P1.get_XF(0,'up',x,Q2)
print 'lhadpf  = ',P2.xfxQ2(901,x,Q2)
print 'fortran = ',flib.get_xf(x,Q2,'up')




