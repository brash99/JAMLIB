#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from corelib import JAMLIB
import lhapdf
from scipy.integrate import quad
from tools import load, save
###################################
def d2calc(interface='JAMLIB'):

  if interface=='JAMLIB':
    jamlib=JAMLIB('JAM15/T3PPDF')
    npos = jamlib.npos
    first = 0
  elif interface=='LHAPDF':
    jamlib=lhapdf.mkPDFs("JAM15_T3PPDF_Ceven")
    npos = 200
    first = 1
    
  def _get_d2(i,Q2,nucleon):
    if interface=='JAMLIB':   
      Du=lambda x: jamlib.get_XF(i,'up',x,Q2)/x
      Dd=lambda x: jamlib.get_XF(i,'dp',x,Q2)/x
    elif interface=='LHAPDF':
      Du=lambda x: jamlib[i].xfxQ2(902,x,Q2)/x
      Dd=lambda x: jamlib[i].xfxQ2(901,x,Q2)/x
    
    if nucleon=='proton':
      D=lambda x: 4/9.*Du(x) + 1/9.*Dd(x)
    if nucleon=='neutron':
      D=lambda x: 1/9.*Du(x) + 4/9.*Dd(x)
    return 2*quad(lambda x: x**2*D(x),0,1)[0]
  
  def get_d2(Q2,nucleon):
    d2=[_get_d2(i+first,Q2,nucleon) for i in  range(npos)]
    d2=[_d2 for _d2 in d2 if np.isnan(_d2)==False]
    return np.mean(d2),np.var(d2)**0.5
  
  Q2=np.arange(1,2,0.1)
  Q2=np.append(Q2,np.arange(2,3,0.2))
  Q2=np.append(Q2,np.arange(3,6.5,0.5))
  d2p=[get_d2(_Q2,'proton') for _Q2 in Q2]
  d2n=[get_d2(_Q2,'neutron') for _Q2 in Q2]
  
def d2plot():
  D=load('d2.dat')

  py.figure(figsize=(7,7))
  ax=py.subplot(111)
  Q2=D['Q2']
  d2p=np.array([x[0] for x in  D['d2p']])
  dd2p=np.array([x[1] for x in  D['d2p']])
  ax.fill_between(Q2,d2p-dd2p,d2p+dd2p\
    ,facecolor='r',edgecolor='r',alpha=0.2)
  ax.plot(Q2,d2p,'r-')

  d2n=np.array([x[0] for x in  D['d2n']])
  dd2n=np.array([x[1] for x in  D['d2n']])
  ax.fill_between(Q2,d2n-dd2n,d2n+dd2n\
    ,facecolor='b',edgecolor='b',alpha=0.2)
  ax.plot(Q2,d2n,'b-')


if __name__=='__main__':

  #int = 'JAMLIB'
  int = 'LHAPDF'
  d2calc(int)   # Note: can also be called w/o --> returns JAMLIB plot
  d2plot()
  py.savefig('d2_'+int+'.pdf')

