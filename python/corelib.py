import sys,os
import numpy as np
import cPickle
from scipy.interpolate import RectBivariateSpline
from tools import load,BAR,logo
 
class JAMLIB(object):
  
  def __init__(self,path):
    logo()
    self.setup_alphaS(order='NLO')
    self.load_tables(path)
    self.dist=path.split('/')[1]

  # alphaS

  def setup_alphaS(self,order='NLO'):

    if   order=='LO':  self.order=0
    elif order=='NLO': self.order=1

    self.aZ  = 0.118 /(4*np.pi)
    self.mc2 = 1.43**2        
    self.mb2 = 4.3**2
    self.mZ2 = 91.187**2
    self.mt2 = 172.9**2

    self.beta=np.zeros((7,3))
    for Nf in range(3,7): 
      self.beta[Nf,0]=11.0-2.0/3.0*Nf 
      self.beta[Nf,1]=102.-38.0/3.0*Nf 
      self.beta[Nf,2]=2857.0/2.0-5033.0/18.0*Nf+325.0/54.0*Nf**2 

    self.Q20=1.0
    # uses alphaS(mZ)--> backwards evolution
    self.ab=self.evolve_a(self.mZ2,self.aZ,self.mb2,5)
    self.at=self.evolve_a(self.mZ2,self.aZ,self.mt2,5)
    self.ac=self.evolve_a(self.mb2,self.ab,self.mc2,4)
    self.a0=self.evolve_a(self.mc2,self.ac,self.Q20,3)

    # we will store all Q2 values of alphaS 
    self.storage={}

  def get_Nf(self,Q2):
    Nf=3
    if Q2>=(4.*self.mc2): Nf+=1
    if Q2>=(4.*self.mb2): Nf+=1
    if Q2>=(4.*self.mt2): Nf+=1
    return Nf

  def beta_func(self,a,Nf):
    betaf = -self.beta[Nf,0]
    if self.order>=1: betaf+=-a*self.beta[Nf,1]
    if self.order>=2: betaf+=-a*self.beta[Nf,2]
    return betaf*a**2

  def evolve_a(self,Q20,a,Q2,Nf):
    # Runge-Kutta implemented in pegasus  
    LR = np.log(Q2/Q20)/20.0
    for k in range(20):
      XK0 = LR * self.beta_func(a,Nf)
      XK1 = LR * self.beta_func(a + 0.5 * XK0,Nf)
      XK2 = LR * self.beta_func(a + 0.5 * XK1,Nf)
      XK3 = LR * self.beta_func(a + XK2,Nf)
      a+= (XK0 + 2.* XK1 + 2.* XK2 + XK3) * 0.166666666666666
    return a

  def get_a(self,Q2):

    if Q2 in self.storage:
      return self.storage[Q2]
    else:
      if self.mt2<=Q2:
        Q20,a0,Nf=self.mt2,self.at,6
      elif self.mb2<=Q2 and Q2<self.mt2: 
        Q20,a0,Nf=self.mb2,self.ab,5
      elif self.mc2<=Q2 and Q2<self.mb2: 
        Q20,a0,Nf=self.mc2,self.ac,4
      elif Q2<self.mc2:
        Q20,a0,Nf=self.mc2,self.ac,3
      a=self.evolve_a(Q20,a0,Q2,Nf)
      self.storage[Q2]=a
      return a

  def get_alphaS(self,Q2):
    return self.get_a(Q2)*4*np.pi

  # interpolation for xF

  def load_tables(self,path):
    SPL=[]
    TAB=[]
    F=os.listdir(path)
    self.npos=len(F)
    bar=BAR('loading tables',len(F))
    for f in F:
      tab=load('%s/%s'%(path,f)) 
      X=tab['X']
      Q2=tab['Q2']
      self.X=X
      self.Q2=Q2
      spl={}
      for k in tab:
        if k=='X' or k=='Q2': continue
        spl[k]={}
        for kk in tab[k]:
          spl[k][kk]=RectBivariateSpline(X,Q2,tab[k][kk])
          tab[k][kk]=np.transpose(tab[k][kk])
      SPL.append(spl)
      TAB.append(tab)
      bar.next()
    bar.finish()
    self.SPL=SPL
    self.TAB=TAB
    self.Xgrid = X
    self.Q2grid = Q2
  
  def get_XF(self,ipos,flav,x,Q2):
    return self.SPL[ipos][self.dist][flav](x,Q2)[0,0]

  def get_XF_TAB(self,ipos,flav,iQ2):
    return self.X,self.TAB[ipos][self.dist][flav][iQ2,:]

  def get_Xgrid(self):
    return self.Xgrid
    
  def get_Q2grid(self):
    return self.Q2grid






