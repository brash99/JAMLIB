#!/usr/bin/env python
import sys,os
import numpy as np
import time
import fnmatch
import cPickle 
from operator import mul

def checkdir(path):
  if not os.path.exists(path): 
    os.makedirs(path)
    return False
  else:
    return True

def save(data,name):  
  f=open(name,"w")
  cPickle.dump(data, f)
  f.close()

def load(name):  
  f=open(name,"r")
  data=cPickle.load(f)
  f.close()
  return data

class BAR(object):

  def __init__(self,msg,size):
    self.msg=msg
    self.size=size
    self.cnt=0

  def next(self):
    sys.stdout.write('\r')
    percentage=int(self.cnt/float(self.size)*100)
    sys.stdout.write('%s [%d%%]' % (self.msg,percentage))
    sys.stdout.flush()
    self.cnt+=1

  def finish(self):
    print 

def logo():
  L=[]
  L.append('################################################ ')
  L.append('                                                 ')
  L.append('     _   _    __  __ _     ___ ____              ')
  L.append('    | | / \  |  \/  | |   |_ _| __ )             ')
  L.append(' _  | |/ _ \ | |\/| | |    | ||  _ \             ')
  L.append('| |_| / ___ \| |  | | |___ | || |_) |            ')
  L.append(' \___/_/   \_\_|  |_|_____|___|____/             ')
  L.append('                                                 ')
  L.append('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ')
  L.append('                                                 ')
  L.append('Authors:                                         ')
  L.append('Nobuo Sato   nsato@jlab.org                      ')
  L.append('Jake Ethier                                      ')
  L.append('Wally Melnitchouk                                ')
  L.append('Alberto Accardi                                  ')
  L.append('################################################ ')
  for l in L: print l


  
