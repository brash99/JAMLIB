#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py

def display_logo(root):

  F=open('%s/DOC/logo'%root)
  L=F.readlines()
  L=[l.rstrip() for l in L]
  F.close()

  for l in L: print l


if __name__=='__main__':

  display_logo('../')

