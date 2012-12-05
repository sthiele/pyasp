# Copyright (c) 2012, Sven Thiele <sthiele78@gmail.com>
#
# This file is part of BioASP.
#
# BioASP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioASP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioASP.  If not, see <http://www.gnu.org/licenses/>.import random

# -*- coding: utf-8 -*-
import re
from bioasp.asp import *



def Readcorrelationmatrix(filename) :
  f = open(filename, 'r')
  lpfacts = TermSet()
  fi = iter(f)
  i=0  
  for s in fi:
    if i==0 :
      #print s,
      lis = s.split('\n')
      genes = lis[0].split('\t')
    else:
      lis = s.split('\n')
      lis = lis[0].split('\t')
      #print lis
      gen = lis[0]
      j=1
      while j<=i:
         #print "norrelation("+gen+","+genes[j-1]+","+lis[j]+")."
        j+=1
      while (j> i) & (j<len(lis)):
        stri ="correlation("+gen+","+genes[j-1]+","+lis[j]+")."
        #print stri
        lpfacts.add(Term('correlation', [gen, genes[j-1],lis[j]]))
        #print "correlation("+gen+","+genes[j-1]+","+lis[j]+")."
        j+=1
    i+=1
    
  return lpfacts
   