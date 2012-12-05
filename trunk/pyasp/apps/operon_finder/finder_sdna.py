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
from bioasp import asp

from bioasp.query import operon_finder as of 

import sys
import os
from bioasp.misc import *
import tempfile
from bioasp.data import op


  
      
def print_genes(answer, dictg) :
  for a in answer:
    if a.pred() == "cgene" :
      print dictg[int(a.arg(0))],
  print

  
def clean_up() :
  if os.path.isfile("parser.out"): os.remove("parser.out")
  if os.path.isfile("asp_py_lextab.py"): os.remove("asp_py_lextab.py")
  if os.path.isfile("asp_py_lextab.pyc"): os.remove("asp_py_lextab.pyc")
  if os.path.isfile("asp_py_parsetab.py"): os.remove("asp_py_parsetab.py")
  if os.path.isfile("asp_py_parsetab.pyc"): os.remove("asp_py_parsetab.pyc") 


if len(sys.argv)==4 :
  genome_string=sys.argv[1]
  metabolism_string=sys.argv[2]
  catalysation_string=sys.argv[3]
  
else:
  print "\nmissing arguments !!!"
  print "usage: "
  print "   python finder.py  genomefile metabolismfile catalysationfile"
  print "\nexample:"
  print "   python finder.py data/genome.txt data/metabolism.txt data/catalyze.txt\n"
  exit(0)
  
instance, dictg, dictr, couples = of.createInstance(genome_string,metabolism_string,catalysation_string) 
inst=instance.to_file()


for s,g in couples: 
  #print "compute best dna strand catalyzing pathway from reaction ",s," to ",g," ..."

  models = of.get_k_best_length(inst, s, g, 30, 5)
 
  if len(models) > 0 :
    print "\nk best genes catalyzing pathway from reaction",dictr[s],"to",dictr[g]
  c=1
  for m in models :
     print " ",str(c)+":",
     c+=1
     print_genes(m, dictg)
     
    
os.unlink(inst)
clean_up()
