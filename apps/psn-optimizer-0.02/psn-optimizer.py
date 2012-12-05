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

from bioasp.query import psnoptimization as pno 
import sys
import os
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.data import psn


#some pretty printers for the boolean model
def print_booleanmodel(model) :
    for i in model:
        if i.pred() == "i" :
           if i.arg(1) == "a" :print '   '+str(i.arg(0))[1:-1]+' 1 '+str(i.arg(2))[1:-1]
           if i.arg(1) == "i" :print '   '+str(i.arg(0))[1:-1]+' -1 '+str(i.arg(2))[1:-1]

def print_graph(model) :
    for a in model :
        if a.pred() == "edge" :
           print str(a.arg(0))[1:-1] +' '+ str(a.arg(2)) + ' ' + str(a.arg(1))[1:-1]
           

def clean_up() :
  if os.path.isfile("parser.out"): os.remove("parser.out")
  if os.path.isfile("parsetab.py"): os.remove("parsetab.py")
  if os.path.isfile("parsetab.pyc"): os.remove("parsetab.pyc")
  if os.path.isfile("asp_py_lextab.py"): os.remove("asp_py_lextab.py")
  if os.path.isfile("asp_py_lextab.pyc"): os.remove("asp_py_lextab.pyc")
  if os.path.isfile("asp_py_parsetab.py"): os.remove("asp_py_parsetab.py")
  if os.path.isfile("asp_py_parsetab.pyc"): os.remove("asp_py_parsetab.pyc")
  if os.path.isfile("sif_parser_lextab.py"): os.remove("sif_parser_lextab.py")
  if os.path.isfile("sif_parser_lextab.pyc"): os.remove("sif_parser_lextab.pyc")
  
def m_quit() :  
  clean_up()
  quit()
         

if len(sys.argv)==3:
  if not os.path.isfile(sys.argv[1]):
    print "File not found:",sys.argv[1]
    m_quit()
  if not os.path.isfile(sys.argv[2]):
    print "File not found:",sys.argv[2]
    m_quit()
    
  net_string=sys.argv[1]
  obs_string=sys.argv[2]

else:
  print "\nmissing arguments !!!"
  print "usage: "
  print "   python psn-optimizer.py networkfile observationsfile"
  print "\nexample:"
  print "   python psn-optimizer.py data/LiverPKNDREAM.sif data/Dataset1.csv\n"
  exit(0)
  
  
print '\nReading network',net_string, '...',
net=psn.readSIF(net_string)
print 'done.'
#print net

print '\nReading observations',obs_string, '...',
obs=psn.readMIDAS(obs_string)
print 'done.'

print '\nCreate instance...\n'
instance = net.union(obs)
#print instance
print '\ncompress network ',obs_string, '...',
compnet = pno.compress(instance)
print 'done.'
print 'compressed network:'
print_graph(compnet)


instance= compnet.union(obs)

#get hypergraph
hyperedges = []
hg = pno.get_hypergraph(instance)
for a in hg:
   if a.pred() == "subset" : 
      #print a
      hn = a.arg(0)[1:-1]
      next = a.arg(2)
      while next != "end" :
	hn= hn+"+"+next.arg(0)[1:-1]
	next = next.arg(2)
      hn = hn+"="+a.arg(4)[1:-1]
      hyperedges.append(hn)
      


print '\nComputing score of the optimal BooleanModel...',
optimum = pno.get_minimal_model_size(instance)
print 'done.'    
print '   The optimal model score is', optimum[0][0],'.' 
#count clause atoms to determine minimal size
minimalsize= 0
for a in optimum[1][0] :
  if a.pred() == "clause":
    minimalsize = minimalsize+1
print '   The minimal size is', minimalsize,'.'    
    
suboptimal = int(optimum[0][0]*1.1)
print '\nComputing all models with score', suboptimal,'and size', minimalsize,'...',
models = pno.get_suboptimal_models(instance, suboptimal, minimalsize)
print 'done.'

matrix = []
for model in models:
   matrixrow={}
   m = model.to_list()
   for i in m:
     if i.pred() == "i" :
       #print i
       j=i.arg(0)
       #print j
       he =str(j.arg(0))[1:-1]
       next = j.arg(2)
       while next != "end" :
	 he= he+"+"+next.arg(0)[1:-1]
	 next = next.arg(2)
         #print he
       he =he+"="+str(i.arg(2))[1:-1]
       #print he      
       if he in hyperedges :
         if i.arg(1) == "a" : matrixrow[he]= 1
         if i.arg(1) == "i" : matrixrow[he]= -1
       else : print "ERROR psn-optimizer line 132: unknown hyperedge",he
   matrix.append(matrixrow)
   
print '\nPrint all models to models.txt.', 
f= open("models.txt",'w') 
f.write(str(len(models))+" Models found with score "+str(suboptimal)+" and size "+str(minimalsize)+"."),

hyperedges.sort(key=len)
for he in hyperedges :
  f.write ("\n"+str(he)+"\t"),
  for model in matrix :
    if he in model: f.write(str(model[he])),
    else : f.write("0"),
f.close()    
   
clean_up()


