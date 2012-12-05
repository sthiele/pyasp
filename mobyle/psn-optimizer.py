#!/usr/bin/python
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
def print_booleanmodel(f,model) :
    for i in model:
        if i.pred() == "i" :
           if i.arg(1) == "a" :f.write( '\n   '+str(i.arg(0))[1:-1]+' 1 '+str(i.arg(2))[1:-1])
           if i.arg(1) == "i" :f.write( '\n   '+str(i.arg(0))[1:-1]+' -1 '+str(i.arg(2))[1:-1])

def print_graph(f,model) :
    for a in model :
        if a.pred() == "edge" :
           f.write('\n'+str(a.arg(0))[1:-1] +' '+ str(a.arg(2)) + ' ' + str(a.arg(1))[1:-1])
           

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
  
  f= open("models.txt",'w')
  
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
f.write('Compressed network:')
print_graph(f,compnet)

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
f.write('\n\nThe optimal model score is '+ str(optimum[0]))


models = pno.get_minimal_models(instance, optimum[0])

#count=0
#for model in models:
   #m = model.to_list()
   #f.write( "\nModel "+str(count)+':')
   #print_booleanmodel(f,m)
   #count+=1

matrix = []
for model in models:
   matrixrow={}
   m = model.to_list()
   for i in m:
     if i.pred() == "i" :
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
       else : print "ERROR psn-optimizer line 122: unknown hyperedge",he
   matrix.append(matrixrow)   
   
   
f.write('\nModels with score '+str(optimum[0])+' in Matrix output:')
hyperedges.sort(key=len)
for he in hyperedges :
  f.write( "\n"+he+"\t"),
  for model in matrix :
    if he in model: f.write(str(model[he])+" "),
    else : f.write("0 "),   
   
   
   
   
f.close()
clean_up()


