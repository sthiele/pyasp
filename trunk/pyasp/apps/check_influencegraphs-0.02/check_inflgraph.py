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

from bioasp.query import influence_graphs as pm 
from bioasp.data import bioquali

import sys
import os
from bioasp.misc import *
import tempfile
from bioasp.data import op

#some pretty printers for predictions, minimal inconsistent cores, etc

def print_predictions(predictions) :
    #predictions.sort()
    predictions = sorted(predictions, key=lambda p: str(p.arg(0)))
    exp = ''
    for p in predictions:
        if p.pred() == "vlabel" :
           if exp!=str(p.arg(0)) :
              print 'Experiment '+str(p.arg(0))+':'
              exp=str(p.arg(0))
           if p.arg(2)== "-1": print('   '+str(p.arg(1))+ ' = - ')
           if p.arg(2)== "1" : print('   '+str(p.arg(1))+ ' = + ')
           if p.arg(2)== "0" : print('   '+str(p.arg(1))+ ' = nc ')
        if p.pred() == "elabel" :
           if p.arg(2) == "-1" : print '   '+str(p.arg(0))+' -> '+str(p.arg(1))+' -'
           if p.arg(2) == "1"  : print '   '+str(p.arg(0))+' -> '+str(p.arg(1))+' +'
                  
def print_mic(mic, net, obs):
  
    nodes = []
    edges = []
    for node in mic: nodes.append(str(node.arg(1)))
    
    predecessors = []
    for e in net:
       if e.pred() == "obs_elabel" :
          #print str(e)
          #print str(e.arg(0)),str(e.arg(1)),str(e.arg(2))
          if str(e.arg(1)) in nodes : 
            predecessors.append(str(e.arg(0)))
            if str(e.arg(2)) == "1" : edges.append( str(e.arg(0))+ " -> " + str(e.arg(1))+ " +")
            if str(e.arg(2)) == "-1" : edges.append(str(e.arg(0))+ " -> " + str(e.arg(1))+ " -")
         #TODO ? edges
    for edge in edges: print('   '+edge)
    for o in obs:
       if o.pred() == "obs_vlabel" :  
          if str(o.arg(1)) in nodes :
              if str(o.arg(2))=="1" :  print '   '+str(o.arg(1))+ " = +"
              if str(o.arg(2))=="-1" :  print '   '+str(o.arg(1))+ " = -"
          if str(o.arg(1)) in predecessors :
              if str(o.arg(2))=="1" :  print '   '+str(o.arg(1))+ " = +"
              if str(o.arg(2))=="-1" :  print '   '+str(o.arg(1))+ " = -"
    

def clean_up() :
  if os.path.isfile("parser.out"): os.remove("parser.out")
  if os.path.isfile("asp_py_lextab.py"): os.remove("asp_py_lextab.py")
  if os.path.isfile("asp_py_lextab.pyc"): os.remove("asp_py_lextab.pyc")
  if os.path.isfile("asp_py_parsetab.py"): os.remove("asp_py_parsetab.py")
  if os.path.isfile("asp_py_parsetab.pyc"): os.remove("asp_py_parsetab.pyc") 
  if os.path.isfile("graph_parser_lextab.py"): os.remove("graph_parser_lextab.py")
  if os.path.isfile("graph_parser_parsetab.py"): os.remove("graph_parser_parsetab.py")

def m_quit() :  
  clean_up()
  quit()

if len(sys.argv)==3 :
  
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
  print "   python check_inflgraph.py  networkfile observationfile"
  print "\nexample:"
  print "   python check_inflgraph.py data/yeast_guelzim.net data/yeast_snf2.obs\n"
  exit(0)
  
print '\nReading network',net_string, '...',
net = bioquali.readGraph(net_string)
print 'done.'

print '\nReading observations',obs_string, '...',
mu = bioquali.readProfile(obs_string)
print 'done.'

print '\nComputing input nodes ...',
inputs = pm.guess_inputs(net)
print 'done.'

print '\nTesting empty network for consistency ...',
consistent = pm.is_consistent(net)
print 'done.'

if consistent: print '   The empty network is consistent.'
else:
  print '   The empty network is inconsistent.'

  empty_net = net.union(inputs)
  print '\nTesting empty network plus input nodes for consistency ...',
  consistent = pm.is_consistent(empty_net)
  print 'done.'
  if consistent: print '   The empty network is consistent.'
  else:
    print '   The empty network is still inconsistent.'
    
    choice= raw_input('\nDo you want to compute the minimal inconsistent cores? y/n ')
    if choice=="y": 
       print '\nComputing minimal inconsistent cores (mic\'s) ...',
       mics = pm.get_minimal_inconsistent_cores(empty_net)
       print 'done.\n'
       count = 1
       oldmic = 0
       for mic in mics:
	  if oldmic != mic: 
	    print 'mic '+str(count)+':'
	    print_mic(mic.to_list(),net.to_list(),[])
	    count+=1
	    oldmic= mic

    print "\nWhich repair mode do you want to perform ?"
    print "[1] flip observed variations"
    print "[2] flip influences"
    print "[3] define network nodes as inputs"
    print "[4] define network nodes as input in an experiment (to use only in case of multiple experiments)"
    print "[5] add influences"
    repairchoice= raw_input('\nChoose repair mode 1-5:')
    repair_options= asp.TermSet()
    print '\nCompute repair options ...',
    if repairchoice=="1":                    
       print 'repair mode: flip observed variations ...',
       repair_options = pm.get_repair_options_flip_obs(empty_net)
       print 'done.'
    if repairchoice=="2":                    
       print 'repair mode: flip influences ...',
       repair_options = pm.get_repair_options_flip_edge(empty_net)
       print 'done.'
    if repairchoice=="3":                    
       print 'repair mode: define network nodes as inputs ...',
       repair_options = pm.get_repair_options_make_node_input(empty_net)
       print 'done.'
    if repairchoice=="4":
       print 'repair mode: define network nodes as input in an experiment ...',
       repair_options = pm.get_repair_options_make_obs_input(empty_net)
       print 'done.'                    
    if repairchoice=="5":
       print 'repair mode: add influence ...',
       repair_options = pm.get_repair_options_add_edges(empty_net)
       print 'done.'                    
   
    if not (repairchoice in ["1","2","3","4","5"]) : print "ERROR", exit(0)

    print '\nCompute minimal numbers of necessary repair operations ...',
    optimum = pm.get_minimum_of_repairs(empty_net,repair_options) 
    print 'done.'    
    
    print '   The data set can be repaired with minimal', optimum[0],'operations.'
    
    do_repair= raw_input('\nDo you want to compute all possible repair sets? Y/n:')
    if do_repair=="Y":

       print '\nComputing all repair sets with size', optimum[0],'...'
       models = pm.get_minimal_repair_sets(empty_net,repair_options,optimum[0])
       print 'done.'
    
       count = 1
       oldmodel = 0
       for model in models:
         if oldmodel != model: 
           oldmodel = model
           repairs = model.to_list()
           print '  repair',count,':'
           for r in repairs : print str(r.arg(0)),
           print ' '
           count+=1
    
    print '\nComputing predictions that hold under all repair sets size', optimum[0],'...',
    model = pm.get_predictions_under_minimal_repair(optimum[0],empty_net,repair_options)
    print 'done.'
    predictions = model.to_list()
    print ( str(len(predictions)) + ' predictions found:')
    print_predictions(predictions)

if consistent: #if network is consistent add data
  net_with_data = net.union(mu).union(inputs)
  print '\nTesting network with data for consistency ...',
  consistent = pm.is_consistent(net_with_data)
  print "done."
  if consistent: 
    print '   The network and data are consistent.'
    
    print '\nComputing predictions under consistency ...',
    model = pm.get_predictions_under_consistency(net_with_data)
    print 'done.'
    predictions = model.to_list()
    #predictions.sort()
    print ( str(len(predictions)) + ' predictions found:')
    print_predictions(predictions)
    
    
  else:
    print '   The network and the data are inconsistent.'
    
    choice= raw_input('\nDo you want to compute the minimal inconsistent cores? y/n ')
    if choice=="y": 
       print '\nComputing minimal inconsistent cores (mic\'s) ...',
       mics = pm.get_minimal_inconsistent_cores(net_with_data)
       print 'done.'
       count = 1
       oldmic = 0    
       for mic in mics:
         if oldmic != mic: 
           print 'mic '+str(count)+':'
           print_mic(mic.to_list(),net.to_list(),mu.to_list())
           count+=1
           oldmic= mic
        

    print "\nWhich repair mode do you want to perform ?"
    print "[1] flip observed variations"
    print "[2] flip influences"
    print "[3] define network nodes as inputs"
    print "[4] define network nodes as input in an experiment (to use only in case of multiple experiments)"
    print "[5] add influences"
    repairchoice= raw_input('\nChoose repair mode 1-5:')
    repair_options= asp.TermSet()
    print '\nCompute repair options ...',
    if repairchoice=="1":
       print 'repair mode: flip observed variations ...',
       repair_options = pm.get_repair_options_flip_obs(net_with_data)
       print 'done.'
    if repairchoice=="2":                    
       print 'repair mode: flip influences ...',    
       repair_options = pm.get_repair_options_flip_edge(net_with_data)
       print 'done.'
    if repairchoice=="3":                    
       print 'repair mode: define network nodes as inputs ...',
       repair_options = pm.get_repair_options_make_node_input(net_with_data)
       print 'done.'
    if repairchoice=="4":                    
       print 'repair mode: define network nodes as input in an experiment ...',
       repair_options = pm.get_repair_options_make_obs_input(net_with_data)
       print 'done.'
    if repairchoice=="5":                    
       print 'repair mode: add influence ...',
       repair_options = pm.get_repair_options_add_edges(net_with_data)
       print 'done.'
   
    if not (repairchoice in ["1","2","3","4","5"]) : print "ERROR", exit(0)

    print '\nCompute minimal numbers of necessary repair operations ...',
    optimum = pm.get_minimum_of_repairs(net_with_data, repair_options)
    print 'done.'
        
    print '   The data set can be repaired with minimal', optimum[0],'operations.'
       
    do_repair= raw_input('\nDo you want to compute all possible repair sets? Y/n:')
    if do_repair=="Y":

       print '\nComputing all repair sets with size', optimum[0],'...',
       models = pm.get_minimal_repair_sets(net_with_data, repair_options, optimum[0])
       print "done."
       count = 1
       oldmodel = 0
       for model in models:
         if oldmodel != model: 
           oldmodel = model
           repairs = model.to_list()
           print "  repair",count,':'
           for r in repairs : print str(r.arg(0)),
           print ' '
           count+=1

    #print 'Computing subset minimal repairs ...'
    #print '\n', pm.subset_minimal_repair_flip_obs(net_with_data)
    
    print '\nComputing predictions that hold under all repair sets size', optimum[0],'...',
    model = pm.get_predictions_under_minimal_repair(net_with_data, repair_options, optimum[0])
    print "done."
    predictions = model.to_list()
    #predictions.sort()
    print ( str(len(predictions)) + ' predictions found:')
    print_predictions(predictions)

clean_up()  
