#!/usr/bin/python
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

def print_predictions(f, predictions) :
    #predictions.sort()
    predictions = sorted(predictions, key=lambda p: str(p.arg(0)))
    exp = ''
    for p in predictions:
        if p.pred() == "vlabel" :
           if exp!=str(p.arg(0)) :
              f.write( 'Experiment '+str(p.arg(0))+':\n')
              exp=str(p.arg(0))
           if p.arg(2)== "-1": f.write('   '+str(p.arg(1))+ ' = - \n')
           if p.arg(2)== "1" : f.write('   '+str(p.arg(1))+ ' = + \n')
           if p.arg(2)== "0" : f.write('   '+str(p.arg(1))+ ' = nc \n')
        if p.pred() == "elabel" :
           if p.arg(2) == "-1" : f.write( '   '+str(p.arg(0))+' -> '+str(p.arg(1))+' -\n')
           if p.arg(2) == "1"  : f.write( '   '+str(p.arg(0))+' -> '+str(p.arg(1))+' +\n')
                  
def print_mic(f,mic, net, obs):
  
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
            if str(e.arg(2)) == "1" : edges.append( str(e.arg(0))+ " -> " + str(e.arg(1))+ " +\n")
            if str(e.arg(2)) == "-1" : edges.append(str(e.arg(0))+ " -> " + str(e.arg(1))+ " -\n")
         #TODO ? edges
    for edge in edges: f.write('   '+edge)
    for o in obs:
       if o.pred() == "obs_vlabel" :  
          if str(o.arg(1)) in nodes :
              if str(o.arg(2))=="1" :  f.write('   '+str(o.arg(1))+ " = +\n")
              if str(o.arg(2))=="-1" :  f.write('   '+str(o.arg(1))+ " = -\n")
          if str(o.arg(1)) in predecessors :
              if str(o.arg(2))=="1" :  f.write('   '+str(o.arg(1))+ " = +\n")
              if str(o.arg(2))=="-1" :  f.write('   '+str(o.arg(1))+ " = -\n")
    

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

if len(sys.argv)==5 :
  if not os.path.isfile(sys.argv[1]):
    print "File not found:",sys.argv[1]
    m_quit()
  if not os.path.isfile(sys.argv[2]):
    print "File not found:",sys.argv[2]
    m_quit()    
  
  net_string=sys.argv[1]
  obs_string=sys.argv[2]
  comp_mics=sys.argv[3]
  repairchoice = sys.argv[4]
  
  f= open("results.txt",'w')
    
else:
  print "\nmissing arguments !!!"
  print "usage: "
  print "   python check_inflgraph.py  networkfile observationfile computemics repairchoice"
  print "\nexample:"
  print "   python check_inflgraph.py data/yeast_guelzim.net data/yeast_snf2.obs 1 2\n"
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

if consistent: f.write('The network without variations is consistent.\n')
else:
  #f.write ('The network without variations is inconsistent.')

  empty_net = net.union(inputs)
  print '\nTesting empty network plus input nodes for consistency ...',
  consistent = pm.is_consistent(empty_net)
  print 'done.'
  if consistent: f.write('The network without variations is consistent.\n')
  else:
    f.write(' The network without variations is inconsistent.\n')
    

    if comp_mics=="1": 
       f.write('\nMinimal inconsistent cores (mic\'s):\n')
       mics = pm.get_minimal_inconsistent_cores(empty_net)

       count = 1
       oldmic = 0
       for mic in mics:
	  if oldmic != mic: 
	    f.write('\nmic '+str(count)+':\n')
	    print_mic(f,mic.to_list(),net.to_list(),[])
	    count+=1
	    oldmic= mic


    repair_options= asp.TermSet()
    f.write( '\nRepair:')
    if repairchoice=="1":                    
       f.write( 'repair mode: flip observed variations ')
       repair_options = pm.get_repair_options_flip_obs(empty_net)
    if repairchoice=="2":                    
       f.write( 'repair mode: flip influences ')
       repair_options = pm.get_repair_options_flip_edge(empty_net)
    if repairchoice=="3":                    
       f.write( 'repair mode: define network nodes as inputs ')
       repair_options = pm.get_repair_options_make_node_input(empty_net)
    if repairchoice=="4":
       f.write( 'repair mode: define network nodes as inputs in an experiment ')
       repair_options = pm.get_repair_options_make_obs_input(empty_net)
    if repairchoice=="5":
       f.write( 'repair mode: add influence ')
       repair_options = pm.get_repair_options_add_edges(empty_net)
   


    print '\nCompute minimal numbers of necessary repair operations ...',
    optimum = pm.get_minimum_of_repairs(empty_net,repair_options) 
    print 'done.'    
    
    f.write( '\nThe data set can be repaired with minimal ', optimum[0],' operations.\n')
    
    #do_repair= raw_input('\nDo you want to compute all possible repair sets? Y/n:')
    #if do_repair=="Y":

       #print '\nComputing all repair sets with size', optimum[0],'...'
       #models = pm.get_minimal_repair_sets(empty_net,repair_options,optimum[0])
       #print 'done.'
    
       #count = 1
       #oldmodel = 0
       #for model in models:
         #if oldmodel != model: 
           #oldmodel = model
           #repairs = model.to_list()
           #print '  repair',count,':'
           #for r in repairs : print str(r.arg(0)),
           #print ' '
           #count+=1
    
    f.write( '\nPredictions that hold under all repair sets size'+str(optimum[0])+"\n")
    model = pm.get_predictions_under_minimal_repair(optimum[0],empty_net,repair_options)
    predictions = model.to_list()
    f.write( ( str(len(predictions)) + ' predictions found\n'))
    print_predictions(f, predictions)

if consistent: #if network is consistent add data
  net_with_data = net.union(mu).union(inputs)
  print '\nTesting network with data for consistency ...',
  consistent = pm.is_consistent(net_with_data)
  print "done."
  if consistent: 
    f.write( '\nThe network and data are consistent.')
    
    f.write( '\nPredictions under consistency:\n')
    model = pm.get_predictions_under_consistency(net_with_data)
    predictions = model.to_list()
    #predictions.sort()
    f.write( str(len(predictions)) + ' predictions found\n')
    print_predictions(f, predictions)
    
    
  else:
    f.write( '\nThe network and the data are inconsistent.\n')
    
    if comp_mics=="1": 
       f.write( '\nMinimal inconsistent cores (mic\'s)')
       mics = pm.get_minimal_inconsistent_cores(net_with_data)
       print 'done.'
       count = 1
       oldmic = 0    
       for mic in mics:
         if oldmic != mic: 
           f.write( '\nmic '+str(count)+':\n')
           print_mic(f, mic.to_list(),net.to_list(),mu.to_list())
           count+=1
           oldmic= mic
        

    repair_options= asp.TermSet()
    f.write( '\nRepair ')
    if repairchoice=="1":
       f.write( 'repair mode: flip observed variations ')
       repair_options = pm.get_repair_options_flip_obs(net_with_data)
    if repairchoice=="2":                    
       f.write( 'repair mode: flip influences ')
       repair_options = pm.get_repair_options_flip_edge(net_with_data)
    if repairchoice=="3":                    
       f.write( 'repair mode: define network nodes as inputs ')
       repair_options = pm.get_repair_options_make_node_input(net_with_data)
    if repairchoice=="4":                    
       f.write( 'repair mode: define network nodes as inputs in an experiment ')
       repair_options = pm.get_repair_options_make_obs_input(net_with_data)
    if repairchoice=="5":                    
       f.write( 'repair mode: add influence ')
       repair_options = pm.get_repair_options_add_edges(net_with_data)


    print '\nCompute minimal numbers of necessary repair operations ...',
    optimum = pm.get_minimum_of_repairs(net_with_data, repair_options)
    print 'done.'
        
    f.write( '\nThe data set can be repaired with minimal '+str(optimum[0])+' operations.\n')
       
    #do_repair= raw_input('\nDo you want to compute all possible repair sets? Y/n:')
    #if do_repair=="Y":

       #print '\nComputing all repair sets with size', optimum[0],'...',
       #models = pm.get_minimal_repair_sets(net_with_data, repair_options, optimum[0])
       #print "done."
       #count = 1
       #oldmodel = 0
       #for model in models:
         #if oldmodel != model: 
           #oldmodel = model
           #repairs = model.to_list()
           #print "  repair",count,':'
           #for r in repairs : print str(r.arg(0)),
           #print ' '
           #count+=1

    #print 'Computing subset minimal repairs ...'
    #print '\n', pm.subset_minimal_repair_flip_obs(net_with_data)
    
    f.write( '\nPredictions that hold under all repair sets size ' + str(optimum[0]) + "\n")
    model = pm.get_predictions_under_minimal_repair(net_with_data, repair_options, optimum[0])
    predictions = model.to_list()
    #predictions.sort()
    f.write (str(len(predictions)) + ' predictions found\n')
    print_predictions(f,predictions)
f.close()
clean_up()  
