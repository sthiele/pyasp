#!/usr/bin/python
# -*- coding: utf-8 -*-
from bioasp.query import network_expansion as pm 
from bioasp.data import sbml
import sys
import os
from bioasp import asp
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.asp import *



#some pretty printers for unproducible target, etc

def print_met(f,predictions) :
    for p in predictions: f.write('\n  '+str(p.arg(0)))
    
def clean_up() :
  if os.path.isfile("parser.out"): os.remove("parser.out")
  if os.path.isfile("parsetab.py"): os.remove("parsetab.py")
  if os.path.isfile("asp_py_lextab.py"): os.remove("asp_py_lextab.py")
  if os.path.isfile("asp_py_lextab.pyc"): os.remove("asp_py_lextab.pyc")
  if os.path.isfile("asp_py_parsetab.py"): os.remove("asp_py_parsetab.py")
  if os.path.isfile("asp_py_parsetab.pyc"): os.remove("asp_py_parsetab.pyc")

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
  if not os.path.isfile(sys.argv[3]):
    print "File not found:",sys.argv[3]
    m_quit()
  if not os.path.isfile(sys.argv[4]):
    print "File not found:",sys.argv[4]
    m_quit()   
  
  
  draftnetfile = sys.argv[1]
  repairnetfile = sys.argv[2]
  seedsfile = sys.argv[3]
  targetsfile =  sys.argv[4]
  
  f= open("analysis.txt",'w')
  
else:
  print "\nmissing arguments !!!"
  print "usage: "
  print "   python network_expansion.py draftnetwork repairnetwork seedfile targetfile"
  print "\nexample:"
  print "   python network_expansion.py data/Esicyc.sbml data/metaCyc.sbml  data/seeds.sbml  data/targets.sbml\n"
  exit(0)

print 'Reading draft network from ',draftnetfile,'...',
draftnet = sbml.readSBMLnetwork(draftnetfile, 'draft')
print 'done.'

print 'Reading targets from ',targetsfile,'...',
targets = sbml.readSBMLtargets(targetsfile)
print 'done.'

print 'Reading seeds from ',seedsfile,'...',
seeds = sbml.readSBMLseeds(seedsfile)
print 'done.'

print '\nChecking draftnet for unproducible targets ...',
model = pm.get_unproducible(draftnet, targets, seeds)
print 'done.'
f.write(str(len(model))+' unproducible targets:')
print_met(f,model.to_list())
unproducible_targets = TermSet()
for a in model :
   target= str(a)[13:]
   t = asp.String2TermSet(target)
   unproducible_targets = unproducible_targets.union(t)
   
print '\nReading repair network from ',repairnetfile,'...',
repairnet = sbml.readSBMLnetwork_with_score(repairnetfile, 'repairnet')
print 'done.'

combined = draftnet
combined = combined.union(repairnet)
print '\nChecking draftnet + repairnet for unproducible targets ...',
model = pm.get_unproducible(combined, targets, seeds)
print 'done.'
f.write('\nwith repair network still '+str(len(model))+' unproducible targets:')
print_met(f,model.to_list())
never_producible = TermSet()
for a in model :
   target= str(a)[13:]
   t = asp.String2TermSet(target)
   never_producible = never_producible.union(t)

reconstructable_targets = TermSet()
for t in unproducible_targets:
   if not (t in never_producible) :      reconstructable_targets.add(t)
f.write('\n'+str(len(reconstructable_targets))+' targets to reconstruct:')
print_met(f,reconstructable_targets)

#essential_reactions = TermSet()
#for t in reconstructable_targets:
   #single_target = TermSet()
   #single_target.add(t)
   #print '\nComputing essential reactions for',t,'...',
   #essentials =  pm.get_intersection_of_completions(draftnet, repairnet, single_target, seeds)
   #print 'done.'
   #print ' ',len(essentials), 'essential reactions found:'
   #print_met(essentials.to_list())
   #essential_reactions = essential_reactions.union(essentials)
#print '\n  Overall',len(essential_reactions), 'essential reactions found.'

#for r in essential_reactions :
   #nr = str(r)[1:].replace(')',',\"draft\")')
   #draftnet.union(asp.String2TermSet(nr))
    
upperbound=0
lowerbound=0
for t in reconstructable_targets:   
    single_target = TermSet()
    single_target.add(t)
    f.write('\n\nMinimal score for a completion for '+str(t))
    optimum =  pm.get_minimal_completion_size(draftnet, repairnet, single_target, seeds)
    f.write(' ='+str(optimum[0]))
    f.write('\nCompletions with score '+str(optimum[0])+' for '+str(t))
    models =  pm.get_optimal_completions(draftnet, repairnet, single_target, seeds,optimum[0]) 
    count=1
    for m in models:
      f.write('\nsolution '+str(count)+':')
      print_met(f,m.to_list())
      count+=1
    upperbound+=optimum[0]
    if optimum[0] > lowerbound : lowerbound=optimum[0]
#print '  Lowerbound for minimal completion for all targets =',lowerbound
#print '  Upperbound for minimal completion for all targets =',upperbound

#clean_up()
##  use upperbound to reduce repairnet
#print '\nComputing all reactions from completions with at most',upperbound,' score ='


#print '\nComputing minimal score for a completion to produce all targets ...',
#optimum =  pm.get_minimal_completion_size(draftnet, repairnet, unproducible_targets, seeds, upperbound)
#print 'done.'
#print '  minimal score =',optimum[0]

#print '\nComputing essential reactions from all completion with score ',optimum[0],'...',
#model =  pm.get_intersection_of_optimal_completions(draftnet, repairnet, unproducible_targets, seeds, optimum[0])
#print 'done.'
#print_met(model.to_list())

#print '\nComputing essential reactions from all completions ...',
#model =  pm.get_intersection_of_completions(draftnet, repairnet, unproducible_targets, seeds)
#print 'done.'
#print_met(model.to_list())


#optimum =  pm.get_minimal_completion_size(draftnet, repairnet2, targets, seeds)
#print 'done' 

f.close()
clean_up()