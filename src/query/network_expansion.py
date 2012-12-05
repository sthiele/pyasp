# -*- coding: utf-8 -*-
import sys
import os
import string
from bioasp import asp
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.asp import *

unproducible_prg =      root + '/query/networkexpansion/unproducible_targets.lp'
ireaction_prg =         root + '/query/networkexpansion/ireactions.lp'
minimal_completion_prg = root + '/query/networkexpansion/card_min_completions_all_targets.lp'
minimal_completion_wb_prg = root + '/query/networkexpansion/card_min_completions_all_targets_with_bounds.lp'
completion_prg =         root + '/query/networkexpansion/completions_all_targets.lp'


def get_mapping_ireaction(termset):
    dict={}
    revdict={}
    for a in termset:
      if a.pred() == "ireaction" :
        if not a.arg(0) in dict:
	  id=len(dict)
          dict[a.arg(0)]=id
          revdict[id]=a.arg(0)
          #print id," ",a.arg(0)
    #print revdict
    #quit()
    return dict, revdict
          
def map_reaction_ids(termset, dict):
    mapped = TermSet()
    for a in termset:
      if a.pred() == "reaction" :
        if a.arg(0) in dict:
          mapped.add(Term('reaction', [str(dict[a.arg(0)]), a.arg(1)]))
        else : mapped.add(a)
        
      elif a.pred() == "xreaction" :
        if a.arg(0) in dict:
          mapped.add(Term('xreaction', [str(dict[a.arg(0)]), a.arg(1)]))
        else : mapped.add(a)
        
      elif a.pred() == "ireaction" :
        if a.arg(0) in dict:
          mapped.add(Term('ireaction', [str(dict[a.arg(0)]), a.arg(1)]))
        else : print "Error: unknown ireaction, networkexpansion.py line 39"
        
      elif a.pred() == "value" :
        if a.arg(0) in dict:
          mapped.add(Term('value', [str(dict[a.arg(0)]), a.arg(1)]))
        else : mapped.add(a)     
        
      elif a.pred() == "product" :
	if  a.arg(1) in dict:
          mapped.add(Term('product', [a.arg(0), str(dict[a.arg(1)]),a.arg(2)]))
        else : mapped.add(a)
        
      elif a.pred() == "reactant" :
	if a.arg(1) in dict:
          mapped.add(Term('reactant', [a.arg(0), str(dict[a.arg(1)]),a.arg(2)]))
        else : mapped.add(a)  
      else :
	 mapped.add(a)
	
    return mapped
    
    
def unmap_reaction_ids(termset, revdict):
    unmapped = TermSet()
    for a in termset:
      if a.pred() == "xreaction" :
	unmapped.add(Term('xreaction', [str(revdict[int(a.arg(0))]), a.arg(1)]))

    return unmapped



def get_unproducible(draft, targets, seeds, exclude=[]):
    prg = [unproducible_prg, draft.to_file(), targets.to_file(), seeds.to_file(), asp.exclude_sol(exclude) ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels=0,collapseTerms=True,collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    os.unlink(prg[4])
    return models[0]

def compute_ireactions(instance):
    prg = [ ireaction_prg, instance.to_file() ]
    #print instance
    solver = asp.GringoClasp()
    models = solver.run(prg,0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]
    
    
def get_minimal_completion_size(draft, repairnet, targets, seeds, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    instance = draft.union(draftfact).union(repairnet).union(targets).union(seeds)
           
    ireactions = compute_ireactions(instance)
    #print ireactions
    
    instance = instance.union(ireactions)

    prg = [minimal_completion_prg , instance.to_file(), asp.exclude_sol(exclude) ]

    options=' --local-restarts --reset-restarts --opt-heu --opt-hier --save-progress --restart-on-model --restarts=16 '
   
    solver = asp.GringoClaspOpt(clasp_options=options)
    optimum = solver.run(prg,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    return optimum
    
def get_minimal_completion_size_with_bounds(draft, repairnet, targets, seeds, upperbound=200, lowerbound=0, exclude=[]):
    
    draftfact = asp.String2TermSet('draft("draft")')
    ubfact = asp.String2TermSet('ub('+str(upperbound)+')')
    lbfact = asp.String2TermSet('lb('+str(lowerbound)+')')
    instance = draft.union(draftfact).union(ubfact).union(lbfact).union(repairnet).union(targets).union(seeds)
           
    ireactions = compute_ireactions(instance)
    #print ireactions
    dict, revdict = get_mapping_ireaction(ireactions)
    
    instance = map_reaction_ids(instance.union(ireactions), dict)

    prg = [minimal_completion_wb_prg , instance.to_file(), asp.exclude_sol(exclude) ]

    options=' --local-restarts --reset-restarts --opt-heu --opt-hier --save-progress --restart-on-model --restarts=16 --opt-val='+str(upperbound)
    #options='-t 2 --opt-hierarch=2 '
    #options=' 0 --sat-prepro --back --restart-on-model --restarts=56 --save-progress --berk-h --reset-restarts --reduce-on-restart --local-restarts --trans-ext=all '
    solver = asp.GringoClaspOpt(clasp_options=options)
    optimum = solver.run(prg)
    os.unlink(prg[1])
    os.unlink(prg[2])
    return optimum
    

def get_intersection_of_optimal_completions(draft, repairnet, targets, seeds, optimum, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    ubfact = asp.String2TermSet('ub('+str(optimum)+')')
    lbfact = asp.String2TermSet('lb('+str(optimum)+')')
    instance = draft.union(draftfact).union(ubfact).union(lbfact).union(repairnet).union(targets).union(seeds)
    
    ireactions = compute_ireactions(instance)
    #print ireactions
    dict, revdict = get_mapping_ireaction(ireactions)
    instance = map_reaction_ids(instance.union(ireactions), dict)
    
    instance.to_file("instance.lp")
    prg = [minimal_completion_wb_prg , instance.to_file(), asp.exclude_sol(exclude) ]
    options='--enum-mode cautious '
    #options='--cautious --opt-all='+str(optimum)
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    
    solution = unmap_reaction_ids(models[0], revdict)
    return solution
    
def get_optimal_completions(draft, repairnet, targets, seeds, optimum, nmodels=0, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    ubfact = asp.String2TermSet('ub('+str(optimum)+')')
    lbfact = asp.String2TermSet('lb('+str(optimum)+')')
    instance = draft.union(draftfact).union(ubfact).union(lbfact).union(repairnet).union(targets).union(seeds)
    
    ireactions = compute_ireactions(instance)
    #print ireactions
    dict, revdict = get_mapping_ireaction(ireactions)
    instance = map_reaction_ids(instance.union(ireactions), dict)    
    
    prg = [minimal_completion_wb_prg , instance.to_file(), asp.exclude_sol(exclude) ]
    #options=' --opt-all='+str(optimum)
    options=''
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    
    solutions= []
    for m in models :
       s = unmap_reaction_ids(m, revdict)
       solutions.append(s)
    return solutions


def get_intersection_of_completions(draft, repairnet, targets, seeds, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    instance = draft.union(draftfact).union(repairnet).union(targets).union(seeds)
    ireactions = compute_ireactions(instance)
    #print ireactions
    instance = instance.union(ireactions)
    #dict, revdict = get_mapping_ireaction(ireactions)
    #instance = map_reaction_ids(instance.union(ireactions), dict)  
    
    prg = [completion_prg , instance.to_file(), asp.exclude_sol(exclude) ]
    options='--enum-mode cautious --opt-ignore'
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    
    #solution = unmap_reaction_ids(models[0], revdict)
    #return solution
    
    return models[0]
    