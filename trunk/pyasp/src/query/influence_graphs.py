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
import sys
import os
import string
from bioasp import asp
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.asp import *

guess_inputs_prg = root + '/query/inflgraph/pm_guess_inputs.gringo'

consistency_prg = root + '/query/inflgraph/pm_consistency.gringo'
prediction_prg = root + '/query/inflgraph/pm_prediction.gringo'
explanation_prg = root + '/query/inflgraph/explanation.gringo'
mic_prg = root + '/query/inflgraph/pm_mic.gringo'
dyn_mic_prg = root + '/query/inflgraph/dynamic_diagnosis.gringo'
reduction_prg = root + '/query/inflgraph/reduction.lp'

repair_options_prg =     root + '/query/inflgraph/compute_repair_options.gringo'
repair_core_prg =        root + '/query/inflgraph/repair_core.gringo'
prediction_core_prg =    root + '/query/inflgraph/prediction_core.gringo'
repair_cardinality_prg = root + '/query/inflgraph/repair_cardinality.gringo'
repair_subset_prg  =     root + '/query/inflgraph/repair_subset.gringo'

def is_consistent(instance):
    '''
    [is_consistent(instance)] returns [True] if there exists a consistent extension
    to the system described by the TermSet object [instance].
    '''
    return get_consistent_labelings(instance,1) != []
    

def get_consistent_labelings(instance,nmodels=0,exclude=[]):
    '''
    [consistent_labelings(instance,nmodels,exclude)] returns a list containing
    [nmodels] TermSet objects representing consistent extensions of the system
    described by the TermSet [instance]. The list [exclude] should contain TermSet objects
    representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    then the list contain all feasible models.
    '''
    inputs = get_reductions(instance)
    prg = [ consistency_prg, inputs.to_file(), instance.to_file(), asp.exclude_sol(exclude) ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models
   
    
#def get_minimal_inconsistent_cores(instance,nmodels=0,exclude=[]):
    #'''
    #[compute_mic(instance,nmodels,exclude)] returns a list containing
    #[nmodels] TermSet objects representing subset minimal inconsistent cores of the system
    #described by the TermSet [instance]. The list [exclude] should contain TermSet objects
    #representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    #then the list contain all feasible models.
    #'''
    #inputs = get_reductions(instance)
    #prg = [ mic_prg, inputs.to_file(), instance.to_file(), asp.exclude_sol(exclude) ] 
    #options=' --heuristic=Vmtf'
    #solver = asp.GringoClaspD(options)
    #models = solver.run(prg,nmodels)
    #os.unlink(prg[1])
    #os.unlink(prg[2])
    #os.unlink(prg[3])
    #return models[0]

def get_minimal_inconsistent_cores(instance,nmodels=0,exclude=[]):
    '''
    [compute_mic(instance,nmodels,exclude)] returns a list containing
    [nmodels] TermSet objects representing subset minimal inconsistent cores of the system
    described by the TermSet [instance]. The list [exclude] should contain TermSet objects
    representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    then the list contain all feasible models.
    '''
    inputs = get_reductions(instance)
    prg = [ dyn_mic_prg, inputs.to_file(), instance.to_file(), asp.exclude_sol(exclude) ] 
    options=' --heuristic=Vmtf'
    solver = asp.GringoClaspD(options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models

def guess_inputs(instance):
    prg = [ guess_inputs_prg, instance.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,0, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]

def get_reductions(instance):
    prg = [ reduction_prg, instance.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,0)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]


def get_repair_options_flip_obs(instance):
    repair_mode = asp.String2TermSet('repair_v')
    instance2 = instance.union(repair_mode)
    prg = [ instance2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]
    

def get_repair_options_flip_edge(instance):
    repair_mode = asp.String2TermSet('repair_e')
    instance2 = instance.union(repair_mode)
    prg = [ instance2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]

def get_repair_options_make_node_input(instance):
    repair_mode = asp.String2TermSet('repair_g')
    instance2 = instance.union(repair_mode)
    prg = [ instance2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]
    
      
def get_repair_options_make_obs_input(instance):
    repair_mode = asp.String2TermSet('repair_i')
    instance2 = instance.union(repair_mode)
    prg = [ instance2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg,1)
    os.unlink(prg[0])    
    return models[0]
    
def get_repair_options_add_edges(instance):
    repair_mode = asp.String2TermSet('repair_a')
    instance2 = instance.union(repair_mode)
    prg = [ instance2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]


def get_minimum_of_repairs(instance,repair_options,exclude=[]):
    inputs = get_reductions(instance)
    instance2 = instance.union(inputs)
    prg = [ instance2.to_file(),repair_options.to_file(), asp.exclude_sol(exclude), repair_core_prg, repair_cardinality_prg ]

    solver = asp.GringoClaspOpt()
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2]) 
    return optimum[0]


def get_minimal_repair_sets(instance, repair_options ,optimum,nmodels=0,exclude=[]):
    inputs = get_reductions(instance)
    instance2 = instance.union(inputs)
    prg = [ instance2.to_file(), repair_options.to_file(), asp.exclude_sol(exclude), repair_core_prg, repair_cardinality_prg ]

    options='--project --opt-all='+str(optimum)
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2])
    return models
    
     
def get_predictions_under_minimal_repair(instance, repair_options, optimum):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [instance], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    inputs = get_reductions(instance)
    instance2 = instance.union(inputs)
    
    prg = [ instance2.to_file(), repair_options.to_file(), prediction_core_prg, repair_cardinality_prg ]

    options='--project --enum-mode cautious --opt-all='+str(optimum)
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return whatsnew(instance,models[0])


def cut_obs_(s):
    return str(s)[4:]

def whatsnew(instance,pred):
    '''
    [whatsnew(instance,pred)] is a TermSet equal to [pred] where all predicates
    vlabel and elabel which have a corresponding obs_vlabel and obs_elabel in
    [instance] have been deleted. This function is meant to see which of the invariants
    are not a direct consequence of the observations.
    '''
    accu = asp.TermSet(pred)
    for t in instance:

        if t.pred() == 'obs_vlabel':
            [_,e,v,s] = t.explode()
            accu.discard(asp.Term('vlabel',[e,v,s]))
        elif t.p('obs_elabel'):
            [_,j,i,s] = t.explode()
            accu.discard(asp.Term('elabel',[j,i,s]))
    return accu


def get_predictions_under_consistency(instance):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [instance], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    inputs = get_reductions(instance)

    prg = [ prediction_prg, inputs.to_file(), instance.to_file(), asp.exclude_sol([]) ]
    solver = asp.GringoClasp(clasp_options='--project --enum-mode cautious')
    models = solver.run(prg,0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return whatsnew(instance,models[0])
    