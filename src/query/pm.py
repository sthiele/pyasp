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
from libsbml import *
import sys
import os
import string
from bioasp import asp
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.asp import *

#influence graphs

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

#influence graphs with no change observations
consistency_with_nc_prg = root + '/query/inflgraph_with_no_change_observations/pm_consistency_with_nc.gringo'
prediction_with_nc_prg = root + '/query/inflgraph_with_no_change_observations/pm_prediction_with_nc.gringo'
repair_options_with_nc_prg =     root + '/query/inflgraph_with_no_change_observations/compute_repair_options_with_nc.gringo'
repair_core_with_nc_prg =        root + '/query/inflgraph_with_no_change_observations/repair_core_with_nc.gringo'
prediction_core_with_nc_prg =    root + '/query/inflgraph_with_no_change_observations/prediction_core_with_nc.gringo'
minimize_inconsistencies_prg =        root + '/query/inflgraph_with_no_change_observations/minimize_inconsistencies.gringo'


#metabolic network expansion
unproducible_prg =      root + '/query/networkexpansion/unproducible_targets.lp'
ireaction_prg =         root + '/query/networkexpansion/ireactions.lp'
minimal_extension_prg = root + '/query/networkexpansion/card_min_extensions_all_targets.lp'
extension_prg =         root + '/query/networkexpansion/extensions_all_targets.lp'

#transkriptionfaktor identification
gene_conflicts_prg = root + '/query/missing_tf_problem/conflicts.lp'
non_conflicting_genes_prg = root + '/query/missing_tf_problem/non_conflicting_genes.lp'
conflicting_genes_prg = root + '/query/missing_tf_problem/conflicting_genes.lp'
tf_consistency_prg = root + '/query/missing_tf_problem/consistency.lp'
scc_prg = root + '/query/missing_tf_problem-prob/scc.lp'

#robustness test
obs_generator_prg = root + '/query/robustness/obs_generator.gringo'

#metabolic networks converter
loop_prg = root + '/query/metabolic_networks/loop_detection.lp'
converter_prg = root + '/query/metabolic_networks/metnet2inflgraph.lp'

def generate_obs(seed):
    prg = [obs_generator_prg ]
    options='--rand-freq=1 --seed='+str(seed)
    solver = asp.GringoClasp(clasp_options=options)
    model = solver.run(prg,1)
    return model
    
def get_conflicts(tfbs, network):
    prg = [gene_conflicts_prg, tfbs, network.to_file() ]
    solver = asp.GringoClasp()
    model = solver.run(prg,1)
    os.unlink(prg[2])
    return model[0]
    
def get_conflicting_genes(tfbs, network):
    prg = [conflicting_genes_prg, tfbs, network.to_file() ]
    solver = asp.GringoClasp()
    model = solver.run(prg,1)
    os.unlink(prg[2])
    return model[0]
    
def get_non_conflicting_genes(tfbs, network):
    prg = [non_conflicting_genes_prg, tfbs, network.to_file() ]
    solver = asp.GringoClasp()
    model = solver.run(prg,1)
    os.unlink(prg[2])
    return model[0]  

def get_scc(tfbs, correlationdata):
    prg = [tfbs, correlationdata, scc_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg,0)
    return models
    
def is_repairable(tfbs, network, unknowntf):
    prg = [tf_consistency_prg, tfbs, network.to_file() ]
    options=' --const n='+str(unknowntf)
    options2='--project --heu=vsids'
    solver = asp.GringoClasp(clasp_options=options2, gringo_options=options)
    model = solver.run(prg,1)
    os.unlink(prg[2])
    return model != []
      
def get_cautious_repairs(tfbs, scc, unknowntf):
    prg = [tf_consistency_prg,tfbs, scc.to_file() ]
    #options=' --restart-on-model --opt-heu --opt-value=54'
    #options=' --restart-on-model --opt-heu'
    #options='--opt-heu --restart-on-model --restarts=3'
    options='--cautious --project --opt-ignore'
    options2= '--const n='+str(unknowntf)
    solver = asp.GringoClasp(clasp_options=options, gringo_options=options2)
    model = solver.run(prg,0)
    os.unlink(prg[2])
    return model[0]    

def get_minimum_repair(tfbs, scc, unknowntf):
    prg = [tf_consistency_prg,tfbs, scc.to_file() ]
    #options=' --restart-on-model --opt-heu --opt-value=54'
    #options=' --restart-on-model --opt-heu'
    #options='--opt-heu --restart-on-model --restarts=3'
    options='--restart-on-model --opt-heu --heu=vsids --project --restarts=3'
    options2= '--const n='+str(unknowntf)
    solver = asp.GringoClaspOpt(clasp_options=options, gringo_options=options2)
    optimum = solver.run(prg)
    os.unlink(prg[2])
    return optimum[0]
    
def get_all_optimal_repairs(tfbs, scc, unknowntf, opt):
    prg = [tf_consistency_prg,tfbs, scc.to_file() ]
    options='--restart-on-model --project --opt-all='+str(opt)
    options2= '--const n='+str(unknowntf)
    solver = asp.GringoClasp(clasp_options=options, gringo_options=options2)
    models = solver.run(prg)
    os.unlink(prg[2])
    return models


def get_unproducible(draft, targets, seeds, exclude=[]):
    prg = [unproducible_prg, draft.to_file(), targets.to_file(), seeds.to_file(), asp.exclude_sol(exclude) ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels=0,collapseTerms=True,collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    os.unlink(prg[4])
    return models[0]

def compute_ireactions(termset):
    prg = [ ireaction_prg, termset.to_file() ]
    #print termset
    solver = asp.GringoClasp()
    models = solver.run(prg,0)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]
    
    
def get_extension_minimum(draftnet, repairnet, targets, seeds, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    instance = draftnet.union(draftfact).union(repairnet).union(targets).union(seeds)
    ireactions = compute_ireactions(instance)
    #print instance
    #print ireactions
    
    prg = [minimal_extension_prg , instance.to_file(), ireactions.to_file(), asp.exclude_sol(exclude) ]
    options=' --local-restarts --reset-restarts --heu=vmtf --save-progress --restart-on-model --restarts=16 '
    #options='-t 2 --opt-hierarch=2 '
    #options=' 0 --sat-prepro --back --restart-on-model --restarts=56 --save-progress --berk-h --reset-restarts --reduce-on-restart --local-restarts --trans-ext=all '
    solver = asp.GringoClaspOpt(clasp_options=options)
    optimum = solver.run(prg)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return optimum[0]
    

def get_intersection_of_all_optimal_extensions(draftnet, repairnet, targets, seeds, optimum, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    instance = draftnet.union(draftfact).union(repairnet).union(targets).union(seeds)
    ireactions = compute_ireactions(instance)
    prg = [minimal_extension_prg , instance.to_file(), ireactions.to_file(), asp.exclude_sol(exclude) ]
    options='--cautious --opt-all='+str(optimum)
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models[0]


def get_intersection_of_all_extensions(draftnet, repairnet, targets, seeds, exclude=[]):
    draftfact = asp.String2TermSet('draft("draft")')
    instance = draftnet.union(draftfact).union(repairnet).union(targets).union(seeds)
    ireactions = compute_ireactions(instance)
    #print ireactions
    prg = [extension_prg , instance.to_file(), ireactions.to_file(), asp.exclude_sol(exclude) ]
    options='--cautious --opt-ignore'
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models[0]
    
def convert(termset):
    prg = [ converter_prg, termset.to_file()]
    solver = asp.GringoClasp()
    model = solver.run(prg,1) 
    os.unlink(prg[1])
    return model[0]
    
    
def get_metabolites_in_loops(termset):
    prg = [loop_prg, termset.to_file()]
    options='--brave --quiet'
    solver = asp.GringoClasp(clasp_options=options)
    model = solver.run(prg,0)
    os.unlink(prg[1])
    return model[0]

def is_consistent(termset):
    '''
    [is_consistent(termset)] returns [True] if there exists a consistent extension
    to the system described by the TermSet object [termset].
    '''
    return consistent_models(termset,1) != []
    
def is_consistent_with_nc(termset):
    '''
    [is_consistent(termset)] returns [True] if there exists a consistent extension
    to the system described by the TermSet object [termset].
    '''
    inputs = compute_reductions(termset)
    prg = [ consistency_with_nc_prg, inputs.to_file(), termset.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels=1,)
    os.unlink(prg[1])
    os.unlink(prg[2])
    return models!= []


def consistent_models(termset,nmodels=0,exclude=[]):
    '''
    [consistent_models(termset,nmodels,exclude)] returns a list containing
    [nmodels] TermSet objects representing consistent extensions of the system
    described by the TermSet [termset]. The list [exclude] should contain TermSet objects
    representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    then the list contain all feasible models.
    '''
    inputs = compute_reductions(termset)
    #print termset
    prg = [ consistency_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models
    
#def get_minimal_explanations_opt(termset,nmodels=0,exclude=[]):
    #inputs = compute_reductions(termset)
    #prg = [ explanation_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ]
    #solver = asp.GringoClaspOpt()
    #optimum = solver.run(prg)
    #os.unlink(prg[1])
    #os.unlink(prg[2])
    #os.unlink(prg[3])
    #return optimum[0][0]

#def get_minimal_explanations(termset,optimum,nmodels=0,exclude=[]):
    #inputs = compute_reductions(termset)
    #prg = [ explanation_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ]
    #options='--opt-all='+str(optimum)
    ##options='--opt-all='+str(optimums[0])+','+str(optimums[1])
    #solver = asp.GringoClasp(clasp_options=options)
    #models = solver.run(prg,nmodels)
    #os.unlink(prg[1])
    #os.unlink(prg[2])
    #os.unlink(prg[3])
    #return models    
    

#def get_predictions_under_minimal_explanations(termset,optimum,nmodels=0,exclude=[]):
    #inputs = compute_reductions(termset)
    #prg = [ explanation_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ]
    #options='--cautious --opt-all='+str(optimum)
    ##options='--opt-all='+str(optimums[0])+','+str(optimums[1])
    #solver = asp.GringoClasp(clasp_options=options)
    #models = solver.run(prg,nmodels)
    #os.unlink(prg[1])
    #os.unlink(prg[2])
    #os.unlink(prg[3])
    #return models[0]
    
    
def compute_mic(termset,nmodels=0,exclude=[]):
    '''
    [compute_mic(termset,nmodels,exclude)] returns a list containing
    [nmodels] TermSet objects representing subset minimal inconsistent cores of the system
    described by the TermSet [termset]. The list [exclude] should contain TermSet objects
    representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    then the list contain all feasible models.
    '''
    inputs = compute_reductions(termset)
    prg = [ mic_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ] 
    options=' --heuristic=Vmtf'
    solver = asp.GringoClaspD(options)
    models = solver.run(prg,nmodels)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return models[0]

def compute_mic_dyn(termset,nmodels=0,exclude=[]):
    '''
    [compute_mic(termset,nmodels,exclude)] returns a list containing
    [nmodels] TermSet objects representing subset minimal inconsistent cores of the system
    described by the TermSet [termset]. The list [exclude] should contain TermSet objects
    representing (maybe partial) solutions that are to be avoided. If [nmodels] equals [0]
    then the list contain all feasible models.
    '''
    inputs = compute_reductions(termset)
    prg = [ dyn_mic_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol(exclude) ] 
    #print inputs
    #print termset
    options=' --heuristic=Vmtf'
    solver = asp.GringoClaspD(options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    #print models
    return models

def guess_inputs(termset):
    prg = [ guess_inputs_prg, termset.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,0, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]

def compute_reductions(termset):
    prg = [ reduction_prg, termset.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,0)
    os.unlink(prg[1])
    assert(len(models) == 1)
    return models[0]


def compute_repair_options_flip_obs(termset):
    repair_mode = asp.String2TermSet('repair_v')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]
    
def compute_repair_options_del_edge(termset):
    repair_mode = asp.String2TermSet('repair_d')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_with_nc_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]

def compute_repair_options_flip_edge(termset):
    repair_mode = asp.String2TermSet('repair_e')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]

def compute_repair_options_make_node_input(termset):
    repair_mode = asp.String2TermSet('repair_g')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]
    
      
def compute_repair_options_make_obs_input(termset):
    repair_mode = asp.String2TermSet('repair_i')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg,1)
    os.unlink(prg[0])    
    return models[0]
    
def compute_repair_options_adding_edges(termset):
    repair_mode = asp.String2TermSet('repair_a')
    termset2 = termset.union(repair_mode)
    prg = [ termset2.to_file(), repair_options_prg ]
    solver = asp.GringoClasp()
    models = solver.run(prg)
    os.unlink(prg[0])    
    return models[0]
    


def compute_card_minimum_of_repairs(instance,repair_options,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    prg = [ termset2.to_file(),repair_options.to_file(), asp.exclude_sol(exclude), repair_core_prg, repair_cardinality_prg ]

    solver = asp.GringoClaspOpt()
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2]) 
    return optimum[0]

    
def compute_card_minimum_inconsistencies_after_delete_with_nc(instance,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), asp.exclude_sol(exclude), minimize_inconsistencies_prg ]

    solver = asp.GringoClaspOpt(clasp_options='--restart-on-model --solu --save-pro --heu=vmtf --opt-heu ')
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return optimum[0]
    
def compute_union_of_minimal_inconsistencies_after_delete_with_nc(instance,optimum,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), asp.exclude_sol(exclude), minimize_inconsistencies_prg ]

    solver = asp.GringoClasp(clasp_options=' --brave --restart-on-model --solu --save-pro --heu=vmtf --project --opt-all='+str(optimum))
    model = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return model[0]
    

def compute_card_minimal_repair_sets(optimum, instance, repair_options ,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    prg = [ termset2.to_file(), repair_options.to_file(), asp.exclude_sol(exclude), repair_core_prg, repair_cardinality_prg ]

    options='--project --opt-all='+str(optimum)
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2])
    return models
    
    
def compute_optimum_of_repairs_with_nc(instance,errornodes):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    #print errornodes
    prg = [ termset2.to_file(),errornodes.to_file(), repair_core_with_nc_prg]

    solver = asp.GringoClaspOpt(clasp_options='--restart-on-model --solu --opt-heu --save-pro  --opt-hier --heu=vmtf')
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return optimum[0]

def compute_union_of_optimal_repair_sets_with_nc(instance,errornodes,optimum):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(),errornodes.to_file(), repair_core_with_nc_prg]

    solver = asp.GringoClasp(clasp_options='--brave --restart-on-model --solu --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    model = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return model[0]
    
#def compute_intersection_of_optimal_repair_sets_with_nc(instance,errornodes,optimum):
    #inputs = compute_reductions(instance)
    #termset2 = instance.union(inputs)
    ##print termset2
    ##print repair_options
    #prg = [ termset2.to_file(),errornodes.to_file(), repair_core_with_nc_prg]

    #solver = asp.GringoClasp(clasp_options='--cautious --restart-on-model --solu --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    #model = solver.run(prg)
    #os.unlink(prg[0])
    #os.unlink(prg[1])
    #return model[0]     
    
def compute_all_optimal_repair_sets_with_nc(instance,errornodes,optimum,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    #print repair_options
    prg = [ termset2.to_file(),errornodes.to_file(), asp.exclude_sol(exclude), repair_core_with_nc_prg]

    solver = asp.GringoClasp(clasp_options='--project --restart-on-model --solu --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    models = solver.run(prg,nmodels=0, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2])
    return models 
    
    
def compute_distribution_over_all_optimal_repair_sets_with_nc(instance,errornodes,optimum,nmodels=0,exclude=[]):
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)

    prg = [ termset2.to_file(),errornodes.to_file(), asp.exclude_sol(exclude), repair_core_with_nc_prg]

    solver = asp.GringoClaspDistr(clasp_options='--project --restart-on-model --solu --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    distribution = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2])
    return distribution 
     
def prediction_under_cardinality_minimal_repair(optimum, instance, repair_options):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    inputs = compute_reductions(instance)
    termset2 = instance.union(inputs)
    
    prg = [ termset2.to_file(), repair_options.to_file(), prediction_core_prg, repair_cardinality_prg ]

    options='--project --cautious --opt-all='+str(optimum)
    #print options
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return whatsnew(instance,models[0])


def prediction_under_cardinality_minimal_repair_with_nc(optimum, instance, repair_options):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    #inputs = compute_reductions(instance)
    termset2 = instance#.union(inputs)

    prg = [ termset2.to_file(), repair_options.to_file(), prediction_core_with_nc_prg]

    options='--project --cautious --opt-all='+str(optimum[0])+","+str(optimum[1])
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return whatsnew(instance,models[0])

def cut_obs_(s):
    return str(s)[4:]

def whatsnew(termset,pred):
    '''
    [whatsnew(termset,pred)] is a TermSet equal to [pred] where all predicates
    vlabel and elabel which have a corresponding obs_vlabel and obs_elabel in
    [termset] have been deleted. This function is meant to see which of the invariants
    are not a direct consequence of the observations.
    '''
    accu = asp.TermSet(pred)
    for t in termset:

        if t.pred() == 'obs_vlabel':
            [_,e,v,s] = t.explode()
            accu.discard(asp.Term('vlabel',[e,v,s]))
        elif t.p('obs_elabel'):
            [_,j,i,s] = t.explode()
            accu.discard(asp.Term('elabel',[j,i,s]))
    return accu


def prediction_under_consistency(termset):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    inputs = compute_reductions(termset)

    prg = [ prediction_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol([]) ]
    solver = asp.GringoClasp(clasp_options='--project --cautious')
    models = solver.run(prg,0)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return whatsnew(termset,models[0])
    
def prediction_under_consistency_with_nc(termset):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    #inputs = compute_reductions(termset)
    inputs =TermSet()

    prg = [ prediction_with_nc_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol([]) ]
    solver = asp.GringoClasp(clasp_options='--project --cautious')
    models = solver.run(prg,0)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return whatsnew(termset,models[0])
