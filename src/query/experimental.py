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
from bioasp.query import influence_graphs as ig
from bioasp.misc import *
import tempfile
from bioasp.data import op
from bioasp.asp import *


#influence graphs with no change observations

graph_compress_prg = root + '/query/inflgraph_with_no_change_observations/compress.lp'
consistency_with_nc_prg = root + '/query/inflgraph_with_no_change_observations/pm_consistency_with_nc.gringo'
prediction_with_nc_prg = root + '/query/inflgraph_with_no_change_observations/pm_prediction_with_nc.gringo'
repair_options_with_nc_prg =     root + '/query/inflgraph_with_no_change_observations/compute_repair_options_with_nc.gringo'
repair_core_with_nc_prg =        root + '/query/inflgraph_with_no_change_observations/repair_core_with_nc.gringo'
prediction_core_with_nc_prg =    root + '/query/inflgraph_with_no_change_observations/prediction_core_with_nc.gringo'
minimize_inconsistencies_prg =        root + '/query/inflgraph_with_no_change_observations/minimize_inconsistencies.gringo'
card_inconsistencies_pierre_prg =        root + '/query/inflgraph_with_no_change_observations/card_inconsistencies_for_pierre.gringo'


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


def is_consistent_with_nc(termset):
    '''
    [is_consistent(termset)] returns [True] if there exists a consistent extension
    to the system described by the TermSet object [termset].
    '''
    inputs = ig.get_reductions(termset)
    prg = [ consistency_with_nc_prg, inputs.to_file(), termset.to_file() ]
    solver = asp.GringoClasp()
    models = solver.run(prg,nmodels=1,)
    os.unlink(prg[1])
    os.unlink(prg[2])
    return models!= []
    
def get_compressed_graph(instance):    
    prg = [ instance.to_file(), graph_compress_prg]
    solver = asp.GringoClasp()
    model = solver.run(prg, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    compressed_graph = TermSet()
    for atom in model[0]:
      if atom.pred() == "compress_obs_elabel" :
        compressed_graph.add(Term('obs_elabel',atom.args())) 
      if atom.pred() == "compress_edge" :
        compressed_graph.add(Term('edge',atom.args())) 
      if atom.pred() == "compress_vertex" :
        compressed_graph.add(Term('vertex',atom.args())) 
      if atom.pred() == "input" :
        compressed_graph.add(atom) 
      if atom.pred() == "obs_vlabel" :
        compressed_graph.add(atom) 
    return compressed_graph    

    
def compute_card_minimum_inconsistencies_after_delete_with_nc(instance,nmodels=0,exclude=[]):
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), asp.exclude_sol(exclude), minimize_inconsistencies_prg ]

    solver = asp.GringoClaspOpt(clasp_options='--restart-on-model --solu --save-pro --heu=vmtf --opt-heu ')
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return optimum[0]
    
def compute_union_of_minimal_inconsistencies_after_delete_with_nc(instance,optimum,nmodels=0,exclude=[]):
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), asp.exclude_sol(exclude), minimize_inconsistencies_prg ]
    solver = asp.GringoClasp(clasp_options=' --brave --restart-on-model --save-pro --heu=vmtf --project --opt-all='+str(optimum))
    model = solver.run(prg,nmodels=0)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return model[0]
    
def compute_card_inconsistencies_for_pierre(instance,nmodels=0,exclude=[]):
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), asp.exclude_sol(exclude), card_inconsistencies_pierre_prg ]

    solver = asp.GringoClaspOpt(clasp_options='--restart-on-model --solu --save-pro --heu=vmtf --opt-heu ')
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return optimum[0]
     
    
def compute_optimum_of_repairs_with_nc(instance,errornodes):
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    #print errornodes
    prg = [ termset2.to_file(), errornodes.to_file(), repair_core_with_nc_prg]

    solver = asp.GringoClaspOpt(clasp_options='--restart-on-model --solu --opt-heu --opt-hier --heu=vmtf') #--save-pro
    optimum = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    print optimum
    return optimum[0]

def compute_union_of_optimal_repair_sets_with_nc(instance, errornodes, optimum):
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)
    #print termset2
    prg = [ termset2.to_file(), errornodes.to_file(), repair_core_with_nc_prg]
    solver = asp.GringoClasp(clasp_options='--brave --restart-on-model --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    model = solver.run(prg,nmodels=0)
    os.unlink(prg[0])
    os.unlink(prg[1])
    return model[0]
    
def compute_all_optimal_repair_sets_with_nc(instance,errornodes,optimum,nmodels=0,exclude=[]):
    inputs = ig.get_reductions(instance)
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
    inputs = ig.get_reductions(instance)
    termset2 = instance.union(inputs)

    prg = [ termset2.to_file(),errornodes.to_file(), asp.exclude_sol(exclude), repair_core_with_nc_prg]

    solver = asp.GringoClaspDistr(clasp_options='--project --restart-on-model --solu --save-pro --heu=vmtf --opt-all='+str(optimum[0])+','+str(optimum[1]))
    distribution = solver.run(prg)
    os.unlink(prg[0])
    os.unlink(prg[1])
    os.unlink(prg[2])
    return distribution 


def prediction_under_cardinality_minimal_repair_with_nc(instance,errornodes,optimum):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    #inputs = get_reductions(instance)
    termset2 = instance#.union(inputs)

    prg = [ termset2.to_file(), errornodes.to_file(), prediction_core_with_nc_prg]

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

def prediction_under_consistency_with_nc(termset):
    '''
    Computes the set of signs on edges/vertices that can be cautiously
    derived from [termset], minus those that are a direct consequence
    of obs_[ev]label predicates
    '''
    #inputs = get_reductions(termset)
    inputs =TermSet()

    prg = [ prediction_with_nc_prg, inputs.to_file(), termset.to_file(), asp.exclude_sol([]) ]
    solver = asp.GringoClasp(clasp_options='--project --cautious')
    models = solver.run(prg,0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[1])
    os.unlink(prg[2])
    os.unlink(prg[3])
    return whatsnew(termset,models[0])
