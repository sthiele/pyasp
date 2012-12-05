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

functions_prg =     root + '/query/psnoptimization/functions.lp'
hyper_prg =         root + '/query/psnoptimization/hyper.lp'
optimization_prg =  root + '/query/psnoptimization/optimization.lp'
output_prg =        root + '/query/psnoptimization/output.lp'
compress_prg =      root + '/query/psnoptimization/compress.lp'

def get_minimal_model_size(instance, p, q):
    prg = [ instance.to_file(), hyper_prg,functions_prg,optimization_prg,output_prg ]
    goptions='--const p=%s --const q=%s' % (p,q)
    coptions=''
    solver = asp.GringoClaspOpt(gringo_options=goptions,clasp_options=coptions)
    optimum = solver.run(prg, collapseAtoms=False)
    os.unlink(prg[0])
    #print optimum[1]
    return optimum
    
    
def get_minimal_models(instance, optimum, p, q):
    prg = [ instance.to_file(), hyper_prg,functions_prg,optimization_prg,output_prg ]
    goptions='--const p=%s --const q=%s' % (p,q)
    coptions=' --opt-all='+str(optimum)
    solver = asp.GringoClasp(gringo_options=goptions,clasp_options=coptions)
    models = solver.run(prg, nmodels=0, collapseTerms=False, collapseAtoms=False)
    os.unlink(prg[0])

    return models 

def get_suboptimal_models(instance, suboptimal, size, p, q):
    prg = [ instance.to_file(), hyper_prg,functions_prg,optimization_prg,output_prg ]
    goptions='--const p=%s --const q=%s --const maxsize=%s' % (p,q,str(size))
    coptions=' --opt-all='+str(suboptimal)
    solver = asp.GringoClasp(gringo_options=goptions,clasp_options=coptions)
    models = solver.run(prg, nmodels=0, collapseTerms=False, collapseAtoms=False)
    os.unlink(prg[0])
    return models 
    
    
def get_compressed_network(instance):
  
    prg = [ instance.to_file(), compress_prg ]
    options=''
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg, nmodels=1, collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[0])
    compnet= TermSet()
    for a in models[0]:
        if a.pred().find('compress_')!= -1 :
           if a.pred() == "compress_edge" :  
              compnet.add(Term('edge',[ a.arg(0), a.arg(1), a.arg(2)]))
           if a.pred() == "compress_vertex":  
              compnet.add(Term('vertex',[ a.arg(0)]))
        else : compnet.add(a)
        
    return compnet
    
def get_hypergraph(instance):
  
    prg = [ instance.to_file(), hyper_prg ]   
    options=''
    solver = asp.GringoClasp(clasp_options=options)
    models = solver.run(prg, nmodels=1, collapseTerms=False, collapseAtoms=False)    
    return models[0]
    
    
def remove(x,pe,pbe,ne,nbe):
    #update the mappings following the signs rules
    for y in pbe[x]:
        if x in pe[y]:
            pe[y].remove(x)
        if x in ne[y]:
            ne[y].remove(x)
        
        pe[y] = pe[y].union(pe[x])
        ne[y] = ne[y].union(ne[x])

    for y in nbe[x]:
        if x in pe[y]:
            pe[y].remove(x)
        if x in ne[y]:
            ne[y].remove(x)
            
        ne[y] = ne[y].union(pe[x])
        pe[y] = pe[y].union(ne[x])

    for y in pe[x]:
        if x in pbe[y]:
            pbe[y].remove(x)
        if x in nbe[y]:
            nbe[y].remove(x)
            
        pbe[y] = pbe[y].union(pbe[x])
        nbe[y] = nbe[y].union(nbe[x])    
    
    for y in ne[x]:
        if x in pbe[y]:
            pbe[y].remove(x)
        if x in nbe[y]:
            nbe[y].remove(x)
            
        nbe[y] = nbe[y].union(pbe[x])
        pbe[y] = pbe[y].union(nbe[x])
    
    #unset the mappings for x
    pe[x] = set()
    pbe[x] = set()
    ne[x] = set()
    nbe[x] = set()


def reduce_network(tc, pe, pbe, ne, nbe, rn = set()):
    rnodes = set()

    #go through the nodes that could be removed.
    for x in tc:
        if len(pbe[x]) + len(nbe[x]) <= 1:
            remove(x,pe,pbe,ne,nbe)
            rnodes.add(x)
            
        elif len(pe[x]) + len(ne[x]) <= 1:
            remove(x,pe,pbe,ne,nbe)
            rnodes.add(x)
            
    if len(rnodes) == 0:
        #fix point
        return rn
    else:
        #run another round of reduction
        return reduce_network(tc.difference(rnodes), pe, pbe, ne, nbe, rnodes.union(rn))

           
def compress(instance):
  
    compressed_net = TermSet()
    
    #mappings of edges forwards and backwards, postite and negative
    pos_edges = {}
    pos_back_edges = {}
    neg_edges = {}
    neg_back_edges = {}


    desi = []
    nodes = []
    for a in instance:
      if a.pred() == "stimuli" :
        desi.append(a.arg(0))
        compressed_net.add(a)
      if a.pred() == "inhibitor" :
        desi.append(a.arg(0))
        compressed_net.add(a)
      if a.pred() == "readout" :
        desi.append(a.arg(0))
        compressed_net.add(a)
      if a.pred() == "vertex" :
	nodes.append(a.arg(0))
    designated = set(desi)
    to_compress = set(nodes).difference(desi)
    #designated = network.types['stimuli']+network.types['inhibitor']+network.types['readout']
    #to_compress = network.nodes.difference(designated)
    #print "designated:", designated
    #print "to_compress:",to_compress
    
    clinks = []
    
    for x in nodes:
    #for x in network.nodes:
        #initialization
        pos_edges[x] = set()
        pos_back_edges[x] = set()
        neg_edges[x] = set()
        neg_back_edges[x] = set()
        
    #read all the edges in the network and build the mappings    
    for a in instance:
      if a.pred() == "edge" :  
        if a.arg(2) == "1":
	  pos_edges[a.arg(0)].add(a.arg(1))
	  pos_back_edges[a.arg(1)].add(a.arg(0))
          #print "+",a
        if a.arg(2) == "-1":
	  neg_edges[a.arg(0)].add(a.arg(1))
	  neg_back_edges[a.arg(1)].add(a.arg(0))
          #print "-",a    
        
    #read all the edges in the network and build the mappings
    #for s,y in network.links:
        ##assuming that there are not hyper-edges, s must be a list with 1 tuple
        #x,r = s[0]
        #if r == 1:
            #pos_edges[x].add(y)
            #pos_back_edges[y].add(x)
        #else:
            #neg_edges[x].add(y)
            #neg_back_edges[y].add(x)
            
            
            

    #reduce the network until reaching a fix point
    rnodes = reduce_network(to_compress,pos_edges,pos_back_edges,neg_edges,neg_back_edges)        
    

    
    #rebuild the list of remaining vertices in ASP format
    for v in rnodes:
        compressed_net.add(Term('vertex', [v]))
    for v in desi:
        compressed_net.add(Term('vertex', [v]))
    
    #rebuild the list of the compressed edges in ASP format
    for s,ts in pos_edges.iteritems():
        for t in ts:
            compressed_net.add(Term('edge', [s,t,1]))
    
    for s,ts in neg_edges.iteritems():
        for t in ts:
            compressed_net.add(Term('edge', [s,t,-1]))
            
            
    
    #rebuild the list of the compressed edges in the BooleanModel format
    #for s,ts in pos_edges.iteritems():
        #for t in ts:
            #clinks.append(([(s,1)],t))
    
    #for s,ts in neg_edges.iteritems():
        #for t in ts:
            #clinks.append(([(s,-1)],t))
            
    #return list of edges and removed nodes.

    return compressed_net
