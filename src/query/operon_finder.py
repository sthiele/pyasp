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
# along with BioASP.  If not, see <http://www.gnu.org/licenses/>.

# -*- coding: utf-8 -*-from bioasp import asp
from bioasp.asp import *


prepro_prg =      root + '/query/operon_finder/preprocess.lp'
length_prg =    root + '/query/operon_finder/density.lp'
convert_prg =   root + '/query/operon_finder/convert.lp'
ksip_prg =      root + '/query/operon_finder/ksip.lp'
filter_prg =      root + '/query/operon_finder/filter_couples.lp'


def createInstance(genome_string,metabolism_string,catalysation_string) :

  instance = TermSet()

  igenome = open(genome_string, "r")
  genome = igenome.read().replace('\r\n','\n').split()
  dictg = {} # Dictionnary for the genes
  revdictg = [0]
  g = 0 # Counter to define the matchings between real names and ASP names
  for line in genome :
    g+=1
    dictg[line] = g   
    revdictg.append(line)
    #print("gene("+str(g)+").")
    instance.add(Term('gene', [g]))
    
  
  imetabolism = open(metabolism_string, "r")
  metabolism = imetabolism.read().replace('\r\n','\n').splitlines()
  dictr = {}# Dictionnary for the reactions
  revdictr = [0]
  r = 0 # Counter to define the matchings between real names and ASP names
  
  for line in metabolism :
    link = line.split()
    if not dictr.has_key(link[0]) :
      r+=1
      dictr[link[0]] = r
      revdictr.append(link[0])
      #print('\t\t\t<node id="%i" label ="'%(r)+ link[0] + '"/>\n')
    if not dictr.has_key(link[1]) :
      r+=1
      dictr[link[1]] = r
      revdictr.append(link[1])
      #print('\t\t\t<node id="%i" label ="'%(r)+ link[1] + '"/>\n')
    #print "redge("+  dictr[link[0]]+","+dictr[link[1]]+")."
    instance.add(Term('redge', [dictr[link[0]],dictr[link[1]]]))

  #print("\t\t</nodes>\n")

  
  icatalyzis = open(catalysation_string, "r")

  catalyzistab=[]
  i=0
  catalyzis = icatalyzis.readlines()
  for line in catalyzis :
    cat = line.replace('\r\n','\n').split()
  
    if dictr.has_key(cat[0]) and dictg.has_key(cat[1]) :
      #print(cat[0]+"="+dictr[cat[0]])
      #print("cat("+dictg[cat[1]]+","+dictr[cat[0]]+").")
      instance.add(Term('cat', [dictg[cat[1]],dictr[cat[0]]]))
      while i < dictr[cat[0]]+1 :
	catalyzistab.append([])
	i+=1
      catalyzistab[dictr[cat[0]]].append(dictg[cat[1]])  
    #else :
      #print "irrelevant gene or reaction",cat[1],cat[0]

  #print len(catalyzistab)      
  #exit(0)
      
  couplestab = []
  for i in range(len(catalyzistab)) :
    #print i
    for j in range(len(catalyzistab))[(i+1):] :
      #print j
      condition = True
      for g1 in catalyzistab[i] :
	for g2 in catalyzistab[j] :
	  if condition :
	    glen = len(dictg)
	    if min(abs(g1-g2),glen-abs(g1-g2)) <= 10 :
	      couplestab.append([i,j])
	      #print "startgoal(",i,",",j,")"
	      #instance.add(Term('startgoal', [i,j]))
	      condition = False    
	      
  return instance, revdictg, revdictr, couplestab
  
  
def filter_couples(couple_facts, instance):
    prg = [filter_prg, instance, couple_facts.to_file()]
    solver = GringoClasp()
    models = solver.run(prg,nmodels=0,collapseTerms=True, collapseAtoms=False)
    os.unlink(prg[2])
    return models[0]
  
 
def get_ksip_instance(instance, pmax):

    pmaxfact = String2TermSet('pmax('+str(pmax)+')')
  
    inst=pmaxfact.to_file() 
    prg = [convert_prg , instance, inst ]

    solver = GringoClasp()
    solution = solver.run(prg,1)
    os.unlink(inst)
    return solution[0] 
    
def preprocess(instance, start, goal):
    startfact = String2TermSet('start('+str(start)+')')
    goalfact = String2TermSet('goal('+str(goal)+')')
    details = startfact.union(goalfact)
    inst=details.to_file() 
    prg = [prepro_prg , instance, inst ]

    solver = GringoClasp()
    solution = solver.run(prg,1)
    os.unlink(inst)
    return solution[0] 
    
    
def get_k_best_length(instance, start, goal, pmax, k=1, exclude=[]):
    solutions = []
    startfact = String2TermSet('start('+str(start)+')')
    goalfact = String2TermSet('goal('+str(goal)+')')
    pmaxfact = String2TermSet('pmax('+str(pmax)+')')
    
    details = startfact.union(goalfact).union(pmaxfact)
        
    #prepro = preprocess(instance, start, goal)
    #if len(prepro)== 0 : return solutions
    #details = prepro[0].union(startfact).union(goalfact).union(pmaxfact)
    
    
    inst=details.to_file()  
    pmin=1
    
    while len(solutions) < k :
      #pminfact = String2TermSet('pmin('+str(pmin)+')')
      
      prg = [length_prg , instance, inst, exclude_sol(solutions) ]
      groptions=' --const pmin='+str(pmin)
      cloptions='--opt-heu --opt-hier --heu=vsids'
      solver = GringoClaspOpt(gringo_options=groptions,clasp_options=cloptions)
      optima = solver.run(prg,collapseTerms=True,collapseAtoms=False)
   
      os.unlink(prg[3])
      #os.unlink(prg[4])
      
      if optima != None :
	solutions.append(optima[1][0])
	pmin=optima[0][0]
	#if len(solutions) < k :
	  #prg = [length_prg , instance, inst, exclude_sol(solutions) ]
	  #options='--opt-heu --opt-hier --opt-all='+str(optima[0][0])
	  #solver = GringoClasp(clasp_options=options)
	  #models = solver.run(prg,k-len(solutions))  
	  #os.unlink(prg[3])
        
	  #solutions = solutions + models   
      else :
        os.unlink(inst)
	return solutions 
	
    os.unlink(inst)
    return solutions
 
def get_k_best_sip(instance, start, goal, k, exclude=[]):
    solutions = []
    #print "prepro ...",
    test = preprocess(instance, start, goal)
    #print len(test.to_list()),
    prepro = test.to_file()
    
    min = "0"
    while len(solutions) < k :
      prg = [ksip_prg, prepro, exclude_sol(solutions) ]
      goptions='--const min='+str(min)
      
      #coptions='--opt-heu=3 --heu=vsids --opt-hier --restart-on-model --shuffle=1,2' # coli:71, smaller:toolo
      #coptions='--quiet=1,1 --opt-hier'
      #print "search ...",
      coptions='--opt-heu --heu=vsids --opt-hier'
      
      #print "search1 ...",
      solver = GringoClaspOpt(gringo_options=goptions,clasp_options=coptions)
      #solver = GringoUnClasp(gringo_options=goptions,clasp_options=coptions)
      
      optima = solver.run(prg,collapseTerms=False,collapseAtoms=False)
   
      os.unlink(prg[2])
      
      if optima != None :
	#print "found."
	solutions.append(optima[1][0])
	min= optima[0][0]
	
	if len(solutions) < k :
	  #print "search2 ...",
	  #prg = [ksip_prg, prepro, exclude_sol(solutions) ]
	  prg = [ksip_prg, prepro, exclude_sol([optima[1][0]]) ]
	  goptions='--const min='+str(min)
	  coptions='--opt-all='+str(min)
	  solver = GringoClasp(gringo_options=goptions,clasp_options=coptions)
	  models = solver.run(prg,k-len(solutions),collapseTerms=False,collapseAtoms=False)  
	  #models = solver.run(prg,0,collapseTerms=False,collapseAtoms=False)  

	  os.unlink(prg[2])
	  solutions += models   
	  #print "found",len(models),"."
	  min+=1
	  
      else :
	#if len(solutions) == 0 : print "not found."
	os.unlink(prepro)
	return solutions 
    os.unlink(prepro)	 
    #if len(solutions) == 0 : print "not found."
    return solutions 
    