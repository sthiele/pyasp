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
from bioasp.asp import *

from bioasp.query import operon_finder as of 

import sys
import os
from bioasp.misc import *
import tempfile
from bioasp.data import op

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
	      
  return instance, revdictg, revdictr, dictr, couplestab

def readcouples(couple_string,dictr, revdictr) :

  icouples = open(couple_string, "r")
  bla =  icouples.read().replace('\r\n','\n').splitlines()
  #dictr = {}# Dictionnary for the reactions
  #revdictr = [0]
  r = len(revdictr)
  #r= 0# Counter to define the matchings between real names and ASP names
  couples = []
  for line in bla :
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
    #print dictr[link[0]],
    #print dictr[link[1]]
    couples.append([dictr[link[0]],dictr[link[1]]])
  return couples, revdictr

     



def print_reactions(answer,dictr,dictg) :
  for a in answer:
    if a.pred() == "fedge" :
      print "\n   ",dictr[int(a.arg(0).arg(1))]+"("+dictg[int(a.arg(0).arg(0))]+")",
      print"->",
      print dictr[int(a.arg(1).arg(1))]+"("+dictg[int(a.arg(1).arg(0))]+")",
      print ": ",a.arg(2),
      
  print
  
      
  
def clean_up() :
  if os.path.isfile("parser.out"): os.remove("parser.out")
  if os.path.isfile("asp_py_lextab.py"): os.remove("asp_py_lextab.py")
  if os.path.isfile("asp_py_lextab.pyc"): os.remove("asp_py_lextab.pyc")
  if os.path.isfile("asp_py_parsetab.py"): os.remove("asp_py_parsetab.py")
  if os.path.isfile("asp_py_parsetab.pyc"): os.remove("asp_py_parsetab.pyc") 


if len(sys.argv)==6 :
  genome_string=sys.argv[1]
  metabolism_string=sys.argv[2]
  catalysation_string=sys.argv[3]
  k=sys.argv[4]
  couple_string =sys.argv[5]
  
  
else:
  print "\nmissing arguments !!!"
  print "usage: "
  print "   python finder_ksip.py  genomefile metabolismfile catalysationfile number_of_sips couplesfile"
  print "\nexample:"
  print "   python finder_ksip.py data/genome.txt data/metabolism.txt data/catalyze.txt 5 couples\n"
  exit(0)
  
instance, dictg, revdictr, dictr, oldcouples = createInstance(genome_string,metabolism_string,catalysation_string) 
inst=instance.to_file()
ksip_instance = of.get_ksip_instance(inst,300)
kinst=ksip_instance.to_file()

couples, revdictr = readcouples(couple_string, dictr, revdictr)

couple_facts = TermSet()
for s,g in couples: 
   couple_facts.add(Term('pair', [str(s),str(g)]))
print "filter queries ...",    
filter_couples = of.filter_couples(couple_facts,kinst)
print "done." 

new_couples = []
for a in filter_couples :
  new_couples.append([int(a.arg(0)), int(a.arg(1))] )

#new_couples=[[829,828]]
for s,g in new_couples: 
  
  #print s,g
  models = of.get_k_best_sip(kinst, s, g, int(k))
 
  if len(models) > 0 :
    print "\n"+str(k)+" shortest paths from reaction",revdictr[s],"to",revdictr[g]
  c=1
  for m in models :
     print " ",str(c)+":",
     c+=1
     print_reactions(m, revdictr,dictg)

     
os.unlink(inst)
os.unlink(kinst)

  
clean_up()
