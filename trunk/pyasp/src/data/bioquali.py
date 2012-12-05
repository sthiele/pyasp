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
import re
from bioasp.asp import *
from bioasp.misc import *

import bioasp.ply.lex as lex
import bioasp.ply.yacc as yacc

import graph_parser
import sif_parser
import profile_parser


def int_of_sign(s):
    if s == '+': return '1'
    elif s == '-': return '-1'
    elif s == '0': return '0'
    elif s == 'nc': return '0'
    else: 
        print s
        assert(False)

def sign_of_int(s):
    if s == '1': return '+'
    elif s == '-1': return '-'
    elif s == '0': return '0'
    else: 
        print s
        assert(False)



GENE_ID = '[-a-zA-Z0-9_:&\(\)/]+'

#tokens = ( 'COMPLEX', 'INFL', 'SIGN', 'GENE',)
#t_GENE = r'[a-zA-Z_][a-zA-Z0-9_]*' 
#operator = ['&','|']

def readGraph(filename):
    p = graph_parser.Parser()
    """
    input: string, name of a file containing a Bioquali-like graph description
    output: asp.TermSet, with atoms matching the contents of the input file
    
    Parses a Bioquali-like graph description, and returns
    a TermSet object.
    Written using original Bioquali
    """
	
    accu = TermSet()
    file = open(filename,'r')
    s = file.readline()
    while s!="":
        try:
            accu = p.parse(s)
        except EOFError:
            break
        s = file.readline()

    return accu

def readSIFGraph(filename):
    p = sif_parser.Parser()
    """
    input: string, name of a file containing a Bioquali-like graph description
    output: asp.TermSet, with atoms matching the contents of the input file
    
    Parses a Bioquali-like graph description, and returns
    a TermSet object.
    Written using original Bioquali
    """
	
    accu = TermSet()
    file = open(filename,'r')
    s = file.readline()
    while s!="":
        try:
            accu = p.parse(s)
        except EOFError:
            break
        s = file.readline()

    return accu


def saveGraph(termset,filename):
    file = open(filename, 'w')
    elabels = termset.filter(lambda t: t.sip('obs_elabel'))
    edges = termset.filter(lambda t: t.sip('edge'))
    inputs = termset.filter(lambda t: t.sip('input'))
    for t in elabels:
        [ _, j, i, s ] = t.explode()
        file.write(unquote(j) + ' -> ' + unquote(i) + ' ' + sign_of_int(s) + '\n')
        edges.discard(Term('edge',[j, i]))
    for t in edges:
        [_, j, i] = t.explode()
        file.write(unquote(j) + ' -> ' + unquote(i) + ' X\n')
    for t in inputs:
        [ _, i] = t.explode()
        file.write(unquote(i) + '\n')
    file.close()


def readProfile_new(filename):
    p = profile_parser.Parser()
    """
    input: string, name of a file containing a Bioquali-like graph description
    output: asp.TermSet, with atoms matching the contents of the input file
    
    Parses a Bioquali-like graph description, and returns
    a TermSet object.
    Written using original Bioquali
    """
	
    accu = TermSet()
    file = open(filename,'r')
    s = file.readline()
    while s!="":
        if s!="\n":
	    try:
		accu = p.parse(s)
	    except EOFError:
		break
	s = file.readline()

    return accu
    
def readProfile(filename):
    #COMPLEX_ID = '[-a-zA-Z0-9_:\(\)/&]+'
    GENE_ID = '[-a-zA-Z0-9_:\(\)/]+'
    SIGN = '[-+0n]c*'
    file = open(filename,'r')
    val_re = '(?P<genid>'+GENE_ID+')( )*=( )*(?P<sign>'+SIGN+')'
    val = re.compile(val_re)
    name_re = '(-\*-)( )*(?P<name>\S+)( )*(-\*-)'
    name = re.compile(name_re)
    obs_name = None
    line_number = 1
    line = file.readline()[:-1]
    # the first line may contain the name
    nm = name.match(line)
    if nm:
        name = nm.group('name')
        line = file.readline()[:-1]
    else:
        name = filename
    name = quote(name)
    accu = TermSet()
    accu.add(Term('exp',[name]))
    # ordinary line
    while line:
        vm = val.match(line)
        if vm:
            vertex = quote(vm.group('genid'))
            accu.add(Term('obs_vlabel',[name, "gen("+vertex+")", int_of_sign(vm.group('sign'))]))
        else:
            print 'Syntax error line:', line_number,'  '+line 
        line = file.readline()[:-1]
        line_number+=1
    return accu

def readProfile_with_nc(filename):
    #COMPLEX_ID = '[-a-zA-Z0-9_:\(\)/&]+'
    GENE_ID = '[-a-zA-Z0-9_:\(\)/]+'
    file = open(filename,'r')
    val_re = '(?P<genid>'+GENE_ID+')( )*=( )*(?P<sign>[-+0n]c*)( )*(?P<input>\(input\))*'
    val = re.compile(val_re)
    name_re = '(-\*-)( )*(?P<name>\S+)( )*(-\*-)'
    name = re.compile(name_re)
    obs_name = None
    line_number = 1
    line = file.readline()[:-1]
    # the first line may contain the name
    nm = name.match(line)
    if nm:
        name = nm.group('name')
        line = file.readline()[:-1]
    else:
        name = filename
    name = quote(name)
    accu = TermSet()
    accu.add(Term('exp',[name]))
    # ordinary line
    while line:
        vm = val.match(line)
        if vm:
            vertex = quote(vm.group('genid'))

            accu.add(Term('obs_vlabel',[name, "gen("+vertex+")", int_of_sign(vm.group('sign'))]))
            if vm.group('input')!=None : 
               accu.add(Term('input',[name, "gen("+vertex+")"]))
        else:
            print 'Syntax error line:', line_number,'  '+line 
        line = file.readline()[:-1]
        line_number+=1
    return accu
    
def saveProfile(termset,condition,filename):
    file = open(filename, 'w')
    file.write('-*- ' + condition + ' -*-\n')
    for t in termset:
        if t.sip('obs_vlabel'):
            [_, c, i, s] = t.explode()
            if condition == unquote(c):
                file.write(unquote(i) + ' = ' + sign_of_int(s) + '\n')
    file.close()

                
