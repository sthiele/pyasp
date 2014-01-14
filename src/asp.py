# Copyright (c) 2012, Sven Thiele <sthiele78@gmail.com>
#
# This file is part of pyasp.
#
# pyasp is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyasp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyasp.  If not, see <http://www.gnu.org/licenses/>.import random

# -*- coding: utf-8 -*-
import re
import tempfile
import os
import subprocess
import threading
from pyasp.misc import *
import pyasp.ply.lex as lex
import pyasp.ply.yacc as yacc
import errno
import json

root = __file__.rsplit('/', 1)[0]

global optimize
# use 0 for debugging
# use 1 for speed (run python once without -O or -OO and then with -O or -OO)
optimize=1
import sys
if sys.flags.optimize in ['-O', '-OO'] and optimize == 0:
  raise RuntimeError("need to use optimize=1 and run python without -O or -OO once to make parsers work")

def debug(s):
  import inspect
  print >>sys.stderr, "DBG %03d: %s" % \
    (inspect.currentframe().f_back.f_lineno,s)

# mutable and not hashable!
class TermSet(set):
    def __init__(self, terms=[]):
        super(TermSet, self).__init__(terms)
        self.score = None
        
    def filter(self,f):
        accu = TermSet()
        for e in self:
            if f(e): accu.add(e)
        return accu
    def to_list(self):
        ret = []
        for t in self: ret.append(t)
        return ret
    def to_file(self,fn=None):
        if fn:
            file = open(fn,'w')
        else:
            fd, fn = tempfile.mkstemp('.lp')
            file = os.fdopen(fd,'w')
        for t in self:
            #print str(t)
            file.write(str(t) + '.\n')
        file.close()
        return fn
    def exclude_rule(self):
        return ':- ' + ','.join(map(str,self)) + '.'
    def __hash__(self):
        raise "TermSet is mutable and not hashable!"

# immutable and hashable
class Term:
    """an ASP term (used for atoms and for function terms)"""
    def __init__(self,predicate,arguments=[]):
        self.predicate = predicate
        self.arguments = arguments

    def nb_args(self): 
        return len(self.arguments)

    def arg(self,n): 
        return self.arguments[n]
    
    def args(self): 
        return self.arguments

    def pred(self): 
        return self.predicate

    def explode(self):
      return [ self.pred() ] + self.arguments

    def __repr__(self):
        if len(self.arguments) == 0:
            return "Term(%s)" % (repr(self.predicate),)
        else:
            return "Term(%s,[%s])" % (repr(self.predicate),",".join(map(repr,self.arguments)))
    
    def __str__(self):
        if len(self.arguments) == 0:
            return self.predicate
        else:
            return self.predicate + "(" + ",".join(map(str,self.arguments)) + ")"
    
    def __hash__(self):
        return tuple([self.predicate] + self.arguments).__hash__()
    
    def __eq__(self,other):
        return self.predicate == other.predicate and self.arguments == other.arguments
    
    def sip(self,s):
        """sip = string in predicate"""
        return (s in self.predicate)
    
    def p(self,s):
        return (s == self.predicate)

class Lexer:
    tokens = (
        'STRING',
        'IDENT',
        'MIDENT',
        'NUM',
        'LP',
        'RP',
        'COMMA',
        'SPACE',
    )

    # Tokens
    t_STRING = r'"((\\")|[^"])*"'
    t_IDENT = r'[a-zA-Z_][a-zA-Z0-9_]*'
    t_MIDENT = r'-[a-zA-Z_][a-zA-Z0-9_]*'
    t_NUM = r'-?[0-9]+'
    t_LP = r'\('
    t_RP = r'\)'
    t_COMMA = r','
    t_SPACE = r'[ \t\.]+'

    def __init__(self):
        global optimize
        self.lexer = lex.lex(object=self,optimize=optimize, lextab='asp_py_lextab')

    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")

    def t_error(self, t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skip(1)

class Parser:
    start = 'answerset'

    def __init__(self,collapseTerms=True,collapseAtoms=False, callback=None):
        """
        collapseTerms: function terms in predicate arguments are collapsed into strings
        collapseAtoms: atoms (predicate plus terms) are collapsed into strings
                       requires that collapseTerms is True

        example: a(b,c(d))
        collapseTerms=True,  collapseAtoms=False: result = Term('a', ['b', 'c(d)'])
        collapseTerms=True,  collapseAtoms=True:  result = 'a(b,c(d))'
        collapseTerms=False, collapseAtoms=False: result = Term('a', ['b', Term('c', ['d'])])
        collapseTerms=False, collapseAtoms=True: invalid arguments
        """
        global optimize
        self.accu = TermSet()
        self.lexer = Lexer()
        self.tokens = self.lexer.tokens
        self.collapseTerms = collapseTerms
        self.collapseAtoms = collapseAtoms
        self.callback = callback

        if collapseAtoms and not collapseTerms:
            raise "if atoms are collapsed, functions must also be collapsed!"
        
        self.parser = yacc.yacc(module=self, optimize=optimize, tabmodule='asp_py_parsetab')

    def p_answerset(self, t):
        """answerset : atom SPACE answerset
                   | atom
        """
        self.accu.add(t[1])

    def p_atom(self, t):
        """atom : IDENT LP terms RP
                | IDENT
                | MIDENT LP terms RP
                | MIDENT
        """
        if self.collapseAtoms:
            if len(t) == 2:
                t[0] = str(t[1])
            elif len(t) == 5:
                t[0] = "%s(%s)" % ( t[1], ",".join(map(str,t[3])) )
        else:
            if len(t) == 2:
                t[0] = Term(t[1])
            elif len(t) == 5:
                t[0] = Term(t[1], t[3])

    def p_terms(self, t):
        """terms : term COMMA terms
                 | term
        """
        if len(t) == 2:
            t[0] = [t[1]]
        else:
            t[0] = [t[1]] + t[3]

    def p_term(self, t):
        """term : IDENT LP terms RP
                | STRING
                | IDENT
                | NUM
        """
        if self.collapseTerms:
            if len(t) == 2:
                t[0] = t[1]
            else:
                t[0] = t[1] + "(" + ",".join(t[3]) + ")"
        else:
            if len(t) == 2:
                if re.match(r'-?[0-9]+', t[1]) != None:
                    t[0] = int(t[1])
                else:
                    t[0] = t[1]
            else:
                t[0] = Term(t[1], t[3])

    def p_error(self, t):
        print "Syntax error at '%s'" % t
        import inspect
        print ''.join(map(lambda x: "  %s:%s\n    %s" % (x[1], x[2], x[4][0]),inspect.stack()))

    def parse(self, line):
        self.accu = TermSet()
        line = line.strip()
    
        if len(line) > 0:
            self.parser.parse(line, lexer=self.lexer.lexer)

        if self.callback:
            self.callback(self.accu)
        
        return self.accu

def filter_empty_str(l):
    return filter(lambda x: x != '',l)

class String2TermSet(TermSet):
    def __init__(self,s):
        re.sub(r'\).\s*',") ",s)
        p = Parser(True, False)
        atoms = p.parse(s)
        TermSet.__init__(self,atoms)
    
def exclude_sol(sols,fn=None):
    if fn:
        file = open(fn,'w')
    else:
        fd, fn = tempfile.mkstemp('.lp')
        file = os.fdopen(fd,'w')
    for s in sols:
        file.write(s.exclude_rule() + '\n')
        
    file.close()
    return fn
        
class GringoClaspBase(object):
    def __init__(self, clasp_bin = root + '/bin/clasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = False):
        self.clasp_bin = clasp_bin
        self.gringo_bin = gringo_bin
        self.clasp_options = clasp_options
        self.gringo_options = gringo_options
        self._clasp = None
        self._gringo = None
        self.clasp_stderr = None
        self.gringo_stderr = None
        self.clasp_noerror_retval = set([10, 20, 30])
        self.gringo_noerror_retval = set([0])
        
        if clasp_options.find('--opt-all') != -1:
            #workaround for backwards compatibility
            optimization = True
        
        self.optimization = optimization
        
    def __parse_witnesses__(self, parser, witnesses, key):
        accu = []
        for answer in witnesses:
            ts = parser.parse(" ".join(answer[key]))
            if answer.has_key('Opt'):
                ts.score = answer['Opt']
                
            accu.append(ts)
            
        return accu
        
    def __get_witnesses_key__(self,result):
        key = 'Value'
        if result['Stats'].has_key('Brave'):
            key = 'Brave'
        elif result['Stats'].has_key('Cautious'):
            key = 'Cautious'
            
        return key
        
    def __ground__(self, programs, additionalProgramText):
        try:
            additionalPrograms = []
            if additionalProgramText != None:
                (fd, fn) = tempfile.mkstemp('.lp')
                file = os.fdopen(fd,'w')
                file.write(str(additionalProgramText))
                file.close()
                additionalPrograms.append(fn)
              
            addoptions = []
            if self.gringo_options:
                addoptions = self.gringo_options.split()

            commandline = filter_empty_str([self.gringo_bin] + addoptions + programs + additionalPrograms)
            self._gringo = subprocess.Popen(commandline, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            
        except OSError, e:
            if e.errno == 2:
                raise Exception('Grounder \'%s\' not found' % self.gringo_bin)
            else: 
                raise e
                
        (grounding,self.gringo_stderr) = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
            raise Exception("got error %d from gringo: '%s'" % \
                (self._gringo.returncode, self.gringo_stderr))
                
        return grounding
                        
    def __solve__(self, nmodels, grounding, opts = [], json=True):
        try:
            addoptions = []
            if self.clasp_options:
                addoptions = self.clasp_options.split()

            if json:
                opts = opts + ['--outf=2']
                
            self._clasp = subprocess.Popen(
                filter_empty_str([self.clasp_bin] + addoptions + [str(nmodels)] + opts),
                stdin=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE)
                
        except OSError, e:
            if e.errno == 2:
                raise Exception('Solver \'%s\' not found' % self.clasp_bin)
            else: 
                raise e
                
        (solving,self.clasp_stderr) = self._clasp.communicate(grounding)
        self.clasp_stderr = solving + self.clasp_stderr
        
        if self._clasp.returncode not in self.clasp_noerror_retval:
            raise Exception("got error %d from clasp: '%s' gringo: '%s'" % \
                (self._clasp.returncode, self.clasp_stderr, self.gringo_stderr))
                    
        self._clasp = None
        self._gringo = None
        
        return solving
                
    def run(self, programs, nmodels = 1, collapseTerms=True, collapseAtoms=True, additionalProgramText=None, callback=None):
        grounding = self.__ground__(programs, additionalProgramText)
        if self.optimization:
            solving = self.__solve__(0, grounding)
        else:
            solving = self.__solve__(nmodels, grounding)
        
        parser = Parser(collapseTerms,collapseAtoms,callback)
        res = json.loads(solving)
        key = self.__get_witnesses_key__(res)
        
        if res['Result'] == "SATISFIABLE":
            if key == 'Value':
                witnesses = res['Witnesses']
            else:
                witnesses = [res['Witnesses'][-1]]

            accu = self.__parse_witnesses__(parser, witnesses, key)            
                
        elif res['Result'] == "OPTIMUM FOUND":
            optimum = res['Witnesses'][-1]
            if nmodels == 1:
                witnesses = [optimum]  
            elif (nmodels != 1) and self.clasp_options.find('--opt-all') == -1:
                optarg = '--opt-all=' + ','.join(map(str, optimum['Opt']))            
                solving = self.__solve__(nmodels, grounding, [optarg])
                optall = json.loads(solving)
            
                if key == 'Value':
                    witnesses = optall['Witnesses']
                else:
                    witnesses = [optall['Witnesses'][-1]]
            else:
                if key == 'Value':
                    witnesses = res['Witnesses']
                else:
                    witnesses = [optimum]
            
            accu = self.__parse_witnesses__(parser, witnesses, key)
            
        else:
            accu = []

        return accu

class GringoClasp(GringoClaspBase):
    pass

class GringoClaspOpt(GringoClaspBase):
    def __init__(self, clasp_bin = root + '/bin/clasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = True):
                       
        super(GringoClaspOpt, self).__init__(clasp_bin, clasp_options, gringo_bin, gringo_options, optimization)
    
class GringoHClasp(GringoClaspBase):
    def __init__(self, clasp_bin = root + '/bin/hclasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = False):
                   
        super(GringoHClasp, self).__init__(clasp_bin, clasp_options, gringo_bin, gringo_options, optimization)
            
    def __parse_witnesses__(self, parser, witnesses, key):
        accu = []
        for answer in witnesses:
            atoms = filter(lambda atom: not atom.startswith('_'), answer[key])
            ts = parser.parse(" ".join(atoms))
            if answer.has_key('Opt'):
                ts.score = answer['Opt']
            accu.append(ts)
        return accu
                   
        
class GringoHClaspOpt(GringoHClasp):
    def __init__(self, clasp_bin = root + '/bin/hclasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = True):
                   
        super(GringoHClaspOpt, self).__init__(clasp_bin, clasp_options, gringo_bin, gringo_options, optimization)

      
class GringoClaspD(GringoClasp):
    def __init__(self, clasp_bin = root + '/bin/claspD', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = False):
                   
        super(GringoClaspD, self).__init__(clasp_bin, clasp_options, gringo_bin, gringo_options, optimization)
        
    def __get_witnesses_key__(self,result):
        if result['Models'].has_key('Brave'):
            key = 'Brave'
        elif result['Models'].has_key('Cautious'):
            key = 'Cautious'
        else:
            key = 'Value'
            
        return key

class GringoUnClasp(GringoClaspBase):
    def __init__(self, clasp_bin = root + '/bin/unclasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = '', 
                       optimization = True):
                   
        super(GringoUnClasp, self).__init__(clasp_bin, clasp_options, gringo_bin, gringo_options, optimization)
        

    def run(self, programs, nmodels = 1, collapseTerms=True, collapseAtoms=True, additionalProgramText=None, callback=None):
        grounding = self.__ground__(programs, additionalProgramText)
        solving = self.__solve__(nmodels, grounding, json=False)

        lines = solving.split('\n')
        parser = Parser(collapseTerms,collapseAtoms,callback)
        if "SATISFIABLE" in lines:
            optimum = parser.parse(lines[-10])
            optimum.score = tuple(map(int, lines[-3][13:].split()))
            accu = [optimum]
        else:
            accu = []

        return accu

class GringoUnClaspOpt(GringoClaspBase):
    def __init__(self, *args, **keywords):
        GringoClaspBase.__init__(self, *args, **keywords)
        self.clasp_bin = root + '/bin/unclasp'

    def run(self, programs, nmodels = 1, collapseTerms=True, collapseAtoms=True, additionalProgramText=None):
        '''
        allows for the following modes of optimization:
        nmodels=0: get all optimal models
        nmodels=1: get first optimal model (-n 0)
        nmodels>1: get first K optimal models
                   (this enumerates all optimal models using --opt-all)
        (in both cases non-optimal models are discarded)

        returns tuple:
        (quality-tuple, list of optimal models)

        the default setting returns the first optimal model
        '''
        # options need to be filtered to work with Popen
        # since gringo/clasp doesn't like whitespace before numbers parameter
        additionalPrograms = []
        if additionalProgramText != None:
          (fd, fn) = tempfile.mkstemp('.lp')
          file = os.fdopen(fd,'w')
          file.write(str(additionalProgramText))
          file.close()
          additionalPrograms.append(fn)
          
        self._gringo = subprocess.Popen(
            filter_empty_str([self.gringo_bin] +
              self.gringo_options.split() + programs + additionalPrograms),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        # get ground program
        (groundout,self.gringo_stderr) = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
          raise RuntimeError("gringo terminated with %d and error %s" %
            (self._gringo.returncode, self.gringo_stderr))
        #
        # now we first do clasp -n 0 to get weights of the optimal model
        # (for getting all models, we can later initialize the solver with this
        # weight and use --opt-all)
        self._clasp = subprocess.Popen(
            filter_empty_str([self.clasp_bin,'--quiet=0'] +
              self.clasp_options.split()),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        # send ground program to clasp and signal end of input
        self._clasp.stdin.write(groundout)
        self._clasp.stdin.close()
        # get resulting answer sets, ignore everything to last (optimal) model
        lastweightline = None
        lastmodelline = None
        l = self._clasp.stdout.readline()
        while l != '' and l[:13] != 'OPTIMUM FOUND':
            #debug(l.strip())
            if l.startswith('Answer'):
                # get model
                lastmodelline = self._clasp.stdout.readline()
                #debug(lastmodelline.strip())
                # get weights
                lastweightline = self._clasp.stdout.readline()
                #debug(lastweightline.strip())
                if lastweightline[:12] != 'Optimization':
                  # not an optimization problem!
                  raise RuntimeError("expected 'Optimization' line, got '%s' " %
                    lastweightline.strip() + "(not an optimization problem?)")
            if l.startswith('UNSATISFIABLE'):
                return None
            l = self._clasp.stdout.readline()

        # clasp execution error handling
        self._clasp.wait()
        self.clasp_stderr = self._clasp.stderr.read()
        if self._clasp.returncode not in self.clasp_noerror_retval:
          raise Exception("got error %d from clasp: '%s' gringo: '%s'" % \
            (self._clasp.returncode, self.clasp_stderr, self.gringo_stderr))

        # received content error handling
        if lastweightline == None or lastmodelline == None:
          claspstdout = self._clasp.stdout.read()
          self.clasp_stderr = self._clasp.stderr.read()
          raise RuntimeError("did not find optimum! stdout='%s' stderr='%s'" %
            (claspstdout, self.clasp_stderr))
        # interpret weight of optimal model
        optimal_weights = lastweightline[13:].split()
        optimal_weights = tuple(map(int,optimal_weights))
        #print >>sys.stderr, "optimal weights %s" % (optimal_weights,)
        #
        # we know the optimal weight -> get optimal model(s) next
        #
        parser = Parser(collapseTerms=collapseTerms,collapseAtoms=collapseAtoms)
        # if we only need 1 model, parse and return it
        if nmodels == 1:
          #print >>sys.stderr, "parsing %s" % (lastmodelline,)
          model = parser.parse(lastmodelline)
          return (optimal_weights, [model])
        # we need more models -> call solver again with preset opt value
        # --opt-val=<optimal_weights> and --opt-all and -n <nmodels>
        # -> we get a maximum of <nmodels> models and they are all optimal
        optarg = '--opt-all=' + ','.join(map(str,optimal_weights))
        args = filter_empty_str([self.clasp_bin, '-n', str(nmodels),
          optarg] + self.clasp_options.split())
        #print >>sys.stderr, "calling clasp again with %s" % (args,)
        self._clasp = subprocess.Popen(args, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # send ground program to clasp and signal end of input
        self._clasp.stdin.write(groundout)
        self._clasp.stdin.close()
        # parse and return all resulting answer sets
        accu = []
        l = self._clasp.stdout.readline()
        while l != '':
            #print >>sys.stderr, "DBG got line '%s'" % (l,)
            if l.startswith('Answer'):
                # get model
                modelline = self._clasp.stdout.readline()
                # get weights
                weightline = self._clasp.stdout.readline()
                assert(weightline == lastweightline)
                #print >>sys.stderr, "parsing modelline %s" % (modelline,)
                model = parser.parse(modelline)
                # build result
                accu.append(model)
            if l.startswith('Models'):
                # statistics start here
                break
            l = self._clasp.stdout.readline()

        # clasp execution error handling
        self._clasp.wait()
        # unfortunately the clasp stats appear on stdout and not on stderr, so we have to pfusch here
        self.clasp_stderr = self._clasp.stdout.read()
        self.clasp_stderr += self._clasp.stderr.read()
        if self._clasp.returncode not in self.clasp_noerror_retval:
          raise Exception("got error %d from clasp: '%s'" % (self._clasp.returncode, self.clasp_stderr))

        self._clasp = None
        self._gringo = None
        return (optimal_weights, accu)
        