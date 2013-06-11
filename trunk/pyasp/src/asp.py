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
  def nb_args(self): return len(self.arguments)
  def arg(self,n): return self.arguments[n]
  def args(self): return self.arguments
  def pred(self): return self.predicate
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

  def __init__(self,collapseTerms=True,collapseAtoms=False):
    '''
    collapseTerms: function terms in predicate arguments are collapsed into strings
    collapseAtoms: atoms (predicate plus terms) are collapsed into strings
                   requires that collapseTerms is True

    example: a(b,c(d))
    collapseTerms=True,  collapseAtoms=False: result = Term('a', ['b', 'c(d)'])
    collapseTerms=True,  collapseAtoms=True:  result = 'a(b,c(d))'
    collapseTerms=False, collapseAtoms=False: result = Term('a', ['b', Term('c', ['d'])])
    collapseTerms=False, collapseAtoms=True: invalid arguments
    '''
    global optimize
    self.accu = TermSet()
    self.lexer = Lexer()
    self.tokens = self.lexer.tokens
    self.collapseTerms = collapseTerms
    self.collapseAtoms = collapseAtoms
    if collapseAtoms and not collapseTerms:
      raise "if atoms are collapsed, functions must also be collapsed!"
    #self.parser = yacc.yacc(module=self, optimize=optimize, tabmodule='asp_py_parsetab', debugfile="asp_py_parser.dbg")
    self.parser = yacc.yacc(module=self, optimize=optimize, tabmodule='asp_py_parsetab')

  def p_answerset(self, t):
    '''answerset : atom SPACE answerset
                 | atom'''
    self.accu.add(t[1])

  def p_atom(self, t):
    '''atom : IDENT LP terms RP
            | IDENT
            | MIDENT LP terms RP
            | MIDENT'''
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
    '''terms : term COMMA terms
             | term'''
    if len(t) == 2:
      t[0] = [t[1]]
    else:
      t[0] = [t[1]] + t[3]

  def p_term(self, t):
    '''term : IDENT LP terms RP
            | STRING
            | IDENT
            | NUM'''
    if self.collapseTerms:
      if len(t) == 2:
        #print "term %s" % (t[1],)
        t[0] = t[1]
      else:
        t[0] = t[1] + "(" + ",".join(t[3]) + ")"
    else:
      if len(t) == 2:
        #print "term %s" % (t[1],)
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
    #print "parsing '%s'" % (line,)
    self.accu = TermSet()
    line = line.strip()
    if len(line) > 0:
      self.parser.parse(line, lexer=self.lexer.lexer)
    return self.accu

def filter_empty_str(l):
    return filter(lambda x: x != '',l)

#class String2Term(Atom):
#  def __init__(self,strg):
#    p = Parser(collapseTerms=True,collapseAtoms=False)
#    ret = p.parse(strg)
#    if len(ret) != 1:
#      raise RuntimeError("expected just 1 atom here, " +
#        "got '%s' for parsing '%s'" % (repr(ret),strg))
#    ret = ret.pop() # take arbitrary element from set
#    #print "ret %s" % (ret,)
#    Atom.__init__(self,ret.pred(),ret.args())

#
class String2TermSet(TermSet):
    def __init__(self,s):
        re.sub(r'\).\s*',") ",s)
        p = Parser(collapseTerms=True,collapseAtoms=False)
        atoms = p.parse(s)
        #print "string2termset: atoms '%s'" % (atoms,)
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
        


#space_re = re.compile('[ \t\n]+')

class GringoClaspBase:
    def __init__(self, clasp_bin = root + '/bin/clasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = ''):
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

    def getClaspStderr(self):
        return self.clasp_stderr

    def getGringoStderr(self):
        return self.gringo_stderr

    def parseStats(self, prefix, stderrcontent):
        #debug('parseStats('+prefix+') for' + stderrcontent)
        valuelines = re.split(r'[\r\n]', stderrcontent)
        stats = {}
        current = ''
        #debug(valuepairs)
        for line in valuelines:
          #debug('"'+line+'"')
          valuepairs = re.findall(r'[^(]+?:\s*[0-9]+[^ )]*', line)
          #debug(valuepairs)
          first = True
          for pair in valuepairs:
            (fkey, value) = re.split(r':', pair)
            key = fkey.strip()
            value = value.strip()
            # get toplevel key
            if first and fkey[0] != ' ':
              current = key
            first = False
            # store values
            intelligent_key = key
            if current != key:
              intelligent_key = current + '/' + key
            stats[prefix+intelligent_key] = value
        #debug(stats)
        return stats

    def parseClaspStats(self):
        return self.parseStats('clasp/', self.getClaspStderr())

    def parseGringoStats(self):
        return self.parseStats('gringo/', self.getGringoStderr())

class GringoClasp(GringoClaspBase):
    def __init__(self, *args, **keywords):
        GringoClaspBase.__init__(self, *args, **keywords)

    def run(self, programs, nmodels = 1, collapseTerms=True, collapseAtoms=True, additionalProgramText=None):
        assert(programs.__class__ == list)
        # options need to be filtered to work with Popen
        # since gringo/clasp doesn't like whitespace before numbers parameter
        try:
            additionalPrograms = []
            if additionalProgramText != None:
              (fd, fn) = tempfile.mkstemp('.lp')
              file = os.fdopen(fd,'w')
              file.write(str(additionalProgramText))
              file.close()
              additionalPrograms.append(fn)
              
            addoptions = []
            if self.gringo_options != None and self.gringo_options != '':
              addoptions = self.gringo_options.split()
              assert(addoptions.__class__ == list)

            commandline = filter_empty_str([self.gringo_bin] + addoptions + programs + additionalPrograms)
            #debug(commandline)
            self._gringo = subprocess.Popen(
              commandline, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #print 'gringo', self.gringo_options.split() + programs
        except OSError, e:
            if e.errno == 2:
                raise Exception('Grounder \'%s\' not found' % self.gringo_bin)
            else: raise
        try:
            addoptions = []
            if self.clasp_options != None and self.clasp_options != '':
              addoptions = self.clasp_options.split()
              assert(addoptions.__class__ == list)

            self._clasp = subprocess.Popen(
                filter_empty_str([self.clasp_bin] + addoptions + [str(nmodels)]),
                stdin=self._gringo.stdout,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE)
            #print self.clasp_options.split() + [str(nmodels)]
        except OSError, e:
            if e.errno == 2:
                raise Exception('Solver \'%s\' not found' % self.clasp_bin)
            else: raise
        accu = []
        l = self._clasp.stdout.readline()
        #f = self._clasp.stderr.readline()
        #print f,'f'
        # take only last set of cautious consequences
        parser = Parser(collapseTerms=collapseTerms,collapseAtoms=collapseAtoms)
        if '--enum-mode cautious' in self.clasp_options:
            lastline = None
            while l[:-1] != '':
                if l.startswith('Cautious consequences'):
                    lastline = self._clasp.stdout.readline()
                l = self._clasp.stdout.readline()
            if lastline != None:
              #print >>sys.stderr, "parsing cautious consequences %s" % (lastline,)
              accu.append(parser.parse(lastline))
        elif '--enum-mode brave' in self.clasp_options:
            lastline = None
            while l[:-1] != '':
                if l.startswith('Brave consequences'):
                    lastline = self._clasp.stdout.readline()
                l = self._clasp.stdout.readline()
            if lastline != None:
              #print >>sys.stderr, "parsing brave consequences %s" % (lastline,)
              accu.append(parser.parse(lastline))
        else:
            while l != '':
                #debug(l)
                if l.startswith('Answer'):
                    l = self._clasp.stdout.readline()
                    #print >>sys.stderr, "parsing %s" % (l,)
                    accu.append(parser.parse(l))
                    #print >>sys.stderr, "accu = %s" % (accu,)
                if l.startswith('Models'):
                    # statistics start here
                    break
                l = self._clasp.stdout.readline()

        # error handling (we can only do it here as gringo may not be finished until here)

        # first clasp
        (claspout,self.clasp_stderr) = self._clasp.communicate()
        # unfortunately the clasp stats appear on stdout and not on stderr, so we have to pfusch here
        self.clasp_stderr = claspout + self.clasp_stderr
        if self._clasp.returncode not in self.clasp_noerror_retval:
          # make sure gringo is dead
          self._gringo.kill()
          self._gringo.wait()
          self.gringo_stderr = self._gringo.stderr.read()
          raise Exception("got error %d from clasp: '%s' gringo: '%s'" % \
            (self._clasp.returncode, self.clasp_stderr, self.gringo_stderr))

        # if clasp terminated successfully, then gringo also terminated
        (gringoout,self.gringo_stderr) = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
          raise Exception("got error %d from gringo: '%s'" % \
            (self._gringo.returncode, self.gringo_stderr))

        self._clasp = None
        self._gringo = None
        return accu

    def run_print(self, programs, printer, nmodels = 1, additionalProgramText=None):
        assert(programs.__class__ == list)
        collapseTerms=False
        collapseAtoms=False
        # options need to be filtered to work with Popen
        # since gringo/clasp doesn't like whitespace before numbers parameter
        try:
            additionalPrograms = []
            if additionalProgramText != None:
              (fd, fn) = tempfile.mkstemp('.lp')
              file = os.fdopen(fd,'w')
              file.write(str(additionalProgramText))
              file.close()
              additionalPrograms.append(fn)
              
            addoptions = []
            if self.gringo_options != None and self.gringo_options != '':
              addoptions = self.gringo_options.split()
              assert(addoptions.__class__ == list)

            commandline = filter_empty_str([self.gringo_bin] + addoptions + programs + additionalPrograms)
            #debug(commandline)
            self._gringo = subprocess.Popen(
              commandline, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #print 'gringo', self.gringo_options.split() + programs
        except OSError, e:
            if e.errno == 2:
                raise Exception('Grounder \'%s\' not found' % self.gringo_bin)
            else: raise
        try:
            addoptions = []
            if self.clasp_options != None and self.clasp_options != '':
              addoptions = self.clasp_options.split()
              assert(addoptions.__class__ == list)

            self._clasp = subprocess.Popen(
                filter_empty_str([self.clasp_bin] + addoptions + [str(nmodels)]),
                stdin=self._gringo.stdout,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE)
            #print self.clasp_options.split() + [str(nmodels)]
        except OSError, e:
            if e.errno == 2:
                raise Exception('Solver \'%s\' not found' % self.clasp_bin)
            else: raise
        accu = []
        l = self._clasp.stdout.readline()
        #f = self._clasp.stderr.readline()
        #print f,'f'
        parser = Parser(collapseTerms=collapseTerms,collapseAtoms=collapseAtoms)
        count=1
        while l != '':

            if l.startswith('Answer'):
	        count += 1
                l = self._clasp.stdout.readline()
                #print count,":",parser.parse(l)
                printer.write(count, parser.parse(l))
                
            if l.startswith('Models'):
                # statistics start here
                break
            l = self._clasp.stdout.readline()

        # error handling (we can only do it here as gringo may not be finished until here)

        # first clasp
        (claspout,self.clasp_stderr) = self._clasp.communicate()
        # unfortunately the clasp stats appear on stdout and not on stderr, so we have to pfusch here
        self.clasp_stderr = claspout + self.clasp_stderr
        if self._clasp.returncode not in self.clasp_noerror_retval:
          # make sure gringo is dead
          self._gringo.kill()
          self._gringo.wait()
          self.gringo_stderr = self._gringo.stderr.read()
          raise Exception("got error %d from clasp: '%s' gringo: '%s'" % \
            (self._clasp.returncode, self.clasp_stderr, self.gringo_stderr))

        # if clasp terminated successfully, then gringo also terminated
        (gringoout,self.gringo_stderr) = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
          raise Exception("got error %d from gringo: '%s'" % \
            (self._gringo.returncode, self.gringo_stderr))

        self._clasp = None
        self._gringo = None
        return count

class GringoHClasp(GringoClaspBase):
    def __init__(self, *args, **keywords):
        GringoClaspBase.__init__(self, *args, **keywords)
        self.clasp_bin = root + '/bin/hclasp'
        self.clasp_noerror_retval = set([0,10,20])

    def run(self, programs, nmodels = 1, collapseTerms=True, collapseAtoms=True, additionalProgramText=None):
        assert(programs.__class__ == list)
        # options need to be filtered to work with Popen
        # since gringo/clasp doesn't like whitespace before numbers parameter
        try:
            additionalPrograms = []
            if additionalProgramText != None:
              (fd, fn) = tempfile.mkstemp('.lp')
              file = os.fdopen(fd,'w')
              file.write(str(additionalProgramText))
              file.close()
              additionalPrograms.append(fn)
              
            addoptions = []
            if self.gringo_options != None and self.gringo_options != '':
              addoptions = self.gringo_options.split()
              assert(addoptions.__class__ == list)

            commandline = filter_empty_str([self.gringo_bin] + addoptions + programs + additionalPrograms)
            #debug(commandline)
            self._gringo = subprocess.Popen(
              commandline, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            #print 'gringo', self.gringo_options.split() + programs
        except OSError, e:
            if e.errno == 2:
                raise Exception('Grounder \'%s\' not found' % self.gringo_bin)
            else: raise
        try:
            addoptions = []
            if self.clasp_options != None and self.clasp_options != '':
              addoptions = self.clasp_options.split()
              assert(addoptions.__class__ == list)

            self._clasp = subprocess.Popen(
                filter_empty_str([self.clasp_bin] + addoptions + [str(nmodels)]),
                stdin=self._gringo.stdout,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE)
            #print self.clasp_options.split() + [str(nmodels)]
        except OSError, e:
            if e.errno == 2:
                raise Exception('Solver \'%s\' not found' % self.clasp_bin)
            else: raise
        accu = []
        
        
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
       
        # error handling (we can only do it here as gringo may not be finished until here)

        # first clasp
        (claspout,self.clasp_stderr) = self._clasp.communicate()
        # unfortunately the clasp stats appear on stdout and not on stderr, so we have to pfusch here
        self.clasp_stderr = claspout + self.clasp_stderr
        if self._clasp.returncode==20: return None
        if self._clasp.returncode not in self.clasp_noerror_retval:
          # make sure gringo is dead
          self._gringo.kill()
          self._gringo.wait()
          self.gringo_stderr = self._gringo.stderr.read()
          raise Exception("got error %d from clasp: '%s' gringo: '%s'" % \
            (self._clasp.returncode, self.clasp_stderr, self.gringo_stderr))

        # if clasp terminated successfully, then gringo also terminated
        (gringoout,self.gringo_stderr) = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
          raise Exception("got error %d from gringo: '%s'" % \
            (self._gringo.returncode, self.gringo_stderr))

        self._clasp = None
        self._gringo = None
        
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
        model = parser.parse(lastmodelline)
        return (optimal_weights, [model])
        
      
class GringoClaspD(GringoClasp):
    def __init__(self, *args, **keywords):
        GringoClasp.__init__(self, *args, **keywords)
        self.clasp_bin = root + '/bin/claspD'
        self.clasp_noerror_retval = set([0])



class GringoClaspOpt(GringoClaspBase):
    def __init__(self, *args, **keywords):
        GringoClaspBase.__init__(self, *args, **keywords)

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
            filter_empty_str([self.clasp_bin, '-n', '0'] +
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
        
        
 
class GringoClaspDistr(GringoClaspBase):
    def __init__(self, clasp_bin = root + '/bin/clasp', clasp_options = '',
                       gringo_bin = root + '/bin/gringo', gringo_options = ''):
        self.clasp_bin = clasp_bin
        self.gringo_bin = gringo_bin
        self.clasp_options = clasp_options
        self.gringo_options = gringo_options
        self._clasp = None
        self._gringo = None

    def run(self, programs):
        # options need to be filtered to work with Popen
        # since gringo/clasp doesn't like whitespace before numbers parameter
        self._gringo = subprocess.Popen(
            filter_empty_str([self.gringo_bin] + self.gringo_options.split() + programs),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self._clasp = subprocess.Popen(
            filter_empty_str([self.clasp_bin] + self.clasp_options.split() + ['--number=0']),
            stdin=self._gringo.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        l = self._clasp.stdout.readline()
        table={}
        countofAS=0
        while l!= '':
            if l[:6] == 'Answer':
                 l = self._clasp.stdout.readline()
                 model = TermSet()
                 atoms=l.split()
                 for atom in atoms :
                    model=model.union(String2TermSet(atom))
                 countofAS += 1
                 #print countofAS
                 for atom  in model :
                    if table.has_key(atom) : table[atom] += 1
                    else : table[atom]=1
            l = self._clasp.stdout.readline()
        self._clasp = None
        self._gringo = None
        opts = l[13:].split()
        keys = table.keys()
        for key in keys : 
           val = table[key]
           table[key] = (100.0/countofAS)*val
        return table       

                
class AbortedException(BaseException):
    def __init__(self):
        self.value = "Solving has been aborted"

class ThreadedGringoClasp(GringoClasp, threading.Thread):
    def __init__(self, *args, **keywords):
        GringoClasp.__init__(self, *args, **keywords)
        threading.Thread.__init__(self)
        self._params = None
        self.answers = None
        self.job_end_cb = None
  
    def set_params(self, *args, **keywords):
        """
        Set the solver parameters. The arguments are the same as GringoClasp.run.
        """
        self._params = (args, keywords)

    def run(self):
        """
        Not to be called directly. Call start() instead. If successful,
        answers will contain the answer sets.
        """
        if self._params is None:
            raise Exception("Parameters not set!")
        self.answers = GringoClasp.run(self, *self._params[0], **self._params[1])
        if '__call__' in dir(self.job_end_cb):
            self.job_end_cb()

    def kill(self):
        """
        Force the solver to abort. This causes the Thread to return with an AbortedException.
        """
        if self._clasp == None or self._gringo == None:
            return
        os.system("kill %s %s" % (self._gringo.pid, self._clasp.pid))
        if '__call__' in dir(self.job_end_cb):
            self.job_end_cb()

class ThreadedGringoClasp2(ThreadedGringoClasp, threading.Thread):
    """
    Computes new answers by excluding already found ones.
    Could use a better name. ;)
    """
    def run(self):
        self._killed = False
        if self._params is None:
            raise Exception("Parameters not set!")
        has_answers = True
        self.answers = []
        excluded_answers = tempfile.NamedTemporaryFile('w')
        self._params[0][0].append(excluded_answers.name)

        # remove 'nmodels' parameter
        self._params = list(self._params)
        self._params[0] = list(self._params[0])
        nmodels = 1
        for param in self._params[0]:
            if isinstance(param, int):
                self._params[0].remove(param)
                nmodels = param
        nmodels_key = None
        for key in self._params[1]:
            if key == 'nmodels':
                nmodels_key = key
        if nmodels_key:
            nmodels = self._params[1].pop(nmodels_key)

        # set 'nmodels' to 1
        self._params[0].append(1)

        # compute answers
        i = 0
        has_answers = True
        while has_answers and not self._killed:
            i += 1
            tmp_answers = GringoClasp.run(self, *self._params[0], **self._params[1])
            has_answers = len(tmp_answers) > 0
            if not has_answers or (nmodels != 0 and i > nmodels):
                break
            self.answers += tmp_answers
            excluded_answers.write(tmp_answers[0].exclude_rule())
            excluded_answers.flush()
        excluded_answers.close()
        if '__call__' in dir(self.job_end_cb):
            self.job_end_cb()
    
    def kill(self):
        """
        Force the solver to abort. This causes the Thread to return with an AbortedException.
        """
        self._killed = True
        if self._clasp == None or self._gringo == None:
            return
        os.system("kill %s %s" % (self._gringo.pid, self._clasp.pid))
        if '__call__' in dir(self.job_end_cb):
            self.job_end_cb()
