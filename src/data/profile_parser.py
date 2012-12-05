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

from bioasp.asp import *
from bioasp.misc import *
class Lexer:
	tokens = (
		'IDENT',
		'ASSIGN',
		'PLUS',
		'MINUS',
		'COMMENT',
	)

	# Tokens
	t_IDENT = r'[a-zA-Z][a-zA-Z0-9_:\-\[\]/]*'
	t_ASSIGN = r'='
	t_PLUS = r'\+'
	t_MINUS = r'-'
	t_COMMENT = r'%+'


	def __init__(self):
		import bioasp.ply.lex as lex
		self.lexer = lex.lex(object = self, optimize=1, lextab='profile_parser_lextab')

	# Ignored characters
	t_ignore = " \t"

	def t_newline(self, t):
		r'\n+'
		t.lexer.lineno += t.value.count("\n")
		
	def t_error(self, t):
		print "Illegal character '%s'" % t.value[0]
		t.lexer.skip(1)


class Parser:
	tokens = Lexer.tokens

	precedence = ( )
		#('left','PLUS','MINUS'),
		#('left','TIMES','DIVIDE'),
		#('right','UMINUS'),
		#)

	def __init__(self):
		self.experiment_name=""
		self.accu = TermSet()
		self.lexer = Lexer()
		import bioasp.ply.yacc as yacc
		#self.parser = yacc.yacc(module=self, tabmodule='calc_parsetab', debugfile="calc_parser.out")
		self.parser = yacc.yacc(module=self, optimize=1)

	def p_statement_expr(self, t):
		'''statement : IDENT ASSIGN value
					| IDENT 
					| COMMENT IDENT'''
		if len(t)==2 : 
			self.accu.add(Term('input', [ self.experiment_name, "gen(\""+t[1]+"\")"]))			
		if len(t)==3 : 
			self.experiment_name = "\""+t[2]+"\""
		if len(t)==4  :
			self.accu.add(Term('obs_vlabel', [self.experiment_name, "gen(\""+t[1]+"\")", t[3]]))
			


	def p_value(self, t):
		'''value : PLUS
		       | MINUS 
		       | IDENT '''
		if t[1] == '-' : t[0] = "-1"
		elif t[1] == '+' : t[0] = "1"
                elif t[1] == 'nc': t[0] = "0"

			
	def p_error(self, t):
		print "Syntax error at '%s'" % t

	def parse(self, line):
		self.parser.parse(line, lexer=self.lexer.lexer)
		return self.accu
