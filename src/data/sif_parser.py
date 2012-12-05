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
		'PLUS',
		'MINUS',
	)

	# Tokens

	t_IDENT = r'[a-zA-Z][a-zA-Z0-9_:\-\[\]/]*'

	t_PLUS = r'1'
	t_MINUS = r'-1'


	def __init__(self):
		import bioasp.ply.lex as lex
		self.lexer = lex.lex(object = self, optimize=1, lextab='sif_parser_lextab')

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

	def __init__(self):
		# dictionary of names
		self.aux_node_counter=0
		#self.names = { }
		self.accu = TermSet()
		self.args = []
		self.lexer = Lexer()
		import bioasp.ply.yacc as yacc
		#self.parser = yacc.yacc(module=self, tabmodule='calc_parsetab', debugfile="calc_parser.out")
		self.parser = yacc.yacc(module=self,optimize=1)

	def p_statement_expr(self, t):
		'''statement : node_expression PLUS node_expression 
					| node_expression MINUS node_expression'''
		if len(t)<3 : 
			self.accu.add(Term('input', [t[1]]))
			print 'input', t[1]
		else :
			#print t[1], t[2], t[3]
			self.accu.add(Term('edge', ["gen(\""+t[1]+"\")","gen(\""+t[3]+"\")"]))
			self.accu.add(Term('obs_elabel', ["gen(\""+t[1]+"\")","gen(\""+t[3]+"\")",t[2]]))
			#print Term('obs_elabel', ["gen(\""+t[1]+"\")","gen(\""+t[3]+"\")",t[2]])


	def p_node_expression(self, t):
		'''node_expression : IDENT'''
		if len(t)<3 : 
				t[0]=t[1]
				#print t[1]
				self.accu.add(Term('vertex', ["gen(\""+t[1]+"\")"]))
		else : t[0] = "unknown"

			
	def p_error(self, t):
		print "Syntax error at '%s'" % t

	def parse(self, line):
		self.parser.parse(line, lexer=self.lexer.lexer)
		return self.accu
