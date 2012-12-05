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
		'ARROW',
		'PLUS',
		'MINUS',
		'AND',
		#'COMPL',
		#'OR',
		#'LP',
		#'RP',
		#'COM',
	)

	# Tokens

	t_IDENT = r'[a-zA-Z][a-zA-Z0-9_:\-\[\]/]*'
	t_ARROW = r'->'
	t_PLUS = r'\+'
	t_MINUS = r'-'
	t_AND = r'\?'
	t_COMPL = r'&'
	t_OR = r'\|'
	t_LP = r'\('
	t_RP = r'\)'
	t_COM = r','

	def __init__(self):
		import bioasp.ply.lex as lex
		self.lexer = lex.lex(object = self, optimize=1, lextab='graph_parser_lextab')

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
		self.aux_node_counter=0
		self.accu = TermSet()
		self.fun_flag = False
		self.args = []
		self.commas = []
		self.lexer = Lexer()
		import bioasp.ply.yacc as yacc
		#self.parser = yacc.yacc(module=self, debug=False, tabmodule='calc_parsetab', debugfile="calc_parser.out")
		self.parser = yacc.yacc(module=self, optimize=1, tabmodule='graph_parser_parsetab',)

	def p_statement_expr(self, t):
		'''statement : node_expression ARROW node_expression value
					| node_expression '''
		if len(t)<3 : 
			self.accu.add(Term('input', [t[1]]))
		else :
			self.accu.add(Term('edge', [t[1],t[3]]))
			if t[4]!="0" : self.accu.add(Term('obs_elabel', [t[1],t[3],t[4]]))
			#self.names[t[0]] = t[2]

	def p_node_expression(self, t):
		'''node_expression : IDENT '''

		if len(t)<3 : 
			
			if self.fun_flag==True :
				t[0]=t[1]
				self.fun_flag=False
			else :
				t[0] = "gen(\""+t[1]+"\")"
				self.accu.add(Term('vertex', ["gen(\""+t[1]+"\")"]))
		elif t[1] == 'opposite_sign' :
			self.aux_node_counter+=1
			aux_node = "aux_node_"+str(self.aux_node_counter)
			t[0] = "opposite_sign"+"("+aux_node+")"
			self.accu.add(Term('edge', [t[3],t[1]+"("+aux_node+")"]))
			self.accu.add(Term('obs_elabel', [t[3],t[1]+"("+aux_node+")","1"]))
		elif t[1] == 'strong_inhibitor' :
			t[0] = "strong_inhibitor"+"("+t[3]+")"
			for i in self.commas : 
				self.accu.add(Term('edge', [i,t[1]+"("+t[3]+")"]))
				self.accu.add(Term('obs_elabel', [i,t[1]+"("+t[3]+")","1"]))
			self.commas= []
		else : t[0] = "unknown"

	def p_value(self, t):
		'''value : PLUS
		       | MINUS
		       | AND '''
		if t[1] == '-' : t[0] = "-1"
		elif t[1] == '+' : t[0] = "1"
		elif t[1] == '?' : t[0] = "0"

	#def p_commalist(self, t):
		#'''commalist : IDENT COM IDENT'''
					  ##| commalist COM IDENT'''
		##t[0] = t[1]+" , "+t[3]
		#t[0] = "gen(\""+t[1]+"\"),gen(\""+t[3]+"\")"
		#self.accu.add(Term('vertex', ["gen(\""+t[1]+"\")"]))
		#self.accu.add(Term('vertex', ["gen(\""+t[3]+"\")"]))
		#self.commas.append("gen(\""+t[1]+"\")")
		#self.commas.append("gen(\""+t[3]+"\")")

	#def p_orexpression2(self, t):
		#'orlist2 : orlist'
		#self.aux_node_counter+=1
		#aux_node = "aux_node_"+str(self.aux_node_counter)
		#self.accu.add(Term('vertex', ["or("+aux_node+")"]))
		#for i in self.args : 
			#self.accu.add(Term('edge', [i,"or("+aux_node+")"]))
			#self.accu.add(Term('obs_elabel', [i,"or("+aux_node+")","1"]))
		#self.args= []
		#t[0]="or("+aux_node+")"
		
	#def p_andexpression2(self, t):
		#'andlist2 : andlist'
		#self.aux_node_counter+=1
		#aux_node = "aux_node_"+str(self.aux_node_counter)
		#self.accu.add(Term('vertex', ["and("+aux_node+")"]))
		#for i in self.args : 
			#self.accu.add(Term('edge', [i,"and("+aux_node+")"]))
			#self.accu.add(Term('obs_elabel', [i,"and("+aux_node+")","1"]))
		#self.args= []
		#t[0]="and("+aux_node+")"

	#def p_orexpression(self, t):
		#'''orlist : IDENT
				#| orlist OR IDENT'''
		#self.fun_flag=True
		#if len(t)<3 : 
			#self.accu.add(Term('vertex', ["gen(\""+t[1]+"\")"]))
			#self.args.append("gen(\""+t[1]+"\")")
			#t[0] = "gen(\""+t[1]+"\")"
		#elif t[2] == '|' : 
			#self.accu.add(Term('vertex', ["gen(\""+t[3]+"\")"]))
			#self.args.append("gen(\""+t[3]+"\")")
			#t[0] = "Unknown"

	#def p_andexpression(self, t):
		#'''andlist : IDENT
					  #| andlist COMPL IDENT'''
		#self.fun_flag=True
		#if len(t)<3 : 
			#self.accu.add(Term('vertex', ["gen(\""+t[1]+"\")"]))
			#self.args.append("gen(\""+t[1]+"\")")
			#t[0] = "gen(\""+t[1]+"\")"
		#elif t[2] == '&' : 
			#self.accu.add(Term('vertex', ["gen(\""+t[3]+"\")"]))
			#self.args.append("gen(\""+t[3]+"\")")
			#t[0] = "Unknown"
			
	def p_error(self, t):
		print "Syntax error at '%s'" % t

	def parse(self, line):
		self.parser.parse(line, lexer=self.lexer.lexer)
		return self.accu
