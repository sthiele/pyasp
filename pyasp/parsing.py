"""
Wrapper around the pyasp.ply submodule.

"""
import inspect

import pyasp.ply.lex as lex
import pyasp.ply.yacc as yacc
from pyasp.constant import OPTIMIZE
from pyasp.term import Term, TermSet


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
        self.lexer = lex.lex(object=self, optimize=OPTIMIZE, lextab='asp_py_lextab')

    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")

    def t_error(self, t):
        print("Illegal character " + str(t.value[0]))
        t.lexer.skip(1)


class Parser:
    start = 'answerset'

    def __init__(self, collapseTerms=True, collapseAtoms=False, callback=None):
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
        self.accu = TermSet()
        self.lexer = Lexer()
        self.tokens = self.lexer.tokens
        self.collapseTerms = collapseTerms
        self.collapseAtoms = collapseAtoms
        self.callback = callback

        if collapseAtoms and not collapseTerms:
            raise ValueError("if atoms are collapsed, functions must"
                             " also be collapsed!")

        self.parser = yacc.yacc(module=self, optimize=OPTIMIZE, tabmodule='asp_py_parsetab')

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
                t[0] = "%s(%s)" % ( t[1], ",".join(map(str, t[3])) )
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
        print("Syntax error at "+str(t))
        print (''.join(map(lambda x: "  %s:%s\n    %s" % (x[1], x[2], x[4][0]),
                           inspect.stack())))

    def parse(self, line):
        self.accu = TermSet()
        line = line.strip()

        if len(line) > 0:
            self.parser.parse(line, lexer=self.lexer.lexer)

        if self.callback:
            self.callback(self.accu)

        return self.accu


def filter_empty_str(l):
    return [x for x in l if x != '']
