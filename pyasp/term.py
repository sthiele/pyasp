"""
Definition of the Term class and associated objects

"""

import os
import re
import tempfile


class TermSet(set):
    """Set of Term instances, offering an API for atoms manipulation"""

    def __init__(self, terms=None):
        super(TermSet, self).__init__([] if terms is None else terms)
        self.score = None

    def filter(self, f):
        accu = TermSet()
        for e in self:
            if f(e): accu.add(e)
        return accu

    def to_list(self):
        ret = []
        for t in self: ret.append(t)
        return ret

    def to_file(self, fn=None):
        if fn:
            file = open(fn, 'w')
        else:
            fd, fn = tempfile.mkstemp('.lp')
            file = os.fdopen(fd, 'w')
        for t in self:
            file.write(str(t) + '.\n')
        file.close()
        return fn

    def exclude_rule(self):
        return ':- ' + ','.join(map(str, self)) + '.'

    @staticmethod
    def from_string(string):
        """Parse given string and return a TermSet instance
        containing found atoms."""
        from pyasp.parsing import Parser
        re.sub(r'\).\s*',") ", string)
        p = Parser(True, False)
        atoms = p.parse(string)
        return TermSet(atoms)


class String2TermSet(TermSet):
    """A TermSet with a constructor taking a string to be parsed into atoms.

    Deprecated, kept for retro-compatibility.
    Equivalent to TermSet.from_string usage.

    """
    def __new__(cls, string):
        return cls.from_string(string)


class Term:
    """an ASP term (used for atoms and for function terms)"""

    def __init__(self, predicate, arguments=None):
        self.predicate = predicate
        self.arguments = [] if arguments is None else arguments

    def nb_args(self):
        return len(self.arguments)

    def arg(self, n):
        return self.arguments[n]

    def args(self):
        return self.arguments

    def pred(self):
        return self.predicate

    def explode(self):
      return [self.pred()] + self.arguments

    def __repr__(self):
        if len(self.arguments) == 0:
            return "Term(%s)" % (repr(self.predicate),)
        else:
            return "Term(%s,[%s])" % (repr(self.predicate),
                                      ",".join(map(repr, self.arguments)))

    def __str__(self):
        if len(self.arguments) == 0:
            return self.predicate
        else:
            return self.predicate + "(" + ",".join(map(str, self.arguments)) + ")"

    def __hash__(self):
        return hash(tuple([self.predicate] + self.arguments))

    def __eq__(self, other):
        return self.predicate == other.predicate and self.arguments == other.arguments

    def sip(self, string):
        """sip = string in predicate ; Equivalent to string in self.predicate"""
        return string in self.predicate

    def p(self, string):
        """Equivalent to self.predicate == string"""
        return string == self.predicate
