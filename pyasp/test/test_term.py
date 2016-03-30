"""
Unit testing about the term.py module, and mostly Term and TermSet classes.

"""

import unittest
from ..term import Term, TermSet


class TestTerm(unittest.TestCase):

    def test_compare_termset(self):
        self.assertEqual(TermSet.from_string('a b ck '),
                         TermSet.from_string('a b ck'))

    def test_termset_from_string(self):
        expected = TermSet((Term('a', []), Term('b', ['c', 'd']),
                            Term('p', ['q'])))
        self.assertEqual(TermSet.from_string('a b(c,d) p(q)'), expected)
