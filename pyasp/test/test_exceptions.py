"""
Unit test about exception management.

"""
import unittest

from .. import asp


class TestGringo4Clasp(unittest.TestCase):

    def test_gringo_error(self):
        gringo = asp.Gringo4()
        with self.assertRaises(EnvironmentError):
            gringo.run([], additionalProgramText='..')

    def test_clasp_error(self):
        clasp = asp.Clasp()
        with self.assertRaises(EnvironmentError):
            clasp.run([], additionalProgramText='..')

    def test_clasp_error(self):
        solver = asp.Gringo4Clasp()
        answer = next(iter(solver.run([], additionalProgramText='a:- not b. b:- not a.')))
        self.assertEqual(answer.__class__, asp.TermSet)
        with self.assertRaises(TypeError):  # try to hash the TermSet
            set([answer])

