"""
Unit test about the whole package.
Testing of the main API using examples found in Potassco documentation.

"""
import unittest

from .. import asp


def build_solver(constants={}, gringo_options='', clasp_options=''):
    """Return an asp.Gringo4Clasp object built using parameters values."""
    # build gringo options by adding the constants
    constants = ' -c '.join(str(k)+'='+str(v) for k, v in constants.items())
    if len(constants) > 0:
        constants = '-c ' + constants
    gringo_options = constants + ' ' + gringo_options

    return asp.Gringo4Clasp(clasp_options=clasp_options,
                            gringo_options=gringo_options)


class TestGringo4Clasp(unittest.TestCase):

    def unified_atoms(self, atoms):
        """Return the same atoms, but unified in order to be compared"""
        return frozenset(atoms)

    def unified_models(self, models):
        """Return the same models, but unified in order to be compared"""
        return set(self.unified_atoms(atoms) for atoms in models)


    def assert_models(self, program, expected_answers, constants={}):
        """Compare models found using given program with expected_answers."""
        solver = build_solver(constants, clasp_options='-n 0')
        models = solver.run([], additionalProgramText=program)
        self.assertEqual(len(models), len(expected_answers))
        self.assertEqual(
            self.unified_models(models),
            self.unified_models(expected_answers)
        )


    def test_solving_simple_neg(self):
        self.assert_models('a:- not b. b:- not a.', ('a', 'b'))

    def test_solving_double_neg(self):
        self.assert_models('a:- not not b. b:- not not a.', (('a', 'b'), ()))


    def test_constants(self):
        for const_value in ('1', 'a', '"t"', '23', '"text with space"'):
            try:
                with self.subTest(constant=const_value):
                    self.assert_models(
                        '#const v=1. p(v).',
                        (('p('+const_value+')',),),
                        constants={'v':const_value},
                    )
            except AttributeError:  # subTest is available only since Python 3.4
                self.assert_models(
                    '#const v=1. p(v).',
                    (('p('+const_value+')',),),
                    constants={'v':const_value},
                )


    def test_constant_with_space(self):
        """This test shows that space is not a valid constant letter"""
        constant = 'text with space'
        with self.assertRaises(Exception):
            solver = build_solver({'v':constant})
            solver.run([], additionalProgramText='#const v=1. p(v).')

    def test_version(self):
        self.assertIsInstance(asp.Gringo.version(), str)
        self.assertIsInstance(asp.Gringo.version_text(), str)
        self.assertIsInstance(asp.Gringo4.version(), str)
        self.assertIsInstance(asp.Gringo4.version_text(), str)
        self.assertIsInstance(asp.Clasp.version(), str)
        self.assertIsInstance(asp.Clasp.version_text(), str)
        self.assertIsInstance(asp.Gringo4Clasp.version(), tuple)
        self.assertIsInstance(asp.Gringo4Clasp.version_text(), tuple)
        self.assertIsInstance(asp.GringoClaspBase.version(), tuple)
        self.assertIsInstance(asp.GringoClaspBase.version_text(), tuple)
