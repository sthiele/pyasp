
from pyasp import Gringo4Clasp


solver = Gringo4Clasp(gringo_options='', clasp_options='2')
print(solver.run([], collapseTerms=True, collapseAtoms=False,
                 additionalProgramText='a:- not b. b:- not a.'))
