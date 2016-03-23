"""
Definition of many constants.

"""

import sys
import pkg_resources  # packaging facilies
from functools import partial

from .info import __pkg_name__, __pkg_version__


# Directories (relative to the package level)
REL_DIR_SOURCES   = ''  # sources are inside the package
REL_DIR_BIN       = 'bin/'

# Packaging: access to the file, inside the package,
#   independently of the installation directory
access_packaged_file = partial(pkg_resources.resource_filename, __pkg_name__)

# Access to the binary files
BIN_GRINGO3  = access_packaged_file(REL_DIR_BIN + 'gringo3')
BIN_GRINGO4 = access_packaged_file(REL_DIR_BIN + 'gringo4')
BIN_CLASP   = access_packaged_file(REL_DIR_BIN + 'clasp')


# Optimization definition
# use 0 for debugging
# use 1 for speed (run python once without -O or -OO and then with -O or -OO)
OPTIMIZE = 1
if sys.flags.optimize in ['-O', '-OO'] and OPTIMIZE == 0:
  raise RuntimeError("need to use OPTIMIZE = 1 and run python without -O"
                     " or -OO once to make parsers work")
