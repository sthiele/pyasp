# Copyright (c) 2012, Sven Thiele <sthiele78@gmail.com>
#
# This file is part of pyasp.
#
# pyasp is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyasp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyasp.  If not, see <http://www.gnu.org/licenses/>.import random

# -*- coding: utf-8 -*-
#from distutils.core import setup
from setuptools import setup
from setuptools.command.install import install as _install
import os
import sys
import site
import logging
import platform
import sysconfig
import subprocess

# version specific imports
if sys.version_info >= (3,):
    from urllib import request
    urlretrieve = request.urlretrieve
else:  # python 2 compatibility
    # TODO: using requests module (http://docs.python-requests.org/en/latest/index.html)
    #  is maybe a good idea
    from urllib import urlretrieve

# import pyasp, get access to binary files
from pyasp import constant
from pyasp import info


BINARIES_BASE_URL = 'http://www.cs.uni-potsdam.de/~sthiele/bioasp/downloads/bin/{}/{}'
BINARIES_NAME = {
    # binary remote name: binary local name
    'clasp-3.1.3': 'clasp',
    'gringo-3.0.5': 'gringo3',
    'gringo-4.5.3': 'gringo4',
}

BASE_URL_PLATFORM_SPECIFIC_SUBPATHS = {
    # (platform, architecture) : path to binaries from BINARIES_BASE_URL
    ('linux', '32'): 'linux-32',
    ('linux', '64'): 'linux-64',
    ('win32', '32'): 'windows-32',
    ('win32', '64'): 'windows-64',
    ('cygwin', '32'): 'windows-32',
    ('cygwin', '64'): 'windows-64',
    ('darwin', '64'): 'macos',  # darwin 32 bits is not supported
}


def system_info():
    """Return (platform name, architecture), or (None, None) if not found"""
    # get architecture and platform information
    arch = platform.architecture()[0][:-3]
    platform_name = None
    for platform_name, _ in BASE_URL_PLATFORM_SPECIFIC_SUBPATHS:
        if sys.platform.startswith(platform_name):
            return platform_name, arch
    return None, None


def binaries_urls(platform_name, arch):
    """Return tuple of binaries name and URL, based on given sys informations.

    If detected system is not supported, an empty iterable is returned.
    Architecture and platforms supported depends of distant binary repository.

    """
    try:
        subpath = BASE_URL_PLATFORM_SPECIFIC_SUBPATHS[platform_name, arch]
    except KeyError:
        logging.getLogger().error(
            'clasp/gringo3/gringo4 binaries are not available for'
            ' platform ' + platform_name + ' under architecture '
            + arch + 'bits.')
        return tuple()  # empty iterable
    # no error: build the tuple of binaries paths
    return tuple((local_name, BINARIES_BASE_URL.format(subpath, remote_name))
                 for remote_name, local_name in BINARIES_NAME.items())


def binaries_directory():
    """Return the installation directory, or None"""
    if '--user' in sys.argv:
        paths = (site.getusersitepackages(),)
    else:
        py_version = '%s.%s' % (sys.version_info[0], sys.version_info[1])
        paths = (s % (py_version) for s in (
            sys.prefix + '/lib/python%s/dist-packages/',
            sys.prefix + '/lib/python%s/site-packages/',
            sys.prefix + '/local/lib/python%s/dist-packages/',
            sys.prefix + '/local/lib/python%s/site-packages/',
            '/Library/Python/%s/site-packages/',
        ))

    # yield the first valid path
    for path in paths:
        # add the package and bin subdir and, if exists, return it
        path = os.path.join(path, '%s/%s' % (info.__pkg_name__, constant.REL_DIR_BIN))
        if os.path.exists(path):
            return path
    logging.getLogger().error(
        'pyasp binaries path not found. You need to download and'
        ' put in place the binaries for gringo3, gringo4 and clasp'
        ' in order to start using pyasp.'
        ' You can find binaries from ' + BINARIES_BASE_URL +
        ' or https://sourceforge.net/projects/potassco/files/,'
        ' or compile them yourself.'
    )
    exit()
    return None


def post_install():
    """Get the binaries online, and give them the execution permission"""
    bin_dir = binaries_directory()
    for binary_name, binary_url in binaries_urls(*system_info()):
        logging.getLogger().info('RETRIEVE: ' + os.path.join(bin_dir, binary_name)
                                 + ' from ' + binary_url)
        bin_path = os.path.join(bin_dir, binary_name)
        urlretrieve(binary_url, bin_path)
        logging.getLogger().info('CHMOD COMMAND: '
                                 + ' '.join(['chmod', '+x', bin_path]))
        subprocess.call(['chmod', '+x', bin_path])


class install(_install):

    def run(self):
        """Call superclass run method, then downloads the binaries"""
        _install.run(self)
        self.execute(post_install, args=[], msg=post_install.__doc__)


setup(
    cmdclass={'install': install},
    name = info.__pkg_name__,
    version = info.__pkg_version__,
    url='http://pypi.python.org/pypi/pyasp/',
    license='GPLv3+',
    description='A convenience wrapper for the ASP tools gringo, gringo4 and clasp.',
    long_description=open('README.rst').read(),
    author='Sven Thiele',
    author_email='sthiele78@gmail.com',
    zip_safe=False,  # not zippable because of the binary retrieving

    package_dir = { 'pyasp' : 'pyasp'},
    package_data = {
        'pyasp' : ['bin/*.txt']
    },
    packages = [
        'pyasp',
        'pyasp.ply'
    ],
)
