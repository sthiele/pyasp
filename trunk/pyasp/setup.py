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
import os
import sys
import platform
import distutils
import site
import urllib.request, urllib.parse, urllib.error
import sysconfig

from setuptools.command.install import install as _install

class install(_install):
    def get_binaries(self, path):
        BASE_URL = "http://www.cs.uni-potsdam.de/~sthiele/bioasp/downloads/bin/"
        architecture = platform.architecture()[0][:-3]
        if sys.platform == 'darwin':
            if architecture == '32':
                print("clasp/gringo3/gringo4 binaries are not yet available for Mac OS 32bits")
                exit()
                
            CLASP_URL = BASE_URL + "macos/clasp-3.1.1"
            GRINGO3_URL = BASE_URL + "macos/gringo-3.0.5"
            GRINGO4_URL = BASE_URL + "macos/gringo-4.4.0"
            
        else:                       
            CLASP_URL = BASE_URL + "linux-%s/clasp-3.1.1" % architecture
            GRINGO3_URL = BASE_URL + "linux-%s/gringo-3.0.5" % architecture
            GRINGO4_URL = BASE_URL + "linux-%s/gringo-4.4.0" % architecture
            
        
        urllib.request.urlretrieve(CLASP_URL, path + "/clasp")
        urllib.request.urlretrieve(GRINGO3_URL, path + "/gringo3")
        urllib.request.urlretrieve(GRINGO4_URL, path + "/gringo4")

        
    def run(self):
        _install.run(self)
        
        if '--user' in sys.argv[-1] :
            userdir = site.getusersitepackages()+"/pyasp/bin" 
            self.get_binaries(userdir)
            cmd="chmod +x "+userdir+"/*"
            print(cmd)
            os.system(cmd)
        else :
            py_version = "%s.%s" % (sys.version_info[0], sys.version_info[1])
            paths = []
            paths.append(sys.prefix+"/lib/python%s/dist-packages/pyasp/bin" % py_version)
            paths.append(sys.prefix+"/lib/python%s/site-packages/pyasp/bin" % py_version)
            paths.append(sys.prefix+"/local/lib/python%s/dist-packages/pyasp/bin" % py_version)
            paths.append(sys.prefix+"/local/lib/python%s/site-packages/pyasp/bin" % py_version)
            paths.append("/Library/Python/%s/site-packages/pyasp/bin" % py_version)

            cmd = None
            for path in paths:
                if os.path.exists(path):
                    self.get_binaries(path)
                    cmd = "chmod +x "+path+"/*"
                    break
                
            if cmd:
                print(cmd)
                os.system(cmd)
            else:
                print("pyasp binaries path not found. You need to download and put in place the binaries for gringo3, gringo4 and clasp in order to start using pyasp.")
                
setup(
    cmdclass={'install': install},
    name = 'pyasp',
    version = '1.4.1',
    url='http://pypi.python.org/pypi/pyasp/',
    license='GPLv3+',   
    description='A convenience wrapper for the ASP tools gringo, gringo4 and clasp.',
    long_description=open('README').read(),
    author='Sven Thiele',
    author_email='sthiele78@gmail.com', 
    
    package_dir = { 'pyasp' : 'src'},
    package_data = {
        'pyasp' : ['bin/*.txt']
    },
    packages = [
        'pyasp', 
        'pyasp.ply'
    ]
)
