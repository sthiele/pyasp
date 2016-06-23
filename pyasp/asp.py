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
from __future__ import print_function
import errno
import json
import os
import subprocess
import tempfile
import threading
import re

from pyasp.term import Term, TermSet
from pyasp.constant import BIN_GRINGO3, BIN_GRINGO4, BIN_CLASP
from pyasp.parsing import filter_empty_str, Parser


class GringoClaspBase(object):
    def __init__(self, clasp_bin=BIN_CLASP, clasp_options='',
                 gringo_bin=BIN_GRINGO3, gringo_options='',
                 optimization=False):
        self.clasp_bin = clasp_bin
        self.gringo_bin = gringo_bin
        self.clasp_options = clasp_options
        self.gringo_options = gringo_options
        self._clasp = None
        self._gringo = None
        self.clasp_stderr = None
        self.gringo_stderr = None
        self.clasp_noerror_retval = set([10, 20, 30])
        self.gringo_noerror_retval = set([0])

        if clasp_options.find('--opt-all') != -1:
            #workaround for backwards compatibility
            optimization = True

        self.optimization = optimization

    def __parse_witnesses__(self, parser, witnesses):
        accu = []
        for answer in witnesses:
            ts = parser.parse(" ".join(answer['Value']))
            if 'Costs' in answer:
                ts.score = answer['Costs']

            accu.append(ts)

        return accu

    def __get_witnesses_key__(self, result):
        key = 'Value'
        if 'Brave' in result['Models']:
            key = 'Brave'
        elif 'Cautious' in result['Models']:
            key = 'Cautious'

        return key

    @staticmethod
    def version_text(gringo_bin=BIN_GRINGO3,
                     clasp_bin=BIN_CLASP):
        """Return the version text"""
        if gringo_bin:
            try:
                commandline = filter_empty_str([gringo_bin] + ['--version'])
                gringo = subprocess.Popen(commandline, stderr=subprocess.PIPE,
                                          stdout=subprocess.PIPE)
            except OSError as e:
                if e.errno == 2:
                    raise OSError('Grounder \'%s\' not found' % gringo_bin)
                else:
                    raise e
            gringo_version, _ = gringo.communicate()
            gringo_version = gringo_version.decode('utf-8')
        else:
            gringo_version = None

        if clasp_bin:
            try:
                clasp = subprocess.Popen(
                    filter_empty_str([clasp_bin] + ['--version']),
                    stdin = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    stdout = subprocess.PIPE
                )

            except OSError as e:
                if e.errno == 2:
                    raise Exception('Solver \'%s\' not found' % clasp_bin)
                else:
                    raise e
            clasp_version, _ = clasp.communicate()
            clasp_version = clasp_version.decode('utf-8')
        else:
            clasp_version = None
        return gringo_version, clasp_version

    @staticmethod
    def version(gringo_bin=BIN_GRINGO3, clasp_bin=BIN_CLASP):
        """Return the version number as string"""
        gringo, clasp = GringoClaspBase.version_text(gringo_bin, clasp_bin)
        return (
            REGEX_VERSION_NUMBER.search(gringo).group() if gringo else None,
            REGEX_VERSION_NUMBER.search(clasp).group() if clasp else None,
        )

    def __ground__(self, programs, additionalProgramText):
        try:
            additionalPrograms = []
            if additionalProgramText != None:
                (fd, fn) = tempfile.mkstemp('.lp')
                file = os.fdopen(fd, 'w')
                file.write(str(additionalProgramText))
                file.close()
                additionalPrograms.append(fn)

            addoptions = []
            if self.gringo_options:
                addoptions = re.findall(r'(?:[^\s,"]|"(?:\\.|[^"])*")+', self.gringo_options)

            commandline = filter_empty_str([self.gringo_bin] + addoptions + programs + additionalPrograms)
            self._gringo = subprocess.Popen(commandline, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

        except OSError as e:
            if e.errno == 2:
                raise OSError('Grounder \'%s\' not found' % self.gringo_bin)
            else:
                raise e

        grounding, self.gringo_stderr = self._gringo.communicate()
        if self._gringo.returncode not in self.gringo_noerror_retval:
            raise EnvironmentError("got error %d from gringo:\n%s" % \
                (self._gringo.returncode, self.gringo_stderr.decode()))

        return grounding

    def __solve__(self, grounding, opts=[], json=True):
        try:
            #opts = opts + ['--stats']
            addoptions = []
            if self.clasp_options:
                addoptions = self.clasp_options.split()

            if json:
                opts = opts + ['--outf=2']
            self._clasp = subprocess.Popen(
                filter_empty_str([self.clasp_bin] + addoptions  + opts),
                stdin = subprocess.PIPE,
                stderr = subprocess.PIPE,
                stdout = subprocess.PIPE)

        except OSError as e:
            if e.errno == 2:
                raise Exception('Solver \'%s\' not found' % self.clasp_bin)
            else:
                raise e

        solving, self.clasp_stderr = self._clasp.communicate(grounding)
        self.clasp_stderr = solving + self.clasp_stderr

        if self._clasp.returncode not in self.clasp_noerror_retval:
            error = "got error %d from clasp:\n%s" % (self._clasp.returncode,
                                                      self.clasp_stderr.decode())
            try:
                if self.gringo_stderr:
                    error += "\n\nFrom gringo:\n%s" % self.gringo_stderr
            except AttributeError:  # no gringo error
                pass
            raise EnvironmentError(error)

        self._clasp = None
        self._gringo = None

        return solving

    def run(self, programs, collapseTerms=True, collapseAtoms=True, additionalProgramText=None, callback=None):
        grounding = self.__ground__(programs, additionalProgramText)

        solving = self.__solve__(grounding)

        parser = Parser(collapseTerms, collapseAtoms, callback)
        res = json.loads(solving.decode())
        key = self.__get_witnesses_key__(res)

        if res['Result'] == "SATISFIABLE":
            if key == 'Value':
                witnesses = res['Call'][0]['Witnesses']
            else:
                witnesses = [res['Call'][0]['Witnesses'][-1]]

            accu = self.__parse_witnesses__(parser, witnesses)

        elif res['Result'] == "OPTIMUM FOUND":
            if key == 'Value':
                numopts = res['Models']['Optimal']
                if numopts == 0 :
                    print('WARNING: OPTIMUM FOUND but zero optimals')
                    witnesses = []
                else :
                    offset = len(res['Call'][0]['Witnesses'])-numopts
                    witnesses = res['Call'][0]['Witnesses'][offset:]

            else:
                witnesses = [res['Call'][0]['Witnesses'][-1]]

            accu = self.__parse_witnesses__(parser, witnesses)

        else:
            accu = []

        return accu


class GringoClasp(GringoClaspBase):
    pass


class Gringo(GringoClaspBase):
    def __init__(self, gringo_bin = BIN_GRINGO3, gringo_options = ''):
        self.gringo_bin = gringo_bin
        self.gringo_options = gringo_options
        self._gringo = None
        self.gringo_stderr = None
        self.gringo_noerror_retval = set([0])

    def run(self, programs, collapseTerms=True, collapseAtoms=True, additionalProgramText=None, callback=None):
        return self.__ground__(programs, additionalProgramText)

    @staticmethod
    def version(gringo_bin=BIN_GRINGO3):
        gringo, _ = GringoClaspBase.version(gringo_bin, clasp_bin=None)
        return gringo

    @staticmethod
    def version_text(gringo_bin=BIN_GRINGO3):
        gringo, _ = GringoClaspBase.version_text(gringo_bin, clasp_bin=None)
        return gringo

class Gringo4Clasp(GringoClaspBase):
    def __init__(self, clasp_bin = BIN_CLASP, clasp_options = '',
                       gringo_bin = BIN_GRINGO4, gringo_options = '',
                       optimization = False):
        self.clasp_bin = clasp_bin
        self.gringo_bin = gringo_bin
        self.clasp_options = clasp_options
        self.gringo_options = gringo_options
        self._clasp = None
        self._gringo = None
        self.clasp_stderr = None
        self.gringo_stderr = None
        self.clasp_noerror_retval = set([10, 20, 30])
        self.gringo_noerror_retval = set([0])

        if clasp_options.find('--opt-all') != -1:
            #workaround for backwards compatibility
            optimization = True

        self.optimization = optimization


class Gringo4(Gringo4Clasp):
    def __init__(self, gringo_bin=BIN_GRINGO4, gringo_options=''):
        self.gringo_bin = gringo_bin
        self.gringo_options = gringo_options
        self._gringo = None
        self.gringo_stderr = None
        self.gringo_noerror_retval = set([0])

    def run(self, programs, collapseTerms=True, collapseAtoms=True, additionalProgramText=None, callback=None):
        return self.__ground__(programs, additionalProgramText)

    @staticmethod
    def version(gringo_bin=BIN_GRINGO4):
        gringo, _ = GringoClaspBase.version(gringo_bin, clasp_bin=None)
        return gringo

    @staticmethod
    def version_text(gringo_bin=BIN_GRINGO4):
        gringo, _ = GringoClaspBase.version_text(gringo_bin, clasp_bin=None)
        return gringo


class Clasp(GringoClaspBase):
    def __init__(self, clasp_bin=BIN_CLASP, clasp_options='',
                       optimization=False):
        self.clasp_bin = clasp_bin
        self.clasp_options = clasp_options
        self._clasp = None
        self.clasp_stderr = None
        self.clasp_noerror_retval = set([10, 20, 30])

        if clasp_options.find('--opt-all') != -1:
            #workaround for backwards compatibility
            optimization = True

        self.optimization = optimization

    def run(self, grounded, collapseTerms=True, collapseAtoms=True,
            additionalProgramText=None, callback=None):
        solving = self.__solve__(grounded)

        parser = Parser(collapseTerms, collapseAtoms, callback)
        res = json.loads(solving.decode())
        key = self.__get_witnesses_key__(res)

        if res['Result'] == "SATISFIABLE":
            if key == 'Value':
                witnesses = res['Call'][0]['Witnesses']
            else:
                witnesses = [res['Call'][0]['Witnesses'][-1]]

            accu = self.__parse_witnesses__(parser, witnesses)

        elif res['Result'] == "OPTIMUM FOUND":
            if key == 'Value':
                numopts = res['Models']['Optimal']
                if numopts == 0 :
                    print('WARNING: OPTIMUM FOUND but zero optimals')
                    witnesses = []
                else :
                    offset = len(res['Call'][0]['Witnesses'])-numopts
                    witnesses = res['Call'][0]['Witnesses'][offset:]

            else:
                witnesses = [res['Call'][0]['Witnesses'][-1]]

            accu = self.__parse_witnesses__(parser, witnesses)

        else:
            accu = []

        return accu


    @staticmethod
    def version(clasp_bin=BIN_CLASP):
        _, clasp = GringoClaspBase.version(gringo_bin=None, clasp_bin=clasp_bin)
        return clasp

    @staticmethod
    def version_text(clasp_bin=BIN_CLASP):
        _, clasp = GringoClaspBase.version_text(gringo_bin=None, clasp_bin=clasp_bin)
        return clasp
