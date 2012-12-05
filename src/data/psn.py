# Copyright (c) 2012, Sven Thiele <sthiele78@gmail.com>
#
# This file is part of BioASP.
#
# BioASP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioASP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioASP.  If not, see <http://www.gnu.org/licenses/>.import random

# -*- coding: utf-8 -*-
import re
from bioasp.asp import *
from bioasp.misc import *

import bioasp.ply.lex as lex
import bioasp.ply.yacc as yacc

import sif2bm_parser
import csv  
import math

def readSIF(filename):
    p = sif2bm_parser.Parser()
	
    accu = TermSet()
    file = open(filename,'r')
    s = file.readline().replace('\r\n','\n')
    while s!="":
        if s!="\n":
	    try:
		accu = p.parse(s)
	    except EOFError:
		break
	s = file.readline().replace('\r\n','\n')

    return accu
    
def readMIDAS(filename):    
    data = MIDASReader(filename)
    lpfacts = data.getTermSet()
    return lpfacts

class MIDASReader(object):
    def __init__(self, filepath):
        self.parseMIDAS(filepath)

    def parseMIDAS(self, filepath):
        with open(filepath, 'rbU') as f:
            reader = csv.DictReader(f)

            #Header without the CellLine column
            species = reader.fieldnames[1:]

            stimulus = map(lambda name: name[3:], filter(lambda name: name.startswith('TR:') and not name.endswith('i'), species))
            inhibitors = map(lambda name: name[3:-1], filter(lambda name: name.startswith('TR:') and name.endswith('i'), species))
            readouts = map(lambda name: name[3:], filter(lambda name: name.startswith('DV:'), species))

            exp_conditions = []
            exp_observations = []
            time_points = []
            for row in reader:
                cond = {}
                for s in stimulus:
                    cond[s] = int(row['TR:' + s])
                for i in inhibitors:
                    #In MIDAS, inhibitors are set to 1 when the inhibitor is present meaning that the node it's inhibited
                    cond[i] = (int(row['TR:' + i + 'i']) + 1) % 2

                obs = {}
                for r in readouts:
                    #Ignore NaN values
                    if not math.isnan(float(row['DV:' + r])):
                        t = int(row['DA:' + r])
                        if not obs.has_key(t):
                            if t not in time_points:
                                time_points.append(t)
                            obs[t] = {}

                        obs[t][r] = float(row['DV:' + r])

                if cond in exp_conditions:
                    index = exp_conditions.index(cond)
                    exp_observations[index].update(obs)
                else:
                    exp_conditions.append(cond)
                    exp_observations.append(obs)
                
            self.stimulus = map(lambda n: n, stimulus)
            self.inhibitors = map(lambda n: n, inhibitors)
            self.readouts = map(lambda n: n, readouts)

            time_points.sort()
            if len(time_points) == 1:
                self.time_point = time_points[0]
            else:
                self.time_point = time_points[1]
                
            self.exp_conditions = exp_conditions
            self.exp_observations = exp_observations    
                    
    def discrete(self, threshold=0.5):
        """
        Returns the 2-discrete experimental observations.
        If in some time-point the value is greater or equal the threshold, the value is set to 1. Otherwise, 0.
        """
        discrete_observations = []
        for o in self.exp_observations:
            discrete_obs = {}
            obs = o[self.time_point]
            for r, value in obs.iteritems():
                if value >= threshold:
                    discrete_obs[r] = 1
                else:
                    discrete_obs[r] = 0
        
            discrete_observations.append(discrete_obs)

        return discrete_observations
                
    def multi_discrete(self, factor):
        """
        Returns the multi-discrete experimental observations.
        Takes the first time-point value that is not 0 times the given factor and then takes the int (floor) value.
        """
        discrete_observations = []
        for o in self.exp_observations:
            discrete_obs = {}
            obs = o[self.time_point]
            for r, value in obs.iteritems():
                discrete_obs[r] = int(value * factor)

            discrete_observations.append(discrete_obs)

        return discrete_observations
             
    def real(self):        
        return [o[self.time_point] for o in self.exp_observations]
    
    def getTermSet(self, factor=100):
        lpfacts = TermSet()
        lpfacts.add(Term('nexp', [len(self.exp_conditions)]))

        for index, exp in enumerate(zip(self.exp_conditions, self.multi_discrete(factor))):  
            cond, obs = exp
            for name,value in cond.iteritems():
                lpfacts.add(Term('exp', [index + 1, '"'+name+'"', value]))
    
            for name,value in obs.iteritems():
                lpfacts.add(Term('obs', [index + 1, '"'+name+'"', value]))
                

        for s in self.stimulus:
            lpfacts.add(Term('stimuli', ['"'+s+'"']))
        for i in self.inhibitors:
            lpfacts.add(Term('inhibitor', ['"'+i+'"'])) 
        for r in self.readouts:
            lpfacts.add(Term('readout', ['"'+r+'"']))  
        
        return lpfacts
