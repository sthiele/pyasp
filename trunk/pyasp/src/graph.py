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

import sys
import re



class Graph:
    def __init__(self):
        self.vertices = set()
        self.pred = {}
        self.succ = {}
        self.elabels = {}

    def __repr__(self):
        buf = str(self.nb_vertices()) + " vertices\n"
        for v in self.vertices:
            for w in self.succ[v]:
                buf += v + " -- " + str(self.elabel(v,w)) + " --> " + w + "\n"
        return buf

    def add_vertex(self,v):
        if v not in self.vertices:
            self.vertices.add(v)
            self.pred[v] = set()
            self.succ[v] = set()

    def remove_vertex(self,v):
        self.vertices.remove(v)
        del self.pred[v]
        del self.succ[v]
        for w in self.vertices:
            if v in self.pred[w]: self.pred[w].remove(v)
            if v in self.succ[w]: self.succ[w].remove(v)
            if (v,w) in self.elabels: del self.elabels[(v,w)]
            if (w,v) in self.elabels: del self.elabels[(w,v)]

    def add_edge(self,v,w,s = None):
        self.add_vertex(v)
        self.add_vertex(w)
        self.pred[w].add(v)
        self.succ[v].add(w)
        if s: self.elabels[(v,w)] = s
                        
    def elabel(self,v,w):
        if w in self.succ[w]: return '0'
        else:
            try: return self.elabels[(v,w)]
            except KeyError: return '?'

    def nb_vertices(self):
        return len(self.vertices)

    def nb_edges(self):
        accu = 0
        for v, succ in self.succ.iteritems():
            accu += len(succ)
        return accu

    def __eq__(self,other):
        return \
            self.vertices == other.vertices \
        and self.pred == other.pred \
        and self.succ == other.succ \
        and self.elabels == other.elabels
    

class ErdosRenyiGraph(Graph):
    """
    Random graph generated under the assumption that edges are present
    with given probability, independently of the others. 

    Arguments:
      n : int, number of vertices
      p : float, probability for an edge to be present
      autoloops : bool (optional), specifies whether loops on one vertex are allowed (default False)
    """
    def __init__(self, n, p, autoloops = False):
        Graph.__init__(self)
        for i in range(n):
            self.add_vertex(str(i))
            for j in range(n):
                if i == j and not autoloops: pass
                if random.random() < p: 
                    if random.random() > 0.5: s = '-'
                    else: s = '-'
                    self.add_edge(str(i),str(j),s)



def rien():
    print sys.prefix


class Profile(dict):
    """
    Basically a vertex-to-sign mapping, which is used as a labelling of the 
    vertices of a graph
    """
    def __init__(self,filename):
        """Taken from bioquali, IRISA Rennes"""
        COMPLEX_ID = '[-a-zA-Z0-9_:\(\)/&]+'
        GENE_ID = '[-a-zA-Z0-9_:\(\)/]+'
        file = open(filename,'r')
        val_re = '(?P<genid>'+GENE_ID+')( )*=( )*(?P<sign>[-+0])'
        val = re.compile(val_re)
        name_re = '(-\*-)( )*(?P<name>\S+)( )*(-\*-)'
        name = re.compile(name_re)
        obs_name = None

        line_number = 1
        line = file.readline()[:-1]
        # the first line may contain the name
        nm = name.match(line)
        if nm:
            self.name = nm.group('name')
            line = file.readline()[:-1]
        # ordinary line
        while line:
            vm = val.match(line)
            if vm:
                vertex = vm.group('genid')
                self[vertex] = vm.group('sign')
            else:
                print 'Syntax error line:', line_number,'  '+line 

            line = file.readline()[:-1]
            line_number+=1
