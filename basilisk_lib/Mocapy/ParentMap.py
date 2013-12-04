#   Mocapy: A parallelized Dynamic Bayesian Network toolkit
#
#   Copyright (C) 2004 Thomas Hamelryck 
#
#   This library is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with BASILISK.  If not, see <http://www.gnu.org/licenses/>.

"""
Fast parent node value/node value lookup for discrete nodes.
"""


from numpy import take, multiply

from MocapyExceptions import MocapyException


class ParentMap:
    """
    L{ParentMap} keeps track of where the data describing the values of a specific 
    node is found in the data sequence.
    """
    def __init__(self, seq, parents_0, parents_1, this_index, weight):
        """
        @param seq: the data sequence
        @type seq: numpy array

        @param parents_0: index of parent value in previous slice 
        5B
        @type parents_0: int

        @param parents_1: index of parent value in current slice 
        @type parents_1: int

        @param this_index: index of node value itself in current slice
        @type this_index: int

        @param weight: sequence weight
        @type weight: float in [0,1]
        """
        # data
        lng=len(seq)
        self.lng=lng
        self.seq=seq.flat
        # self index
        # Cache the indices
        indices=[]
        # Step size
        if len(seq.shape)==1:
            step=1
        else:
           step=multiply.reduce(seq.shape[1:])
        for l in range(0, lng):
            index=[]
            if l>0:
                for pi in parents_0:
                    index.append((l-1)*step+pi)
            for pi in parents_1:
                index.append(l*step+pi)
            index.append(l*step+this_index)
            indices.append(index)
        self.this_index=this_index
        self.step=step
        self.indices=indices
        self.weight=weight

    def __setitem__(self, l, v):
        """
        Set data value of a node at position l.

        @param l: sequence index
        @type l: int

        @param v: data value
        @type v: array of type 'd'
        """
        index=self.indices[l][-1]
        self.seq[index]=v

    def __getitem__(self, l):
        """
        Return the parent/node values at position 'l'.

        @param l: sequence position
        @type l: int

        @return: parent/node values at position 'l' 
        @rtype: array, consisting of [parent value 1, parent value 2,..., node value]
        """
        index=self.indices[l]
        ptv=take(self.seq, index).astype('i')
        return ptv

    def __getstate__(self):
        """
        Dummy method to avoid pickling data.
        """
        return {}

    def __len__(self):
        return self.lng

    def get_weight(self):
        """
        Return the sequence weight.

        @return: sequence weight
        @rtype: float in [0,1]
        """
        return self.weight

    def replace_seq(self, seq):
        """ 
        Replace the sequence - used in undo operations.
        """
        if not len(seq)==self.lng:
            raise MocapyException, "Replacing sequence has wrong shape"
        self.seq=seq.flat


class ContinuousParentMap(ParentMap):
    def __init__(self, seq, parent, this_index, dim, weight):
        self.dim=dim
        ParentMap.__init__(self, seq, [], [parent], this_index, weight)

    def __getitem__(self, l):
        i,j=self.indices[l]
        parent_value=int(self.seq[i])
        node_value=self.seq[j:j+self.dim]
        return (parent_value, node_value)

    def __setitem__(self, l, v):
        index=self.indices[l][-1]
        self.seq[index:index+self.dim]=v
 
if __name__=="__main__":

    from numpy import *

    seq=array(((1,2,3),
               (3,4,5)), 'd')

    parents_0=[0,1]
    parents_1=[0,1]

    f=ParentMap(seq, parents_0, parents_1, 2, 0.0)

    print f[0]
    print f[1]

    f[0]=30
    f[1]=60

    print f[0]
    print f[1]

