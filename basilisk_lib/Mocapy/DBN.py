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
The DBN class, plus some loading/initialisation code.
"""

import random
import pickle
from types import IntType

from numpy import zeros, sum, log, exp
import numpy.random

from MocapyExceptions import *


def mocapy_seed(a, b):
    """
    Initialize the random number generators that are used
    by Mocapy. These are numpy.random and Python's random 
    module.

    @param a: seed for numpy.random
    @type a: int

    @param b: seed for Python's random module
    @type b: int
    """
    numpy.random.seed(a)
    random.seed(b)


def load_dbn(fname):
    """
    Load a pickled DBN object.

    @param fname: filename of the pickled DBN object
    @type fname: string

    @return: the unpickled DBN object
    @rtype: L{DBN} object
    """
    fp=open(fname, 'r')
    dbn=pickle.load(fp)
    fp.close()
    return dbn


class DBN:
    """
    This class stores the Node objects and performs the necessary initialisations. 
    In addition, it checks wether the DBN definition makes sense. The class can also
    be used to sample a sequence from the DBN, and to calculate the loglik of a 
    sequence (when all node values are known).
    """
    def __init__(self, nodes_0, nodes_1, name="DBN"):
        """
        @param nodes_0: nodes in slice 0.
        @type nodes_0: list of L{Node} objects.

        @param nodes_1: nodes in slice 1.
        @type nodes_1: list of L{Node} objects.

        @param name: optional name of the L{DBN}.
        @type name: string
        """
        # Flag - used to check proper initialisation
        self.is_constructed=0
        # Nodes in slice 0
        self.nodes_0=nodes_0
        # Nodes in slices 1-T
        self.nodes_1=nodes_1
        if len(nodes_0)!=len(nodes_1):
            raise MocapyDBNException, "Slices need to have the same sumber of nodes."
        # Number of nodes/slice
        self.nr_nodes=len(self.nodes_0)
        self.name=name
        # Set the node indices
        data_index=0
        for n in range(0, self.nr_nodes):
            n0=nodes_0[n]
            n1=nodes_1[n]
            # Identifier of the node
            n0.set_node_index(n)
            n1.set_node_index(n)
            # Where is the node data found?
            n0.set_data_index(data_index)
            n1.set_data_index(data_index)
            assert(n0.output_size==n1.output_size)
            # Move to next free index in data
            data_index+=n1.output_size
        # Map node name to index
        self.index_map=self._make_index_map(nodes_0, nodes_1)
        self.total_output_size=data_index

    def _make_index_map(self, nodes_0, nodes_1):
        """
        Create a dicationary that maps the name of a node to its index.

        @param nodes_0: nodes in slice 0.
        @type nodes_0: list of L{Node} objects.

        @param nodes_1: nodes in slice 1.
        @type nodes_1: list of L{Node} objects.

        @rtype: dictionary that maps node names to (index, slice) tuples
        """
        index_map={}
        for i in range(0, len(nodes_0)):
            node=nodes_0[i]
            if node.name!="":
                index_map[node.name]=(i, 0)
            node=nodes_1[i]
            if node.name!="":
                index_map[node.name]=(i, 1)
        return index_map
    
    def _map_to_indices(self, parent_i, child_i):
        """
        Return node indices associated with the two given node names 'parent'
        and 'child'. If a name is an integer, it is returned as is.

        @param parent_i: parent index or name
        @type parent_i: string or int

        @param child_i: child index or name
        @type child_i: string or int

        @return: (parent index, child index)
        @rtype: (int, int)
        """
        if not type(parent_i)==IntType:
            if not self.index_map.has_key(parent_i):
                raise MocapyDBNException, "Unknown node name: %s" % str(parent_i)
            else:
                parent_i, slice=self.index_map[parent_i]
        if not type(child_i)==IntType:
            if not self.index_map.has_key(child_i):
                raise MocapyDBNException, "Unknown node name: %s" % str(child_i)
            else:
                child_i, slice=self.index_map[child_i]
        return parent_i, child_i
            
    # Special

    def __repr__(self):
        """
        @return: "<DBN I{name} #nodes=I{number of nodes}>".
        @rtype: string
        """
        s="<DBN %s #nodes=%i>" % (self.name, self.nr_nodes)
        return s

    # Public

    def get_nodes(self):
        """
        Return nodes in slice l=0, and in slice l>0
        as two node lists.

        @return: (node list at l=0, node list at l>0)
        @rtype: ([L{Node},...,L{Node}], [L{Node},...,L{Node}])
        """
        return self.nodes_0, self.nodes_1

    def get_node_by_name(self, name):
        """
        Return the node object that is associated with the given name.

        @param name: name of the node
        @type name: string

        @return: node associated with name
        @rtype: L{Node}
        """
        index, slice=self.index_map[name]
        if slice==0:
            return self.nodes_0[index]
        else:
            return self.nodes_1[index]

    def sample_sequence(self, length):
        """
        Return a sampled sequence sequence from the DBN
        with specified length.

        @param length: length of data sequence
        @type length: int

        @return: a sequence with the hidden nodes sampled
        @rtype: array, shape=(sequence length, output size)
        """
        assert(self.is_constructed)
        seq=zeros(([length, self.total_output_size]), 'd')
        for node in self.unique_nodes:
            parentmap=node.get_parentmap(seq, weight=1)
            node.set_parentmap(parentmap)
        ll=0
        for l in range(0, length):
            if l==0:
                # Use nodes in slice 0
                node_list=self.nodes_0
            else:
                # Use nodes in slices 1-T
                node_list=self.nodes_1
            for node in node_list:
                node.sample(l)
                slice_log_ll, family=node.get_slice_log_likelihood(l)
                ll+=slice_log_ll
        return seq, ll/length

    def add_inter(self, parent_i, child_i):
        """
        Add an edge between slices, from parent to child.

        @param parent_i: parent node index or name 
        @type parent_i: int or string

        @param child_i: child node index or name 
        @type child_i: int or string
        """
        assert(not self.is_constructed)
        # Map names to indices if necessary
        parent_i, child_i=self._map_to_indices(parent_i, child_i)
        # slice 0
        parent_0=self.nodes_0[parent_i]
        # slice 1
        parent_1=self.nodes_1[parent_i]
        child_1=self.nodes_1[child_i]
        # Add children to parents
        parent_0.add_inter_child(child_1)
        # Note that child_0 cannot have parents
        child_1.add_inter_parent(parent_0.data_index, parent_0.node_size)
        # Parents tied?
        if not parent_0 is parent_1:
            parent_1.add_inter_child(child_1)

    def add_intra(self, parent_i, child_i):
        """
        Add an edge inside a slice, from parent to child.

        @param parent_i: parent node index or name 
        @type parent_i: int or string

        @param child_i: child node index or name 
        @type child_i: int or string
        """
        assert(not self.is_constructed)
        # Map names to indices if necessary
        parent_i, child_i=self._map_to_indices(parent_i, child_i)
        # slice 0
        parent_0=self.nodes_0[parent_i]
        child_0=self.nodes_0[child_i]
        # slice 1
        parent_1=self.nodes_1[parent_i]
        child_1=self.nodes_1[child_i]
        # Add children to parents
        parent_0.add_intra_child(child_0)
        if not parent_0 is parent_1:
            parent_1.add_intra_child(child_1)
        # Add parents to children
        child_0.add_intra_parent(parent_0.data_index, parent_0.node_size)
        if not child_0 is child_1:
            child_1.add_intra_parent(parent_0.data_index, parent_0.node_size)

    def construct(self):
        """
        Initialize the DBN data structures based on 
        the added nodes and edges. After calling this
        method no edges can be added.
        """
        assert(not self.is_constructed)
        # construct nodes in slice 0
        self.unique_nodes=[]
        for n in range(0, self.nr_nodes):
            # Slice 0
            node_0=self.nodes_0[n]
            node_1=self.nodes_1[n]
            if node_0 is node_1:                 
                # Tied node
                node_0.construct(slice="TIED")
                self.unique_nodes.append(node_0)
            else:
                node_0.construct(slice="START")
                node_1.construct(slice="END")
                self.unique_nodes.append(node_0)
                self.unique_nodes.append(node_1)
        self.all_nodes=self.nodes_0+self.nodes_1
        self.is_constructed=1

    def calc_ll(self, seq):
        """
        Calculate the loglik for a sequence (with all node values known).

        @return: LogLik (normalized for number of slices)
        @rtype: float
        """
        ll=0
        seq_len=seq.shape[0]
        for node in self.unique_nodes:
            parentmap=node.get_parentmap(seq, weight=1)
            node.set_parentmap(parentmap)
        for l in range(0, seq_len):
            if l==0:
                nodes=self.nodes_0
            else:
                nodes=self.nodes_1
            for n in range(0, self.nr_nodes):
                node=nodes[n]
                node_log_ll, fv=node.get_slice_log_likelihood(l)
                ll+=node_log_ll
        return ll/seq_len

    def get_parameter_count(self):
        """
        Return total number of parameters.

        @return: total number of parameters
        @rtype: int
        """
        pc=0
        for node in self.unique_nodes:
            pc+=node.par_count
        return pc

    def save(self, fname):
        """
        Save a persistent copy of the DBN object (ie. a pickled object).
        This persistent copy can be loaded again with 'load_dbn'.

        @param fname: filename
        @type fname: string
        """
        fp=open(fname, 'w')
        pickle.dump(self, fp)
        fp.close()


