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
Abstract base class for all L{Node} classes.
"""

from numpy import log, exp

try:
    import mpi
except ImportError:
    # Replace mpi module with dummy module
    # to run on a single processor without MPI
    import DummyMPI
    mpi=DummyMPI

from MocapyExceptions import *


class Node:
    """
    Abstract base class for all L{Node} objects.
    """
    def __init__(self, output_size, node_type):
        """
        @param output_size: dimension of output vector
        @type output_size: int

        @param node_type: type of node (ie. DISCRETE)
        @type node_type: string
        """
        self.output_size=output_size
        # Number of values the node can adopt
        # Children in same slice
        self.children_1=[]
        # Children in next slice
        self.children_2=[]
        # Parent indices in previous slice
        self.parents_0=[]
        # Sizes of parents in previous slice
        self.parents_0_sizes=[]
        # Parents indices in same slice
        self.parents_1=[]
        # Sizes of parents in same slice
        self.parents_1_sizes=[]
        self.fixed=0
        self.is_constructed=0
        self.node_type=node_type

    # Public

    def set_node_index(self, ni):
        """
        Set the index of the node in the node list.

        @param ni: index of node (ie. 0=first node in slice)
        @type ni: int
        """
        self.node_index=ni

    def set_data_index(self, di):
        """
        Set the index of where to find the node data in the 
        sequence data.

        @param di: index of node data in a slice
        @type di: int
        """
        self.data_index=di

    def add_intra_child(self, node):
        """
        Add an intra-slice child.

        @param node: child node (same slice) of current node 
        @type node: L{Node}
        """
        self.children_1.append(node)

    def add_inter_child(self, node):
        """
        Add a child in the previous slice.

        @param node: child node (previous slice) of current node
        @type node: L{Node}
        """
        self.children_2.append(node)

    def add_inter_parent(self, data_index, node_size):
        """
        Add a parent in the previous slice.
        
        @param data_index: index of parent value in previous slice
        @type data_index: int

        @param node_size: output size of parent
        @type node_size: int
        """
        self.parents_0.append(data_index)
        self.parents_0_sizes.append(node_size)

    def add_intra_parent(self, data_index, node_size):
        """
        Add an intra-slice parent.

        @param data_index: index of parent value in same slice
        @type data_index: int

        @param node_size: output size of parent
        @type node_size: int
        """
        self.parents_1.append(data_index)
        self.parents_1_sizes.append(node_size)

    def get_parentmap(self, seq, weight):
        """
        Translate a sequence into a FamilyMap and return 
        the FamilyMap object. The FamilyMap object is a
        cache for fast parent and child lookup.

        @param seq: data sequence
        @type seq: numpy array

        @return: family map data structure which implements fast child
            and parent value lookup
        @rtype: FamilyMap object
        """
        raise NotImplementedError

    def set_parentmap(self, parentmap):
        """
        Load a L{ParentMap} object into the  L{Node} object. The L{ParentMap} object 
        is a cache for fast parent and child value lookup. Loading a 
        L{ParentMap} object thus means loading new sequence data.

        @param parentmap: load a L{ParentMap} object into the L{Node} object
        @type parentmap: L{ParentMap}
        """

        self.parentmap=parentmap
        self.weight=parentmap.get_weight()
        self.seq_len=len(parentmap)

    def construct(self, slice):
        """
        Initialize the data structures of the node.
        'slice' indicates to which slice the node belongs
        (ie. slice l=0 or slice l>0).

        @param slice: slice (0=slice at l=0, 1=slice at l>0)
        @type slice: int
        """
        if self.is_constructed==1:
            raise MocapyException, "a Node should only be constructed once."
        self.slice=slice
        self.is_constructed=1

    def sample(self, l):
        """
        Sample the node at slice l. 
        This is done by using P(node|parents), ie. not 
        taking the children into account.
        
        @param l: sequence position
        @type l: int
        """
        raise NotImplementedError

    def update_ess(self, family_values):
        """
        After sampling a hidden node, the ESS are updated using
        this function. This is done using the values of the node itself and its
        parents (ie. the family values).
        
        @param family_values: values of parents and node itself
        @type family_values: numpy array
        """
        raise NotImplementedError

    def save_ess(self):
        """
        Calculate the ESS for a whole sequence, save it, and re-initialize
        the ESS data structure. This function is called after a whole sequence
        has been sampled.
        """ 
        raise NotImplementedError

    def blanket_sample(self, l):
        """
        Sample the node value at slice l based on its
        Markov blanket. The Markov blanket consists of
        the parent and child nodes of the node.

        @param l: sequence position
        @type l: int
        """
        raise NotImplementedError

    def get_slice_log_likelihood(self, l):
        """
        Return the log likelihood of the node in slice l.
        @param l: sequence position
        @type l: int

        @return: log likelihood of nodes in slice l
        @rtype: float
        """
        raise NotImplementedError

    def get_ess_list(self):
        """
        Return the list of ESS values.
        This method is used to gather all ESS calculated
        on different nodes when running in parallel. The
        ESS are put in one list, and passed on to the 
        'do_M_step' method.

        @return: list of ESS
        @rtype: list
        """
        raise NotImplementedError

    def do_M_step(self, ess_list):
        """
        Update the parameters of the node using the ESS values.

        @param ess_list: list of ESS
        @type ess_list: list
        """
        raise NotImplementedError

    def get_parameters(self):
        """
        Return the relevant parameters of the node.
        """
        raise NotImplementedError

    def fix(self, flag):
        """
        Fix the node (flag=1), or unfix the node (flag=0).

        @param flag: fix flag
        @type flag: int (0 or 1)
        """
        self.fixed=flag
