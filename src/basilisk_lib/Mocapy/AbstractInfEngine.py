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
Abstract base class for inference engines. 
"""

from numpy import zeros, array, sum, log, exp
from ArrayMisc import argmax


class AbstractInfEngine:
    """
    Abstract base class for inference engines.
    """

    def __init__(self, dbn, seq, weight):
        """
        @param dbn: DBN definition
        @type dbn: L{DBN} object

        @param seq: sequence data
        @type seq: numpy array

        @param weight: sequence weight
        @type weight: float
        """
        self.dbn=dbn
        self.seq=seq
        self.weight=weight
        self.seq_len=seq.shape[0]
        self.nodes_0, self.nodes_1=dbn.get_nodes()
        self.output_size=dbn.total_output_size

    def __len__(self):
        """
        Returns sequence length.

        @return: sequence length
        @rtype: int
        """
        return self.seq_len

    def _make_parentmap_list(self, seq, weight):
        """
        Construct a list of L{ParentMap} objects from the 
        sequence. 

        @param seq: sequence
        @type seq: numpy array

        @param weight: sequence weight
        @type weight: float
        
        @return: parent map for each node
        @rtype: list of L{ParentMap} objects
        """
        l=[]
        for node in self.dbn.unique_nodes:
            parentmap=node.get_parentmap(seq, weight)
            l.append(parentmap)
        return l

    def _set_parentmap(self, parentmap_list):
        """
        Put previously calculated L{ParentMap} objects (ie. for fast node value
        lookup in sequences) into the nodes.

        @param parentmap_list: list of parent maps
        @type parentmap_list: list of L{ParentMap} objects
        """
        for i in range(0, len(self.dbn.unique_nodes)):
            node=self.dbn.unique_nodes[i]
            parentmap=parentmap_list[i]
            assert(len(parentmap)==self.seq_len)
            node.set_parentmap(parentmap)

    def get_viterbi(self):
        """
        Calculate the Viterbi path. This method can act on a set of samples 
        (resulting in a stochastic Viterbi path calculation) or on all possible
        hidden node values for each slice (resulting in a deterministic Viterbi
        path calculation). The former is done by L{InfEngineMCMC}, the latter 
        by L{InfEngineHMM}.

        @return: (viterbi path, Loglik) tuple
        @rtype: (numpy array, float)
        """
        dbn=self.dbn
        seq=self.seq
        nr_slices=self.nr_slices
        seq_len=self.seq_len
        nodes_0=self.nodes_0
        nodes_1=self.nodes_1
        slices=self.slices
        core_seq=zeros(seq.shape, 'd')
        parentmap_list=self._make_parentmap_list(core_seq, self.weight)
        self._set_parentmap(parentmap_list)
        viterbi=zeros((seq_len, self.output_size), 'd')
        gain=zeros((seq_len, nr_slices), 'd')
        path=zeros((seq_len-1, nr_slices), 'i')
        # Initialization
        # l=0
        for i in range(0, nr_slices):
            slice=slices[0,i]
            core_seq[0]=slice
            ll=0.0
            for node in nodes_0:
                slice_log_lik, fv=node.get_slice_log_likelihood(0)
                ll+=slice_log_lik
            gain[0,i]=ll
        # Propagation
        # l>0
        for l in range(1, seq_len):
            for j in range(0, nr_slices):
                slice=slices[l,j]
                core_seq[l]=slice
                m=zeros(nr_slices, 'd')
                for i in range(0, nr_slices):
                    core_seq[l-1]=slices[l-1,i]
                    ll=0
                    for node in nodes_1:
                        slice_log_lik, fv=node.get_slice_log_likelihood(l)
                        ll+=slice_log_lik
                    m[i]=gain[l-1,i]+ll
                index=argmax(m)
                gain[l,j]=m[index]
                path[l-1,j]=index
        # Termination
        index=argmax(gain[seq_len-1,:])
        viterbi[seq_len-1,:]=slices[seq_len-1, index]
        # Backtracking
        for l in range(seq_len-2,-1,-1):
            # Here's some Numpy weirdness! 
            # path[a,b] returns an array of size (1,)!!!
            # This is so if it is of type int, not if type float!
            index=path[l,index][0]
            viterbi[l,:]=slices[l, index]
        # REMARK: remember that dbn.calc_ll puts the Viterbi sequence 
        # IN THE NODES so we make a copy of the Viterbi path to avoid 
        # side effects
        #ll=dbn.calc_ll(array(viterbi))
        #print "M=",max(gain[seq_len-1,:])/seq_len
        ll=max(gain[seq_len-1])/seq_len
        return viterbi, ll


if __name__=="__main__":

    pass
