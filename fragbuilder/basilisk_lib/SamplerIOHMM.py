#	 BASILISK: A generative, probabilistic model of side chains in proteins
#
#	 Copyright (C) 2010		Tim Harder and Jes Frellsen 
#
#	 BASILISK is free software: you can redistribute it and/or modify
#	 it under the terms of the GNU General Public License as published by
#	 the Free Software Foundation, either version 3 of the License, or
#	 (at your option) any later version.
#
#	 BASILISK is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
#	 GNU General Public License for more details.
#
#	 You should have received a copy of the GNU General Public License
#	 along with BASILISK.	 If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

from Mocapy.MocapyExceptions import MocapyException
from random import random
from bisect import bisect
from numpy import zeros, cumsum, transpose, sum, exp, log

# those make the comparision a little easier
one = 0.9
zero = 0.1

class SamplerFwBt:
    """
    Sampler for a multiple-input/multiple-output HMM with only one
    hidden DISCRETE node using forward backtrack.
    """
    def __init__(self, dbn, hidden_node):
        # Some copies for persistence.
        self.dbn=dbn

        # List of node objects
        nodes_0=self.nodes_0=dbn.nodes_0
        nodes_1=self.nodes_1=dbn.nodes_1
        self.nr_nodes=dbn.nr_nodes
        self.total_output_size=len(self.nodes_0)

        self.hidden_node = hidden_node

        # The position index in the cpd of the hidden node in slice 1
        # for the hidden node in the previous slice
        self.hidden_node_index_cpd = self.nodes_1[hidden_node].parents_0.index(hidden_node)
        self.hidden_node_size = self.nodes_1[hidden_node].node_size

    def sweep(self, mismask, start, end):
        # Make some local copies for speed.
        nodes_0 = self.nodes_0
        nodes_1 = self.nodes_1
        nr_nodes = self.nr_nodes
        hidden_node = self.hidden_node
        hidden_node_size = self.hidden_node_size
        seq_len = len(mismask)
        hnic = self.hidden_node_index_cpd
        
        assert(nodes_0[hidden_node].node_size==nodes_1[hidden_node].node_size)
        node_size = int(nodes_0[hidden_node].node_size)
        
        slice_count=(end-start)
        trans_probs = []
        
        # [BASILISK - SPECIFIC CODE]
        # if all the child nodes are unobserved (backbone indendent sampling)
        # we can save some time and just do ancestral sampling here. 
        if (int( sum(mismask[:,2])) == seq_len ) :
            for l in xrange(0, seq_len):
                if l==0:
                    node=nodes_0[hidden_node]
                else:
                    node=nodes_1[hidden_node]
                # Sample the hidden node
                node.sample(l)

                # Sample the children of the hidden node
                for child in node.children_1:
                    child.sample(l)   
            return slice_count
        # [/BASILISK - SPECIFIC CODE]
        
        # Calculate the forward array
        forward = zeros((slice_count, node_size))
        for i,l in enumerate(xrange(start, end)):
            # For the first slice
            if l==0 or i==0:
                if l==0:
                    node=nodes_0[hidden_node]
                else:
                    node=nodes_1[hidden_node]

                # Get the values of the parrents
                pv=tuple(node.parentmap[l][:-1])

                # Get the transistion probability (and cache it for backtrack)
                trans_prob = node.cpd[pv]
                trans_probs.append(trans_prob)

                forward[i] = node.cpd[pv]

            # For all following slices
            else:
                node=nodes_1[hidden_node]

                # Get the values of the parents except the value of the hidden node in previous slice
                pv = list(node.parentmap[l][:-1])
                pv[hnic] = slice(None)
                pv = tuple(pv)

                # Get the transistion probability (and cache it for backtrack)
                trans_prob = node.cpd[pv]
                trans_probs.append(trans_prob)

                # Multiply transistionsvalues g->h by the
                # forwardvalues [i-1, g] and sum over nodevalues of g
                forward[i] = sum(transpose(trans_prob)*forward[i-1], axis=1)

            # Multiply the forward values by the probability of the observed children
            # Note that the value of the hidden node is preserved in prev_node_value
            for child in node.children_1:
                if mismask[l,child.node_index] < one:
                    if child.node_type=="DISCRETE":
                        child_value = child.parentmap[l][-1]
                        child_dist = child.cpd[:,child_value]
                        forward[i,:] *= child_dist  
                    else:
                        org_node_value = node.parentmap[l][-1]
                        for j in xrange(hidden_node_size):
                            node.parentmap[l] = j
                            forward[i,j] *= exp( child.get_slice_log_likelihood(l)[0] )
                        node.parentmap[l] = org_node_value

            # If the is the end slice but not the end of the DBN, take the next hidden value into account
            if l==end-1 and l!=seq_len-1:
                node_next = nodes_1[hidden_node]
                pv_next = list(node_next.parentmap[l+1])
                pv_next[hnic] = slice(None)
                forward[i] *= node_next.cpd[tuple(pv_next)]
                                    
            # Scale forward to be a probability distribution (as
            # described in 'Biological Sequence Analysis' by R Durbin
            # et al., p. 78, Cambridge University Press, 1998.
            s=sum(forward[i])
            forward[i]/=s
                    
        # Do the backtrack
        for l in xrange(end-1, start-1,-1):
            i = l-start

            if l==0:
                node=nodes_0[hidden_node]
            else:
                node=nodes_1[hidden_node]

            # Choose whether to sample from the forward probability
            # distribution or take the previous choice into account
            if l==end-1:
                distribution = forward[i]
            else:
                proportional = forward[i] * trans_probs[i+1][:,choice]
                distribution = proportional / sum(proportional)

            # Sample the hidden node from the forward probability distribution
            cum_forward = cumsum(distribution)
            r=random()
            choice=bisect(cum_forward, r)
            node.parentmap[l]=choice

            # Sample the children
            for child in node.children_1:
                if mismask[l,child.node_index] > zero:
                    child.sample(l)

        return slice_count

