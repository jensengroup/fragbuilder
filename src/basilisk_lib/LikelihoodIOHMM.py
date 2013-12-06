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

from numpy import zeros, exp, log, transpose, sum, alltrue

one = 0.99

# Likelihood calculator for IOHMM using the Forward algorithm
class LikelihoodIOHMM_Fw:
    """
    Likelihood calculator for IOHMM using the Forward algorithm
    """
    def __init__(self, dbn, hidden_node):
        self.dbn=dbn
        # List of node objects
        nodes_0=self.nodes_0=dbn.nodes_0
        nodes_1=self.nodes_1=dbn.nodes_1

        self.hidden_node = hidden_node
        self.hidden_node_index_cpd = 0
        self.hidden_node_size = nodes_1[hidden_node].node_size

        self.forward = None
        self.parentmap_list = None

        # Assert that the hidden node has same size in both slices
        assert(nodes_0[hidden_node].node_size==nodes_1[hidden_node].node_size)

    def set_data_in_dbn(self, seq):
        dbn = self.dbn
        parentmap_list = ([],[])

        # Put the sequence into the dbn
        for node in dbn.unique_nodes:
            parentmap=node.get_parentmap(seq, weight=1)
            node.set_parentmap(parentmap)
            parentmap_list[0].append(parentmap)

    def calc_ll(self, seq, mismask, start=None, end=None, ignore_child_mism=True, set_data_in_dbn=True):
        # Make some local copies for speed.
        dbn = self.dbn
        nodes_0 = self.nodes_0
        nodes_1 = self.nodes_1
        hidden_node = self.hidden_node
        hnic = self.hidden_node_index_cpd
        hidden_node_size = self.hidden_node_size

        # Assert that seq and mismask has same shape
        assert(seq.shape == mismask.shape)

        # Put the data into the dbn
        if set_data_in_dbn:
            self.set_data_in_dbn(seq)
        
        # Set the start- and end-variable correct
        if start == None:
            start = 0

        if end == None:
            end = seq.shape[0]

        # Set the rest of the variables
        slice_count = (end-start)
        scales = zeros(slice_count)

        # Calculate the forward array
        forward = zeros((slice_count, hidden_node_size))
        for i,l in enumerate(xrange(start, end)):
            # For the first slice
            if l==0 or i==0:
                if l==0:
                    node=nodes_0[hidden_node]
                    nodes = nodes_0
                else:
                    node=nodes_1[hidden_node]
                    nodes = nodes_1

                # Get the values the parents
                pv = tuple(node.parentmap[l][:-1])

                # Get the transistion probability
                trans_prob = node.cpd[pv]
                forward[i] = trans_prob
            
            # For all following slices
            else:
                node=nodes_1[hidden_node]
                nodes = nodes_1

                # Get the values the parents
                pv = list(node.parentmap[l][:-1])

                # Slice out the previous value of the hidden node
                pv[hnic] = slice(None)
                pv = tuple(pv)

                # Get the transistion probability
                trans_prob = node.cpd[pv]

                # Multiply transistionsvalues g->h by the
                # forwardvalues [i-1, g] and sum over nodevalues of g
                forward[i] = sum(transpose(trans_prob)*forward[i-1], axis=1)

            # Multiply the forward values by the probability of the observed children
            # Note that the value of the hidden node is preserved in prev_node_value
            prev_node_value = node.parentmap[l][-1]
            for child in node.children_1:
                if ignore_child_mism or mismask[l,child.node_index] < one:
                    for j in xrange(hidden_node_size):
                        node.parentmap[l] = j
                        forward[i,j] *= exp( child.get_slice_log_likelihood(l)[0] )
            node.parentmap[l] = prev_node_value

            # Scale forward to be a probability distribution (this is
            # allowed according to BSA sec. 3.6 (p. 78)
            s=sum(forward[i])
            forward[i]/=s
            scales[i] = s


        self.forward = forward

        return sum(log(scales)), slice_count
