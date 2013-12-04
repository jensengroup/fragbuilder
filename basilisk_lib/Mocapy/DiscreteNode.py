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
Discrete node.
"""
from types import FloatType
from bisect import bisect
from random import randint

from numpy import *
from numpy.random import random


try:
    import mpi
except ImportError:
    # Replace mpi module with dummy module
    # to run on a single processor without MPI
    import DummyMPI
    mpi=DummyMPI

from MocapyExceptions import *
from ParentMap import ParentMap
from Node import Node


# This avoids log underflow
_MIN_TRANSITION=1e-50
_LOG_MIN_TRANSITION=log(_MIN_TRANSITION)


def normalize_cpd(cpd):
    """
    Normalize a CPD matrix (ie. so the sum of its values
    along the last axis is one).

    @param cpd: conditional probability distribution of a discrete node
    @type cpd: numpy array

    @return: normalized CPD
    @rtype: numpy array
    """
    s=sum(cpd+_MIN_TRANSITION, -1)
    if len(s.shape)==0:
        # CPD is a 0D array - sum is a float
        return cpd/s
    else:
        # CPD is an ND array - sum is an array
        return cpd/s[..., newaxis]

def make_random_cpd(sizes, no_zeros=1):
    """
    Return a random CPD.

    @param sizes: size of the CPD, the shape of a numpy array
    @type sizes: tuple of int

    @param no_zeros: if 1, make sure the CPD has no 0 entries
    @type no_zeros: int (0 or 1)

    @return: random CPD
    @rtype: numpy array
    """
    cpd=random(sizes)
    # avoid zero entries
    if no_zeros:
        cpd=clip(cpd, 0.1, 2)
    cpd=normalize_cpd(cpd)
    return cpd


class DiscreteNode(Node):
    """
    Discrete node.
    """
    def __init__(self, node_size, user_cpd=None, prior=None, name=""):
        """
        @param node_size: number of values the node can adopt
        @type node_size: int

        @param user_cpd: optional user supplied CPD table
        @type user_cpd: numpy array

        @param prior: a prior object
        @type prior: object from L{DiscretePriors}
        """
        # Number of values the node can adopt
        self.node_size=node_size
        # Use custom CPD if present
        self.user_cpd=user_cpd
        self.arange_node_size=arange(0, node_size)
        self.prior=prior
        self.name=name
        self.counts=None
        # Parameter count
        Node.__init__(self, output_size=1, node_type="DISCRETE")

    # Special

    def __repr__(self):
        return "<DiscreteNode name=%s size=%i\n%s\n>" % \
                (self.name, self.node_size, self.par2str())

    # Private

    def _test_cpd(self, cpd, cpd_shape):
        """
        Check if user supplied CPD table has the right shape.
        If the shape is wrong, a L{MocapyDiscreteException} is raised.

        @param cpd: CPD table
        @type cpd: numpy array

        @param cpd_shape: CPD table shape
        @type cpd_shape: tuple of ints
        """
        if len(cpd_shape)!=len(cpd.shape):
            raise MocapyDiscreteException, "User CPD has wrong dimensions."
        for i in range(0, len(cpd_shape)):
            dim1=cpd_shape[i]
            dim2=cpd.shape[i]
            if dim1!=dim2:
                raise MocapyDiscreteException, "User CPD has wrong dimensions."

    # Public

    def get_ess_list(self):
        return self.ess_list

    def get_parentmap(self, seq, weight=1):
        return ParentMap(seq, self.parents_0, self.parents_1, 
                         self.data_index, weight)

    def construct(self, slice):
        self.cpd_shape=(self.parents_0_sizes+self.parents_1_sizes+[self.node_size])
        # Parameter count
        self.par_count=multiply.reduce(self.cpd_shape)
        if not self.user_cpd is None:
            # Use user-provided CPD
            self.set_cpd(self.user_cpd)
            del self.user_cpd
        else:
            # Use random CPD
            cpd=make_random_cpd(self.cpd_shape)
            self.set_cpd(cpd)
        # Get ready for E-step
        self.ess=zeros(self.cpd_shape, 'd')
        self.ess_list=[zeros(self.cpd_shape, 'd')]
        Node.construct(self, slice)

    def set_cpd(self, new_cpd=None):
        if mpi.rank==0:
            assert(not new_cpd is None)
            cpd=normalize_cpd(new_cpd)
            # Get rid of very small values that lead to log overflow
            cpd=clip(cpd, _MIN_TRANSITION, 1000) 
            # Calculate log(cpd) for speed reasons
            log_cpd=log(cpd)
            cum_cpd=cumsum(cpd, -1)
            # Avoid bisect errors
            cum_cpd[...,-1]=2.0
            self.cpd, self.log_cpd, self.cum_cpd=mpi.bcast((cpd, log_cpd, cum_cpd))
            if __debug__:
                self._test_cpd(self.cpd, self.cpd_shape)
        else:
            self.cpd, self.log_cpd, self.cum_cpd=mpi.bcast()

    def sample(self, l):
        # Cast to tuple to avoid fancy indexing in numpy
        pv=tuple(self.parentmap[l][:-1])
        if len(pv)>0:
            cum_cpd=self.cum_cpd[pv]
        else:
            cum_cpd=self.cum_cpd
        #choice=very_simple_sampler(self.arange_node_size, cum_cpd)
        r=random()
        choice=bisect(cum_cpd, r)
        self.parentmap[l]=choice

    def update_ess(self, par_this_values):
        self.ess[par_this_values]+=1

    def save_ess(self):
        # Take sequence weight into account!
        self.ess_list[0]=self.ess_list[0]+self.weight*self.ess
        # Get ready for new ESS
        self.ess=zeros(self.cpd_shape, 'd')

    def blanket_sample(self, l):
        # local copies
        seq_len=self.seq_len
        cpd=self.cpd
        log_cpd=self.log_cpd
        children_1=self.children_1
        children_2=self.children_2
        parentmap=self.parentmap
        # Get parent+this values
        ptv=self.parentmap[l]
        # Probabilities
        p=zeros(self.node_size, 'd')
        # Loop over all proposed new states
        for state in range(0, self.node_size):
            # Set value of node in slice_1 to state
            ptv[-1]=state
            parentmap[l]=state
            # Not that we are working in LOG space!
            p[state]=log_cpd[tuple(ptv)]
            log_lik_product=0
            # Children in current slice
            for child in children_1:
                lik, cptv=child.get_slice_log_likelihood(l)
                log_lik_product+=lik
            # Children in next slice
            if l<(seq_len-1):
                # Skip this when we reach end of sequence  
                # (because there are no children in the next slice)
                for child in children_2:
                    lik, cptv=child.get_slice_log_likelihood(l+1)
                    log_lik_product+=lik
            p[state]+=log_lik_product
        # Get rid of potentially very large values in p
        p=p-max(p)
        # Get rid of potentially very small values in p
        p=clip(p, _LOG_MIN_TRANSITION, 1000)
        # Get out of LOG space
        p=exp(p)
        # Normalize probabilities
        cpd=p/sum(p)
        # Make cumulative cpd 
        cum_cpd=cumsum(cpd, -1)
        # Make sure cum_cpd sums to 1
        cum_cpd=cum_cpd/cum_cpd[-1]
        # Avoid problem with bisect
        cum_cpd[-1]=2
        # Choose a node value
        r=random()
        choice=bisect(cum_cpd, r)
        parentmap[l]=choice

    def get_slice_log_likelihood(self, l):
        ptv=tuple(self.parentmap[l])
        p=self.log_cpd[ptv]
        return p, ptv 

    def do_M_step(self, ess_list):
        if mpi.rank==0:
            # sum ESS
            ess=ess_list.pop()
            for a in ess_list:
                ess+=a
            # Save the real counts
            self.counts=array(ess)
            if self.prior:
                ess=self.prior.apply_prior(ess)
            self.set_cpd(ess)   
        else:
            self.set_cpd()
        # Get ready for new E-step
        self.ess_list=[zeros(self.cpd_shape, 'd')]

    def get_counts(self):
        return self.counts

    def get_parameters(self):
        """
        Return CPD table.

        @return: the node's CPD table
        @rtype: numpy array
        """
        return self.cpd

    def par2str(self, precision=2):
        s=array2string(self.cpd, precision=precision, max_line_width=200, 
                       suppress_small=1)
        return s
