#   An extension to: Mocapy
#
#   Copyright (C) 2007 Jes Frellsen, Ida Moltke and Martin Thiim
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
Von Mises node.
"""

from math import pi, cos, sin, acos
from numpy import zeros, concatenate, array, sum, array2string, cov,\
                  log, exp, isnan
from numpy.random import random, uniform

try:
    import mpi
except ImportError:
    # Replace mpi module with dummy module
    # to run on a single processor without MPI
    import DummyMPI
    mpi=DummyMPI

from Node import Node
from ParentMap import ContinuousParentMap
from ArrayMisc import norm
from VonMises import VMSampler, VMDens, mod2pi, to_radian
from VonMisesEstimator import estimate_kappa
from MocapyExceptions import MocapyVMException


class VMNode(Node):
    "Von Mises node for observations on a 2-dimensional sphere"
    "described by a direction in the interval from 0 to 2pi"
    def __init__(self, node_size, user_mus=None, user_kappas=None, name=""):
        """
        @param node_size: number of mixture components 
        @type node_size: int
        
        @param user_mus: user supplied mean directions
        @type user_mus: numpy array, shape=nsize x 1 

        @param user_kappas: user supplied Kappa parameters
        @type user_kappas: numpy array, shape=nsize
        """
        # Set different variable values
        self.node_size=node_size
        self.user_mus=user_mus
        self.user_kappas=user_kappas
        self.ess_r, self.ess_w=self._get_empty_ess()
        self.ess_list=[self._get_empty_ess()]
        self.name=name
        # Parameter count
        self.par_count=node_size*3
        # Initialize Node base class
        Node.__init__(self, output_size=1, node_type="VM")

    def __repr__(self):
        s="<VMNode name=%s size=%i\n%s\n>" % (self.name, self.node_size, self.par2str())
        return s

    def _get_empty_ess(self):
        """
        Return initialized ESS data structures.
        """
        ess_r=zeros((self.node_size, 2), 'd')
        ess_w=zeros(self.node_size, 'd')
        return [ess_r, ess_w]

    def _make_rnd_kms(self, node_size):
        """
        Make random Mu's and Kappa's.
        """
        mus=zeros((node_size,1), 'd')
        for i in range(0, node_size):
            mus[i]=uniform(0,2*pi)
        kappas=random(node_size)            
        return mus, kappas

    def _make_vm_list(self, mus, kappas, node_size):
        """
        Initialize the list of L{VMDens} objects.
        These objects calculate P(obs|mu, kappa), and
        are here pre-constructed for speed.
        """
        vm_list=[]
        samplers=[]
        for i in range(0, node_size):
            m=mus[i]
            k=kappas[i]
            vm=VMDens(k, m)
            vm_list.append(vm)
            sampler=VMSampler(k, m)
            samplers.append(sampler)
        return vm_list, samplers

    # Public
    def get_ess_list(self):
        return self.ess_list

    def construct(self, slice):
        Node.construct(self, slice)
        self.parent_index=self.parents_1[0]
        assert(len(self.parents_0)==0)
        assert(len(self.parents_1)==1)
        parent_size=self.parents_1_sizes[0]
        assert(self.node_size==parent_size)
        self.mus_shape=(self.node_size,1)
        self.kappas_shape=(self.node_size),
        # Check or construct means
        rnd_mus, rnd_kappas=self._make_rnd_kms(self.node_size)

        # If the user has specified some values these are used        
        if not self.user_mus is None:
            self.mus=self.user_mus
            for i in range(0, self.node_size):
                # Make sure we have angles in [0,2pi[
                self.mus[i]=mod2pi(self.mus[i])
            assert(self.mus.shape==self.mus_shape)
            del self.user_mus
        else:
            # Otherwise the means are set to random values
            self.mus=rnd_mus

        # Check or construct kappa values
        if not self.user_kappas is None:
            self.kappas=self.user_kappas
            assert(self.kappas.shape==self.kappas_shape)
            del self.user_kappas
        else:
            self.kappas=rnd_kappas
        self.vm_list, self.samplers=self._make_vm_list(self.mus, self.kappas, self.node_size)

    def get_parentmap(self, seq, weight):
        return ContinuousParentMap(seq, self.parent_index, self.data_index,
                                   1, weight)

    def update_ess(self, family_values):
        # Get which node value the parent has
        parent_value=family_values[0]
        node_value=family_values[1]
        # Add ess info on the current node value to the current info on
        # the mixture element corresponding to the parent value
        self.ess_r[parent_value]+=array([cos(node_value), sin(node_value)])
        self.ess_w[parent_value]+=1

    def save_ess(self):
        if self.weight != 1:
            raise MocapyVMException, "Sequence weighting is not supported"
        ess=self.ess_list[0]
        ess[0]+=self.weight*self.ess_r
        ess[1]+=self.weight*self.ess_w
        self.ess_r, self.ess_w=self._get_empty_ess()

    def sample(self, l):
        # Get which node value the parent node has in slice l
        parent_value, node_value=self.parentmap[l]

        # Sample from the corresponding von mises sampler
        sampler=self.samplers[parent_value]
        s=sampler()

        # Store the sampled value as the current node value

        self.parentmap[l]=s
        return s

    def blanket_sample(self,l):
        return self.sample(l)

    def get_slice_log_likelihood(self, l):
        # Get the value of the parent node in slice l
        # and the value of the node in slice l
        parent_value, node_value=self.parentmap[l]

        # Return log of the density at the node value
        # in the von mises distribution which corresponds to
        # the value of the parent node
        vm_dens=self.vm_list[parent_value]
        ll=vm_dens(node_value, log_space=1)
        return ll, (parent_value, node_value)

    def do_M_step(self, ess_list):
        # CPU 0
        if mpi.rank==0:
            mus=zeros((self.node_size,1), 'd')
            kappas=zeros(self.node_size, 'd')
            ess_r, ess_w=ess_list.pop()
            for i in range(0, len(ess_list)):
                r, w=ess_list[i]
                ess_r+=r
                ess_w+=w
            for i in range(0, self.node_size):
                r=ess_r[i]
                w=ess_w[i]
                norm_r=norm(r)

                print "do_M_step (%i): r=%s, w=%s, norm_r=%s" % (i,r,w,norm_r)

                if w == 0:
                    print "WARNING: There is no data to estimate mu and kappa from"
                    kappas[i]=self.kappas[i]
                    mus[i]=self.mus[i]
                else:
                    mu_coords=r/norm_r
                    kappas[i]=estimate_kappa(w,r,mu_coords)
                    mus[i]=to_radian(mu_coords)

            vm_list, samplers=self._make_vm_list(mus, kappas, self.node_size)
            self.mus, self.kappas, self.vm_list, self.samplers=\
                    mpi.bcast((mus, kappas, vm_list, samplers))
        # Other CPUs
        else:
            self.mus, self.kappas, self.vm_list, self.samplers=mpi.bcast()
        self.ess_list=[self._get_empty_ess()]

    def get_parameters(self):
        """
        Return the VM parameters.

        @return: Mu's, Kappa's
        @rtype: tuple of 2 numpy arrays
        """
        return self.mus, self.kappas

    def par2str(self, precision=2):
        s=""
        for i in range(0, self.node_size): 
            mu=self.mus[i]
            k=self.kappas[i]
            s=s+"Mu %i: %f\nKappa %i: %f\n\n" % (i, mu, i, k)
        return s[:-2]
