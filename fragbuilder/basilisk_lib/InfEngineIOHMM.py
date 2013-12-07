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

from Mocapy.AbstractInfEngine import AbstractInfEngine
from Mocapy.MocapyExceptions import MocapyException
from numpy import array


class InfEngineIOHMM(AbstractInfEngine):
    """
    Inference engines for inference in IO HMMs.
    """

    def __init__(self, dbn, sampler, seq, mismask, weight=1, support_undo=False):
        AbstractInfEngine.__init__(self, dbn, seq, weight)
        self.sampler=sampler
        self.mismask=mismask
        self.start=0
        self.end=len(seq)
        self.seq=seq

        self.support_undo=support_undo
        if support_undo:
            self.prev_seq=array(seq)
        else:
            self.prev_seq=None

        if weight!=1:
            raise MocapyException, "Weight different from 1 not supported"

        # List of Family objects - one for each node
        # The Family objects tie the data and the nodes together
        self.parentmap_list=self._make_parentmap_list(seq, weight)

    def get_sample_generator(self, mismask=None):
        """
        Return a sample generator (which is a Python generator object).
        The sample generator's 'next' method will return a sampled sequence. 
        """
        if mismask is None:
            mismask=self.mismask
        sampler=self.sampler

        # Making sure we start with the original node values
        self._set_parentmap(self.parentmap_list) 

        # Start the real thing
        while 1:
            if self.support_undo:
                self.prev_seq=array(self.seq)
            slice_count=sampler.sweep(mismask=mismask, start=self.start, end=self.end)
            sampled_seq=array(self.seq)
            yield (sampled_seq, slice_count)


    def set_start_end(self, start, end):
        """
 
        @param start: start position in the sequence 
        @type start: int

        @param end: end position in the sequence
        @type end: int
        """
        if not 0<=start<self.seq_len or not start<end<=self.seq_len:
            raise MocapyException, "Start:end out of bounds"
        self.start=start
        self.end=end

    def undo(self):
        """
        Undo last sampling step.
        """
        if self.prev_seq is None:
            raise MocapyException, "Undo functionality was not enabled"
        for pm in self.parentmap_list:
            pm.replace_seq(self.prev_seq)
        self.seq = self.prev_seq
        self.prev_seq = None

    def replace_seq(self, seq):
        if self.support_undo:
            self.prev_seq = self.seq
        for pm in self.parentmap_list:
            pm.replace_seq(seq)
        self.seq = seq

    def get_viterbi(self):
        raise NotImplementedError

