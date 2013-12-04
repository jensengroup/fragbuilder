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
#
#	 load and interact with the dbn 
#

# import the shipped Mocapy
from Mocapy import *
import basilisk_parameter 

from InfEngineIOHMM import *
from LikelihoodIOHMM import *
from SamplerIOHMM import *

# 
# turning the chi angle + the amino acid
# information into an index 
#
# note that 0 and 1 are reserved for the phi and psi
# backbone angles, which are not amino acid dependent
aa_angle2index = { 
 '1-x1' :  2, 
 '2-x1' :  3,  '2-x2' :  4, 
 '3-x1' :  5,  '3-x2' :  6,  '3-x3' : 7,  
 '4-x1' :  8,  '4-x2' :  9,
 '6-x1' : 10,  '6-x2' : 11, 
 '7-x1' : 12,  '7-x2' : 13,
 '8-x1' : 14,  '8-x2' : 15,  '8-x3' : 16, '8-x4' : 17,
 '9-x1' : 18,  '9-x2' : 19, 
'10-x1' : 20, '10-x2' : 21, '10-x3' : 22,
'11-x1' : 23, '11-x2' : 24, 
'12-x1' : 25, '12-x2' : 26, 
'13-x1' : 27, '13-x2' : 28, '13-x3' : 29,
'14-x1' : 30, '14-x2' : 31, '14-x3' : 32, '14-x4' : 33, 
'15-x1' : 34,  
'16-x1' : 35, 
'17-x1' : 36, 
'18-x1' : 37, '18-x2' : 38,
'19-x1' : 39, '19-x2' : 40
}

#
# numbers of chi angles we expect for 
# different amino acid types
aa_chi_len = [0,1,2,3,2
			 ,0,2,2,4,2
			 ,3,2,2,3,4
			 ,1,1,1,2,2]


class basilisk_dbn:

	################################################################################
	def __init__(self, load_dbn=True):
		"""
		constructor, initialzing a couple values
		"""
		self.is_initialized = False
		self._dbn = None
		# lets see if we shall load the dbn right away
		if (load_dbn) :
			self._dbn = self.load_dbn_from_parameters()
			self.is_initialized = True

	################################################################################	
	def load_dbn_from_parameters(self) :
		"""
		loads a basilisk dbn from a parameter file	
		The basilisk DBN has the following network graph
		
			I		I		I		I
			|		|		|		|
			H  ---	H  ---  H  ---  H
			| 		|		|		|
		    VM		VM		VM		VM
		
		where I is an discrete index node, identifieing the
		angle to be modeled in the according slice, H is the 
		hidden with 30 states and VM is the vonMises node 
		that does actually model the angles in continous space.
		
		"""		 
		#
		# initialising all the nodes we need
		# note the parameters are in the cpds, the mus and kappas. 
		index=DiscreteNode(node_size=41, name='index', user_cpd=basilisk_parameter.index_cpd)

		h0=DiscreteNode(node_size=30, name='h0', user_cpd=basilisk_parameter.h0_cpd)
		h1=DiscreteNode(node_size=30, name='h1', user_cpd=basilisk_parameter.h1_cpd)

		angle = VMNode(node_size=30, name='angle', user_mus=basilisk_parameter.angle_mus, user_kappas=basilisk_parameter.angle_kappas)

		# assemble the nodes in the first and the following slices
		start_nodes=[index,h0, angle]
		end_nodes=[index,h1, angle]
		
		# set up a new dbn 
		dbn=DBN(start_nodes, end_nodes)
		
		# add the edges to the graph 
		dbn.add_intra('h1', 'angle')
		dbn.add_intra('index', 'h1')
		dbn.add_inter('h0', 'h1')
		
		dbn.construct()
		
		return dbn;

	################################################################################	
	def get_log_likelihood(self, residue_type, chis, phi=-5., psi=-5.) :
		"""
		Method to determine  the likelihood of a given set of angels.
		"""
		if not self.is_initialized :
			sys.stderr.write("ERROR: DBN has not been initialzed properly. Cannot determine likelihood.\n")
			return None
		
		# there is no sidechains for alanine and glycine
		if residue_type == 0 or residue_type == 5:
			return 0
		
		data, mism = self.get_data_mism_likelihood(residue_type, chis, phi, psi)
		ll_fw = LikelihoodIOHMM_Fw (self._dbn, 1)
		ll = ll_fw.calc_ll(data, mism, ignore_child_mism=False, set_data_in_dbn=True)[0]
		
		#
		# substract the backbone contribution 
		if phi > -5. and psi > -5. :
			for i in range(2, len(mism)) :
				mism[i][2] = 1;
			bb_ll = ll_fw.calc_ll(data, mism, ignore_child_mism=False, set_data_in_dbn=False)[0]
			ll-=bb_ll
			
		return ll
		
	################################################################################	
	def get_sample(self, residue_type, phi=-5., psi=-5., no_ll=False) :
		"""
		draws a sample from the dbn 
		input is the amino acid type and possibly the
		backbone angles
		"""
		if not self.is_initialized :
			sys.stderr.write("ERROR: DBN has not been initialzed properly. Cannot draw a sample.\n")
			return None
		
		# 
		data, mism = self.get_data_mism_sampling(residue_type, phi, psi)
		chis = []
		bb = []
		
		# setting up the sampler
		sampler = SamplerFwBt(self._dbn, 1)
		inf_engine = InfEngineIOHMM(self._dbn, sampler, data, mism, support_undo=False)
		sg = inf_engine.get_sample_generator()
		
		seq, slices = sg.next()
		for i in range(2, slices) :
			chis.append(seq[i][2])
		
		bb.append(seq[0][2])
		bb.append(seq[1][2])

		ll = -1.
		if not no_ll :
			ll = self.get_log_likelihood(residue_type, chis, bb[0], bb[1])

		return chis, bb, ll
	

		
	################################################################################		
	def get_data_mism_sampling(self, residue_type, phi=-5., psi=-5.) :
		"""
		assembles a data and mismask array for sampling. 
		optionally setting the phi and psi angles
		"""
		data = []
		mism = []
		
		if (phi > -4.) :
			data.append( [0,0,phi] ) 		#phi
			mism.append( [0,1,0] )
		else :
			data.append( [0,0,0.] ) 
			mism.append( [0,1,1] )	
		
		if (psi > -4.) :
			data.append( [1,0,psi] ) 		#phi
			mism.append( [0,1,0] )
		else :
			data.append( [1,0,0.] ) 
			mism.append( [0,1,1] )	
			
		n_x = aa_chi_len[residue_type]
		# append the rest of the indexes
		for i in range(1, n_x+1) :			# one indexed loop 
			index = aa_angle2index["%d-x%d" %(residue_type,i)]
			data.append( [index,0,0.] )		# chi angles
			mism.append( [0,1,1] )
		
		data = array(data)
		mism = array(mism)
		
		return data, mism

	################################################################################		
	def get_data_mism_likelihood(self, residue_type, chis, phi=-5., psi=-5.) :
		"""
		assembles a data and mismask array for the likelihood calculation. 
		optionally setting the phi and psi angles
		"""
		data = []
		mism = []
		
		if (phi > -4.) :
			data.append( [0,0,phi] ) 		#phi
			mism.append( [0,1,0] )
		else :
			data.append( [0,0,0.] ) 
			mism.append( [0,1,1] )	
		
		if (psi > -4.) :
			data.append( [1,0,psi] ) 		#phi
			mism.append( [0,1,0] )
		else :
			data.append( [1,0,0.] ) 
			mism.append( [0,1,1] )	
			
		n_x = aa_chi_len[residue_type]
		
		assert(n_x == len(chis)) 
		
		# append the rest of the indexes
		for i in range(1, n_x+1) :			# one indexed loop 
			index = aa_angle2index["%d-x%d" %(residue_type,i)]
			data.append( [index,0,chis[i-1]] )		# chi angles
			mism.append( [0,1,0] )
		
		data = array(data)
		mism = array(mism)
		
		return data, mism
		
		
		

