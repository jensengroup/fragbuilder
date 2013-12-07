#  BASILISK: A generative, probabilistic model of side chains in proteins
#
#  Copyright (C) 2010   Tim Harder and Jes Frellsen
#
#  BASILISK is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BASILISK is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BASILISK.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################
#
#  a couple functions to make life easier.
#


import sys

#
# lets see if we can find biopython
try:
    import Bio.PDB
except ImportError :                    
    """
    FRAGBUILDER: DOES NOT STRICTLY REQUIRE BioPython
    """
    # sys.stderr.write("This program requires biopython to run. Please make sure you have biopython installed. \nSee http://www.biopython.org for download and help.\n\n")
    # sys.exit()

from os.path import isfile


############################################################################
#
# that should really be in python .. well .. here is the workaround
def is_numeric(obj) :
	return (is_float(obj)) 

def is_int(obj) :
	try :
		x = int(obj)
	except ValueError:
		return False
	else :
		# we need that, because python 
		# just rounds a float to an int
		if x == obj :
			return True
	return False

def is_float(obj) :
	try :
		x = float(obj)
	except ValueError:
		return False
	else :
		return True
	return False



############################################################################
#
#  a couple ways to transform aminoacid labels 
aa1="ACDEFGHIKLMNPQRSTVWY"
aa3=["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
     "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

d1_to_index={}
dindex_to_1={}
d3_to_index={}
dindex_to_3={}

for i in range(0, 20):
    n1=aa1[i]
    n3=aa3[i]
    d1_to_index[n1]=i
    dindex_to_1[i]=n1
    d3_to_index[n3]=i
    dindex_to_3[i]=n3

def index_to_one(index):
    "Amino acid index to single letter (eg 0 to A)"
    return dindex_to_1[index]

def one_to_index(s):
    "Amino acid single letter to index (eg A to 0)"
    return d1_to_index[s]

def index_to_three(i):
    "Amino acid index to three letter (eg 0 to ALA)"
    return dindex_to_3[i]

def three_to_index(s):
    "Amino acid three letter to index (eg ALA to 0)"
    return d3_to_index[s]

def three_to_one(s):
    "Amino acid three letter to single letter (eg ALA to A)"
    i=d3_to_index[s]
    return dindex_to_1[i]

def one_to_three(s):
    "Amino acid single letter to three letter (eg A to ALA)"
    i=d1_to_index[s]
    return dindex_to_3[i]


########################################                                      
def load_pdb_file(filename) :               
	"""                             
	Loads a PDB file as Bio.PDB structure object and returns it.
	"""                        
	
	if isfile( filename) :                                          	
		parser = Bio.PDB.PDBParser()                                
		structure = parser.get_structure(filename, filename)
	else :                                                            
		sys.stderr.write("No such file .. cannot load PDB. (%s)\n" % filename)
		sys.exit()
	return structure

