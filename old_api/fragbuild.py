## FragBuild 0.1 Copyright (c) 2012 Anders S. Christensen
## API file, fragbuild.py
import math
import numpy
import string
import sys
import os
import openbabel
import pybel

from Bio.PDB import *
from BackboneDBN import TorusDBN
from basilisk_lib import *

alphabet            = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
lower_case_alphabet = alphabet.lower()
numbers             = "0123456789"

regular_chars = alphabet + lower_case_alphabet + numbers + "._-"

#to_rad              = math.pi/180.0
#to_degree           = 180.0/math.pi

direct_g09 = "/home/andersx/pbs/direct_g09.sh"
fragbuild_dir = "/home/andersx/programs/fragbuilder/"



smiles = dict([('HN', "[$([H]N(C=O)C)]"     ),
               ('NH', "[$([N](C=O)C)]"      ),
               ('XX', "[$([N])]"      ),
               ('CO', "[$([C](NC)=O)]"      ),
               ('OC', "[$([O]=C(NC)C)]"     ),
               ('CA', "[$([C](NC)C(=O))]"   ),
               ('HA', "[$([H]C(NC)C(=O))]"  ),
	       ('CB', "[$([C]C(NC)C(=O))]"  )])


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



def run_g09(filename):
	os.system(direct_g09 + " " + filename)


def to_rad(degrees):
	degrees = float(degrees)
	return degrees/180.0*math.pi


def read_xyz(filename):
	return pybel.readfile("xyz", filename).next()


def ReadCoordzAwesome(Filename):
        ThisAwesomeMol = []
        XYZFile = open(Filename,"r")
        for Line in XYZFile:
                if (len(string.split(Line)) == 4):
                        Type = string.split(Line)[0]
                        X = float(string.split(Line)[1])
                        Y = float(string.split(Line)[2])
                        Z = float(string.split(Line)[3])
                        Coords = numpy.array([X,Y,Z])
                        Element = [Type, Coords]
                        ThisAwesomeMol.append(Element)
#       print ThisAwesomeMol
        return ThisAwesomeMol

 
# Rotamers From Lowell et al.,
# "The Penultimate Rotamer Library",
# PROTEINS: Structure, Function, and Genetics 40:389 - 408 (2000)

class LeftMethylCap:
	Charge 		= 0
	Index 		= 23
	Filename	= fragbuild_dir + "residues/left.xyz"
	ResName 	= "X"
	BB		= [6, 3, 1]
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightMethylCap:
	Charge 		= 0
	Index		= 24
	Filename	= fragbuild_dir + "residues/right.xyz"
	ResName 	= "X"
	BB		= [1, 3, 4]
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class LeftNeutralCap:
	Charge 		= 0
	Index 		= 25
	Filename	= fragbuild_dir + "residues/leftneutral.xyz"
	ResName 	= "N"
	BB		= []
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightNeutralCap:
	Charge 		= 0
	Index 		= 26
	Filename	= fragbuild_dir + "residues/rightneutral.xyz"
	ResName 	= "N"
	BB		= [1]
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class LeftChargedCap:
	Charge 		= 1
	Index 		= 27
	Filename	= fragbuild_dir + "residues/leftcharged.xyz"
	ResName 	= "Z"
	BB		= []
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightChargedCap:
	Charge 		= -1
	Index 		= 28
	Filename	= fragbuild_dir + "residues/rightcharged.xyz"
	ResName 	= "Z"
	BB		= [1]
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)



class Alanine:
	Charge 		= 0
	Index 		= 1
	Filename	= fragbuild_dir + "residues/ala.xyz"
	ResName 	= "A"
	BB		= [1, 3, 6]
	SC		= []
	Rotamer		= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Arginine:
	Charge 		= 1
	Index		= 20
	Filename	= fragbuild_dir + "residues/arg.xyz"
	ResName 	= "R"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13],[5,9,13,16],[9,13,16,18]]
	Rotamer 	= [[62, 180, 65, 85],
			  [62, 180, 65, -175],
			  [62, 180, 180, 85],
			  [62, 180, 180, 180],
			  [62, 180, 180, -85],
			  [62, 180, -65, 175],
			  [62, 180, -65, -85],
			  [-177, 65, 65, 85],
			  [-177, 65, 65, -175],
			  [-177, 65, 180, 85],
			  [-177, 65, 180, 180],
			  [-177, 180, 65, 85],
			  [-177, 180, 65, -175],
			  [-177, 180, 65, -105],
			  [-177, 180, 180, 85],
			  [-177, 180, 180, 180],
			  [-177, 180, 180, -85],
			  [-177, 180, -65, 105],
			  [-177, 180, -65, 175],
			  [-177, 180, -65, -85],
			  [-67, 180, 65, 85],
			  [-67, 180, 65, -175],
			  [-67, 180, 65, -105],
			  [-67, 180, 180, 85],
			  [-67, 180, 180, 180],
			  [-67, 180, 180, -85],
			  [-67, 180, -65, 105],
			  [-67, 180, -65, 175],
			  [-67, -167, -65, -85],
			  [-62, -68, 180, 85],
			  [-62, -68, 180, 180],
			  [-62, -68, 180, -85],
			  [-62, -68, -65, 175],
			  [-62, -68, -65, -85]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Aspargine:
	Charge 		= 0
	Index		= 16
	Filename	= fragbuild_dir + "residues/asn.xyz"
	ResName 	= "N"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,11]]
	Rotamer 	= [[62,	-10],
			  [62,	30],
			  [-174, -20],
			  [-177, 30],
			  [-65,	-20],
			  [-65,	-75],
			  [-65,	120]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Aspartate:
	Charge 		= -1
	Index		= 14
	Filename	= fragbuild_dir + "residues/asp.xyz"
	ResName 	= "D"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,11]]
	Rotamer 	= [[62, -10],
			  [62, 30],
			  [-177, 0],
			  [-177, 65]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Cysteine:
	Charge 		= 0
	Index		= 8
	Filename	= fragbuild_dir + "residues/cys.xyz"
	ResName 	= "C"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9]]
	Rotamer 	= [[62],[-177],[-65]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glutamate:
	Charge 		= -1
	Index		= 15
	Filename	= fragbuild_dir + "residues/glu.xyz"
	ResName 	= "E"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13],[5,9,13,15]]
	Rotamer 	= [[62,  180, -20],
			  [70,   -80, 0],
			  [-177,  65, 10],
			  [-177, 180, 0],
			  [-80,  -50, -25],
			  [-65,   85, 0],
			  [-67,  180, -10],
			  [-65,  -65, -40]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glutamine:
	Charge 		= 0
	Index		= 17
	Filename	= fragbuild_dir + "residues/gln.xyz"
	ResName 	= "Q"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13],[5,9,13,14]]
	Rotamer 	= [[62, 180, 20],
			  [70,   -75, 0],
			  [-177,  65, -100],
			  [-177,  65, 60],
			  [-177, 180, 0],
			  [-65,	  85, 0],
			  [-67,	 180, -25],
			  [-65,	 -65, -40],
			  [-65,	 -65, 100]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glycine:
	Charge 		= 0
	Index		= 5
	Filename	= fragbuild_dir + "residues/gly.xyz"
	ResName 	= "G"
	BB		= [1, 3, 4]
	SC		= []
	Rotamer 	= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Histidine:
	Charge 		= 1
	Index		= 18
	Filename	= fragbuild_dir + "residues/hip.xyz"
	ResName 	= "H"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,12]]
	Rotamer 	= [[62, -75],
			  [62, 80],
			  [-177, -165],
			  [-177, -80],
			  [-177, 60],
			  [-65,	-70],
			  [-65,	165]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Isoleucine:
	Charge 		= 0
	Index		= 3
	Filename	= fragbuild_dir + "residues/ile.xyz"
	ResName 	= "I"
	BB		= [1, 3, 6]
	SC		= [[1, 3, 5, 9],[3, 5, 9, 16]]
	Rotamer 	= [[ 62, 100],
			  [  62, 170],
			  [-177,  66],
			  [-177, 165],
			  [ -65, 100],
			  [ -65, 170],
			  [ -57, -60]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Leucine:
	Charge 		= 0
	Index		= 4
	Filename	= fragbuild_dir + "residues/leu.xyz"
	ResName 	= "L"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13]]
	Rotamer 	= [[ 62,  80],
			  [-177,  65],
			  [-172, 145],
			  [ -85,  65],
			  [ -65, 175]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Lysine:
	Charge 		= 1
	Index		= 19
	Filename	= fragbuild_dir + "residues/lys.xyz"
	ResName 	= "K"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13],[5,9,13,16],[9,13,16,19]]
	Rotamer 	= [[62, 180, 68, 180],
			  [62, 180, 180, 65],
			  [62, 180, 180, 180],
			  [62, 180, 180, -65],
			  [62, 180, -68, 180],
			  [-177, 68, 180, 65],
			  [-177, 68, 180, 180],
			  [68, -65, 9, 10],
			  [-177, 180, 68, 65],
			  [-177, 180, 68, 180],
			  [-177, 180, 180, 65],
			  [-177, 180, 180, 180],
			  [-177, 180, 180, -65],
			  [-177, 180, -68, 180],
			  [-177, 180, -68, -65],
			  [-90, 68, 180, 180],
			  [-67, 180, 68, 65],
			  [-67, 180, 68, 180],
			  [-67, 180, 180, 65],
			  [-67, 180, 180, 180],
			  [-67, 180, 180, -65],
			  [-67, 180, -68, 180],
			  [-67, 180, -68, -65],
			  [-62, -68, 180, 65],
			  [-62, -68, 180, 180],
			  [-62, -68, 180, -65],
			  [-62, -68, -68, 180]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Methionine:
	Charge 		= 0
	Index		= 9
	Filename	= fragbuild_dir + "residues/met.xyz"
	ResName 	= "M"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,13],[5,9,13,14]]
	Rotamer 	= [[62,	180, 75],
			   [62,	180, -75],
			   [-177, 65, 75],
			   [-177, 65, 180],
			   [-177, 180, 75],
			   [-177, 180, 180],
			   [-177, 180, -75],
			   [-67, 180, 75],
			   [-67, 180, 180],
			   [-67, 180, -75],
			   [-65, -65, 103],
			   [-65, -65, 180],
			   [-65, -65, -70]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Phenylalanine:
	Charge 		= 0
	Index		= 11
	Filename	= fragbuild_dir + "residues/phe.xyz"
	ResName 	= "F"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,11]]
	Rotamer 	= [[62, 90],
			   [-177, 80],
			   [-65, -85],
			   [-65, -30]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Proline:
	Charge 		= 0
	Index		= 10
	Filename	= fragbuild_dir + "residues/pro.xyz"
	ResName 	= "P"
	BB		= [1, 3, 9]
	SC		= []
	Rotamer 	= []
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Serine:
	Charge 		= 0
	Index		= 6
	Filename	= fragbuild_dir + "residues/ser.xyz"
	ResName 	= "S"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9]]
	Rotamer 	= [[62],
			   [-177],
			   [-65]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Threonine:
	Charge 		= 0
	Index		= 7
	Filename	= fragbuild_dir + "residues/thr.xyz"
	ResName 	= "T"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9]]
	Rotamer 	= [[62],
			   [-175],
			   [-65]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Tryptophan:
	Charge 		= 0
	Index		=  13
	Filename	= fragbuild_dir + "residues/trp.xyz"
	ResName 	= "W"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,11]]
	Rotamer 	= [[  62,  -90],
			   [  62,   90],
			   [-177, -105],
			   [-177,   90],
			   [ -65,  -90],
			   [ -65,   -5],
			   [ -65,   95]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Tyrosine:
	Charge 		= 0
	Index		= 12
	Filename	= fragbuild_dir + "residues/tyr.xyz"
	ResName 	= "Y"
	BB		= [1, 3, 6]
	SC		= [[1,3,5,9],[3,5,9,11]]
	Rotamer 	= [[62,		90],
			   [-177, 80],
			   [-65, -85],
			   [-65, -30]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Valine:
	Charge 		= 0
	Index		= 2
	Filename	= fragbuild_dir + "residues/val.xyz"
	ResName 	= "V"
	BB		= [1, 3, 6]
	SC		= [[1, 3, 5, 8]]
	Rotamer 	= [[63],
		          [175],
			  [-60]]
	def __init__(self):
		self.Mol = pybel.readfile("xyz", self.Filename).next()
		self.AwesomeMol = ReadCoordzAwesome(self.Filename)

 
# Unit vector normal to two vectors
def normalvector(v1, v2):
        J = cross(v1, v2)
        return J/(sqrt(dot(J, J)))

# Vector perpendicular to (v1-v2) and (v3-v2) centered in v2.
def perp_vector(v1, v2, v3):
        J = numpy.cross(v2 - v3, v2 - v1) + v2
        J = (J - v2) /(math.sqrt(numpy.dot(J - v2, J - v2))) + v2
        return J


# Bond angle between three coordinates
def bondangle(a,b,c):
	# In case numpy.dot() returns larger than 1
	# and we cannot take acos() to that number
	acos_out_of_bound = 1.0
	v1 = a - b
	v2 = c - b
	v1 = v1 / math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
	v2 = v2 / math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
	dot_product = numpy.dot(v1,v2)

	if dot_product > acos_out_of_bound:
		dot_product = acos_out_of_bound 
	if dot_product < -1.0 * acos_out_of_bound:
		dot_product = -1.0 * acos_out_of_bound 

	return math.acos(dot_product)


# Distance between two vectors
def distance(a,b):
        v1 = a - b
        return math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


# Length of a vector
def length(v1):
        return math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


# Dihedral angle between three points in space
def dihedral(a,b,c,d):
        b1 = b - a
        b2 = c - b
        b3 = d - c
        rb2 =math.sqrt(numpy.dot(b2,b2))
        b2xb3 = numpy.cross(b2,b3)
        b1xb2 = numpy.cross(b1,b2)
        dihedral = math.atan2( numpy.dot(rb2*b1,b2xb3), numpy.dot(b1xb2,b2xb3) )
        return dihedral


# Rotate a set of vectors
def rotate(V, J, T):
        x = V[0]
        y = V[1]
        z = V[2]
        u = J[0]
        v = J[1]
        w = J[2]
        a = (u*(u*x + v*y + w*z) + (x * (v*v + w*w) - u *(v*y + w*z))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(-w*y + v*z)*math.sin(T))/(u*u + v*v + w*w)
        b = (v*(u*x + v*y + w*z) + (y * (u*u + w*w) - v *(u*x + w*z))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(w*x - u*z)*math.sin(T))/(u*u + v*v + w*w)
        c = (w*(u*x + v*y + w*z) + (z * (u*u + v*v) - w *(u*x + v*y))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(-v*x + u*y)*math.sin(T))/(u*u + v*v + w*w)
        return numpy.array([a, b, c])


# Dictionary of amino acids and their constructors
aa_dictionary = dict([('G', Glycine       ),
		      ('A', Alanine       ),
                      ('S', Serine        ),
                      ('T', Threonine     ),
                      ('C', Cysteine      ),
                      ('V', Valine        ),
                      ('L', Leucine       ),
                      ('I', Isoleucine    ),
                      ('M', Methionine    ),
                      ('P', Proline       ),
                      ('F', Phenylalanine ),
                      ('Y', Tyrosine      ),
                      ('W', Tryptophan    ),
                      ('D', Aspartate     ),
                      ('E', Glutamate     ),
                      ('N', Aspargine     ),
                      ('Q', Glutamine     ),
                      ('H', Histidine     ),
                      ('K', Lysine        ),
                      ('R', Arginine      )])


def initialize_bb_angles(residues):

	raise NotImplementedError

def RotateAwesomeMol(AwesomeMol, Axis, Center, Angle):
#       print AwesomeMol
#       print Axis
#       print Center
#       print Angle
        for Atom in AwesomeMol:
                #Translate Atom, so Center corresponds to (0,0,0).
                Atom[1] = Atom[1] - Center
                #Rotate!
                Atom[1] = rotate(Atom[1], Axis - Center, Angle)
                Atom[1] = Atom[1] + Center

        return AwesomeMol


def get_charge_sum(sequence):

	charge = 0
	for i in sequence:
		charge += aa_dictionary[i].Charge
	return charge


class peptide:

	bb_atoms  = []
	sc_atoms  = []

	bb_angles = []
	sc_angles = []


	def _init_residues(self, sequence):
	
	        residues = [LeftMethylCap()]
	
	        for i in sequence:
	                residues.append(aa_dictionary[i]())
	
	        residues.append(RightMethylCap())
	
		return residues
	
	def _get_residue_length(self, residue):
		i = 0
		for atom in residue.Mol:
			i += 1
		return i
	
	def _get_all_bb_angles(self):

		angles = []
		for res_nr, residue in enumerate(self._residues):
	
			if res_nr == 0 or res_nr == len(self._residues) - 1:
				continue

			offset_prev = 0
			offset_this = 0
			offset_next = 0
			for j, residue2 in enumerate(self._residues):
				atoms_in_residue = self._get_residue_length(residue2)
				if j < res_nr - 1:
					offset_prev += atoms_in_residue
				if j < res_nr:
					offset_this += atoms_in_residue
				if j < res_nr + 1:
					offset_next += atoms_in_residue

			omega = False
			phi = False
			psi = False

			# Special case for first residue with 
			if res_nr == 1 and self._state_nterm != "methyl":
		                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
		                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
		                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
		                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)
	
				omega = False
				phi = False
				psi = self._molecule.OBMol.GetTorsion(NH1, CA1, CO1, NH2) #PSI
	
			# For standard residues
			else:
		                CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
		                CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
		                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
		                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
		                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
		                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)
	
				omega = self._molecule.OBMol.GetTorsion(CA0, CO0, NH1, CA1) #Omega
				phi = self._molecule.OBMol.GetTorsion(CO0, NH1, CA1, CO1) #PHI
				psi = self._molecule.OBMol.GetTorsion(NH1, CA1, CO1, NH2) #PSI

			angles.append([res_nr, [phi, psi, omega]])

		return angles


	def _get_all_sc_angles(self):
		return 0		


	def __init__(self, init_sequence):
		self._sequence = init_sequence.upper()

		#default values for n- and c-term state
		self._state_nterm = "methyl"
		self._state_cterm = "methyl"

		self._molecule = self._assemble_peptide(self._sequence)
		self._charge   = get_charge_sum(self._sequence)
		self._length   = len(init_sequence)
		self._residues = self._init_residues(init_sequence)


#		print "DBG RESIDUES:", self._residues

		self._get_freeze_string()

	# Create a new peptide with identical backbone angles.
	# For use when mutating residues or terminals.
	def _reassemble_peptide(self, sequence):
		
		all_bb_angles = self._get_all_bb_angles()
#		print all_bb_angles
		self._molecule = self._assemble_peptide(self._sequence)
		
#		self._set_all_angles(all_angles)



	def get_smiles(self):
		return self._molecule

	def set_nterm(self, state):
		if state in ["methyl", "charged", "neutral"]:
			self._state_nterm = state
			if state == "methyl":
				self._residues[0] = LeftMethylCap()
			elif state == "charged":
				self._residues[0] = LeftChargedCap()
			elif state == "neutral":
				self._residues[0] = LeftNeutralCap()

			# self._molecule = self._reassemble_peptide(self._sequence)
			self._reassemble_peptide(self._sequence)

		else:
			print "ERROR: Unsupported n-term,", state

		
	
	def set_cterm(self, state):
		if state in ["methyl", "charged", "neutral"]:
			self._state_cterm = state
			if state == "methyl":
				self._residues[-1] = RightMethylCap()
			elif state == "charged":
				self._residues[-1] = RightChargedCap()
			elif state == "neutral":
				self._residues[-1] = RightNeutralCap()
			self._reassemble_peptide(self._sequence)
		else:
			print "ERROR: Unsupported c-term,", state

	def write_xyz(self, filename):
		self._molecule.write("xyz", filename, overwrite=True)


	def get_length(self):
		return self._length

	def set_residue_bb_angles(self, res_nr, angles):
		if len(angles) == 3:
			phi   = angles[0]
			psi   = angles[1]
			omega = angles[2]
		elif len(angles) == 2:
			phi   = angles[0]
			psi   = angles[1]
			omega = math.pi
		else: 
			print "ERROR: Angles must be of type: [phi, psi, omega] or [phi, psi]"
			sys.exit(1)

		if res_nr < 1 or res_nr > len(self._residues) - 2:
                        print "ERROR: Error in set_residue_bb_angles. User supplied index:", res_nr, "Allowed range: 1 -", len(self._residues) - 2
                        sys.exit(1)
		
		offset_prev = 0
		offset_this = 0
		offset_next = 0
#		print self._residues

		for i, residue in enumerate(self._residues):
			atoms_in_residue = self._get_residue_length(residue)
			#print i, atoms_in_residue, self._residues[i]
			if i < res_nr - 1:
				offset_prev += atoms_in_residue
			if i < res_nr:
				offset_this += atoms_in_residue
			if i < res_nr + 1:
				offset_next += atoms_in_residue
#		print offset_prev, offset_this, offset_next
#		print "lol", self._residues[res_nr].Filename
		if res_nr == 1 and self._state_nterm != "methyl":
	                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
	                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
	                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
	                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

			self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI
		else:
	                CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
	                CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
	                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
	                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
	                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
	                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

			self._molecule.OBMol.SetTorsion(CA0, CO0, NH1, CA1, omega) #Omega
			self._molecule.OBMol.SetTorsion(CO0, NH1, CA1, CO1, phi) #PHI
			# If c-term cap is not methyl, there is no NH2 atom
			if NH2 != None:
				self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI

#		print "It worked!"


	def set_residue_omega(self, res_nr, angle):

		# if self.
                CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
                CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
                # CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
                # NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

		self._molecule.OBMol.SetTorsion(CA0, CO0, NH1, CA1, omega) #Omega
		# self._molecule.OBMol.SetTorsion(CO0, NH1, CA1, CO1, phi) #PHI
		# self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI


	def set_residue_phi(self, res_nr, angle):
                # CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
                CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
                # NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

		# self._molecule.OBMol.SetTorsion(CA0, CO0, NH1, CA1, omega) #Omega
		self._molecule.OBMol.SetTorsion(CO0, NH1, CA1, CO1, phi) #PHI
		# self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI


	def set_residue_psi(self, res_nr, angle):
                # CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
                # CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

		# self._molecule.OBMol.SetTorsion(CA0, CO0, NH1, CA1, omega) #Omega
		# self._molecule.OBMol.SetTorsion(CO0, NH1, CA1, CO1, phi) #PHI
		self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI




	def set_residue_chi_angles(self, res_nr, angles):
		offset = 0

		for i, residue in enumerate(self._residues):
			atoms_in_residue = self._get_residue_length(residue)
			if i < res_nr:
				offset += atoms_in_residue

		for i, angle in enumerate(angles):
			chi_atom1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][0] + offset)
			chi_atom2 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][1] + offset)
			chi_atom3 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][2] + offset)
			chi_atom4 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][3] + offset)

			self._molecule.OBMol.SetTorsion(chi_atom1, chi_atom2, chi_atom3, chi_atom4, angle)


	def set_random_chi_angles(self, res_nr):
		raise NotImplementedError


	def set_random_bb_angles(self, res_nr):
		raise NotImplementedError


	_nprocs = 1
	_mem_in_mb = "400mb"
	_method_and_basis = "pm6"
	_comment = "No Title"
	_solvation = ""

	def set_solvent(self, solvent_type):
		if solvent_type == "":
			self._solvation = ""
		else:
			self._solvation = "scrf=(solvent=" + solvent_type + ")"

	def set_method_and_basis_set(self, command):
		self._method_and_basis = command


	def set_memory(self, amount):
		self._mem_in_mb = amount


	def set_nprocs(self, nprocs):
		self._nprocs = str(nprocs)


	def set_comment(self, comment):
		self._comment = comment


	def write_g09_opt_file(self, filename):
		output_stream  = "%mem = "
		output_stream += self._mem_in_mb + "\n"
		output_stream += "%nprocs= " + str(self._nprocs)
		output_stream += "\n#t opt=modredundant " 
		output_stream += self._method_and_basis
		output_stream += " " + self._solvation
		output_stream += "\n\n" + self._comment+ "\n\n " 
		output_stream += str(self.get_charge()) + " 1\n" 

		for atom in self._molecule: 
			output_stream += atom.type[0] + " " + str(atom.coords[0]) + " " + str(atom.coords[1]) + " " + str(atom.coords[2]) + "\n"

		output_stream += self._get_freeze_string()
		output_stream += "\n\n"

		g09_file = open(filename, "w")
		g09_file.write(output_stream)
		g09_file.close()


        def get_sequence(self):
		return self._sequence


	def get_charge(self):
		return self._charge


	# A list of dihedral angles to be frozen in G09 optimization.
	def _get_freeze_string(self):

		backbone_chain = []

		offset = 0	

		for i, residue in enumerate(self._residues):
			atoms_in_residue = self._get_residue_length(residue)
			for atom in residue.BB:
				backbone_chain.append(atom + offset)
			offset += atoms_in_residue

		backbone_chain.sort()

		freeze_string = "\n"

		for i in range(len(backbone_chain)-3):
	        	freeze_string += "D " 
			for j in range(4):
				freeze_string += str(backbone_chain[i+j]) + " "
			freeze_string +="F\n"

		offset = 0

                for i, residue in enumerate(self._residues):
                        atoms_in_residue = self._get_residue_length(residue)
                        for dihedral in residue.SC:
				freeze_string += "D "       
				for atom_index in dihedral:
					freeze_string += str(atom_index + offset) + " "
				freeze_string +="F\n"
                        offset += atoms_in_residue

		
		return freeze_string


	#Backbone angles
		return 0

	def print_sequence(self):
		print self._sequence

	
	# Returns the peptide as an OBMol with phi/psi angles -120/140 degrees
	# which corresponds to nicely to a straight B-sheet strand.
	def _assemble_peptide(self, sequence):
	
		residues = [LeftMethylCap()]
	
		for i in sequence:
			residues.append(aa_dictionary[i]())
	
		residues.append(RightMethylCap())
	
		Fragment = residues
	
	        TotalAtoms = 0
	        FragmentCharge = 0
	        for Residue in Fragment:
	                TotalAtoms = TotalAtoms + len(Residue.Mol.atoms)
	                FragmentCharge = FragmentCharge + Residue.Charge
	
	        for i in range(len(Fragment)-1):
	                #1.1 Get coordinate of connecting point.
	                V3 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 1]-1][1]
	                V2 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 2]-1][1]
	                V1 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 3]-1][1]
	                ConnectionPoint = V3 + (V2 - V1)/length(V2 - V1)*1.4
	                ConnectionVector = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[0]-1][1] - ConnectionPoint
	
	                #1.2 Translocate 
	                for Atom in Fragment[i+1].AwesomeMol:
	                        Atom[1] = Atom[1] - ConnectionVector
	
	                #2.1 Get rotation
	                V4 = V3 - V2 + ConnectionPoint
	                Axis1 = perp_vector(Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1], ConnectionPoint, V4 )
	                Angle1 = - bondangle(Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1], ConnectionPoint, V4 )
	                Center1 = ConnectionPoint
	
	                #2.2 Rotate all coordinates in the AwesomeMol
	                Fragment[i+1].AwesomeMol = RotateAwesomeMol(Fragment[i+1].AwesomeMol, Axis1, Center1, Angle1)
	
	                #2.3 Get other rotation around the dihedral
	                Axis2 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1] - ConnectionPoint
	                Axis2 = Axis2/length(Axis2) + ConnectionPoint
	                D1 = V3
	                D2 = ConnectionPoint
	                D3 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1]
	                D4 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[2]-1][1]
			# If next residue is proline, a little tweak is necessary.
			# using try/except, since there may be no next residue, lawl.
			try:
				if Fragment[i+1].ResName == "P":
			                Angle2 = math.pi - dihedral(D1, D2, D3, D4) + 90.0/180.0*math.pi
				else:
			                Angle2 = math.pi - dihedral(D1, D2, D3, D4)
			# If there is no next residue
			except:
				Angle2 = math.pi - dihedral(D1, D2, D3, D4)	
	
	                Center2 = ConnectionPoint
	
	                #2.4 Rotate all coordinates in the AwesomeMol
	                Fragment[i+1].AwesomeMol = RotateAwesomeMol(Fragment[i+1].AwesomeMol, Axis2, Center2, Angle2)

		# So far the peptide has been assembled with methyl-caps. Now I replace with charged/neutral caps
		# if the n-term and c-term states are not set to "methyl"

		# A little houskeeping here. Adjust total number of atoms and 
		if self._state_nterm == "methyl":
			temp = 0 # do nothing
		elif self._state_nterm == "neutral":
			TotalAtoms = TotalAtoms - len(Fragment[0].Mol.atoms) + 1
			Fragment[0] = LeftNeutralCap()

			NH = Fragment[1].AwesomeMol[Fragment[1].BB[0]-1][1]
			CA = Fragment[1].AwesomeMol[Fragment[1].BB[1]-1][1]
			CO = Fragment[1].AwesomeMol[Fragment[1].BB[2]-1][1]

			Fragment[0].AwesomeMol[0][1] = NH - ( (CO - CA)* 0.6  + perp_vector(NH, CA, CO) - CA ) * 0.8
			
		elif self._state_nterm == "charged":
			TotalAtoms = TotalAtoms - len(Fragment[0].Mol.atoms) + 2
			Fragment[0] = LeftChargedCap()

			NH = Fragment[1].AwesomeMol[Fragment[1].BB[0]-1][1]
			CA = Fragment[1].AwesomeMol[Fragment[1].BB[1]-1][1]
			CO = Fragment[1].AwesomeMol[Fragment[1].BB[2]-1][1]

			Fragment[0].AwesomeMol[0][1] = NH - ( (CO - CA)* 0.6  + perp_vector(NH, CA, CO) - CA ) * 0.8
			Fragment[0].AwesomeMol[1][1] = NH - ( (CO - CA)* 0.6  - perp_vector(NH, CA, CO) + CA ) * 0.8


		if self._state_cterm == "methyl":
			temp = 0 # do nothing

		elif self._state_cterm == "neutral":
			TotalAtoms = TotalAtoms - len(Fragment[-1].Mol.atoms) + 2
			Fragment[-1] = RightNeutralCap()

			NH = Fragment[-2].AwesomeMol[Fragment[-2].BB[0]-1][1]
			CA = Fragment[-2].AwesomeMol[Fragment[-2].BB[1]-1][1]
			CO = Fragment[-2].AwesomeMol[Fragment[-2].BB[2]-1][1]

			Fragment[-1].AwesomeMol[0][1] = CO - (NH - CA) * 0.8
			Fragment[-1].AwesomeMol[1][1] = CO + (CO - NH) * 0.8



		elif self._state_cterm == "charged":
			TotalAtoms = TotalAtoms - len(Fragment[-1].Mol.atoms) + 1
			Fragment[-1] = RightChargedCap()

			NH = Fragment[-2].AwesomeMol[Fragment[-2].BB[0]-1][1]
			CA = Fragment[-2].AwesomeMol[Fragment[-2].BB[1]-1][1]
			CO = Fragment[-2].AwesomeMol[Fragment[-2].BB[2]-1][1]

			Fragment[-1].AwesomeMol[0][1] = CO - (NH - CA) * 0.8

	
		import uuid
		uid = uuid.uuid4()
		temp_xyz = uid.hex  + ".xyz"
	        file_out = open(temp_xyz, "w")
	        file_out.write(str(TotalAtoms) + "\n")
	        file_out.write("\n")
	
#		print Fragment
	        for Residue in Fragment:
	                for Atom in Residue.AwesomeMol:
	                        file_out.write(Atom[0] + "  " +  str(Atom[1][0]) + "  " + str(Atom[1][1]) + "  " + str(Atom[1][2]) + "\n")
	        file_out.write("\n")
	        file_out.close()

	        mol = pybel.readfile("xyz", temp_xyz).next()
	        for i in range(len(Fragment)-2):

	                ThisResidue = i + 1
	                NextOffset = 0
	                for j in range(ThisResidue+1):
	                        NextOffset = NextOffset + len(Fragment[j].Mol.atoms)
	                Offset = 0
	                for j in range(ThisResidue):
	                        Offset = Offset + len(Fragment[j].Mol.atoms)
	                PrevOffset = 0
	                for j in range(ThisResidue-1):
	                        PrevOffset = PrevOffset + len(Fragment[j].Mol.atoms)
			try:	
		                CA0 = mol.OBMol.GetAtom(Fragment[ThisResidue - 1].BB[-2] + PrevOffset)
		                CO0 = mol.OBMol.GetAtom(Fragment[ThisResidue - 1].BB[-1] + PrevOffset)
		                NH1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[0] + Offset)
		                CA1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[1] + Offset)
		                CO1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[2] + Offset)
		                NH2 = mol.OBMol.GetAtom(Fragment[ThisResidue + 1].BB[0] + NextOffset )
	
		                mol.OBMol.SetTorsion(CA0, CO0, NH1, CA1, math.pi) #Omega
		                mol.OBMol.SetTorsion(CO0, NH1, CA1, CO1, -120.0/180.0*math.pi)
		                mol.OBMol.SetTorsion(NH1, CA1, CO1, NH2, 140.0/180.0*math.pi)
			except:
				temp = 0
#				print self._state_nterm
		
		os.remove(temp_xyz)
#		print Fragment
	
		return mol

#self._hbonds = []


	def get_residue_nr(atom_nr):
		correct_residue_nr = 0
		nr_of_atoms_up_to_this_residues = 0

		for i, residue in enumerate(self._residues):
			nr_of_atoms_up_to_this_residues += self._get_residue_length(residue)
			if atom_nr <= nr_of_atoms_up_to_this_residues:
				return 
		print "ERROR: Atom nr", atom_nr, "not found.", nr_of_atoms_up_to_this_residues," in molecule."
		exit(1)

	def print_molecule(self):
		print self._molecule	


	def add_hbond(self, res_no, atom_type, bonding_partner):
#		if atom_type == "HA":
#			
#		elif atom_type == "HN":
#		elif atom_type == "OC"


		raise NotImplementedError


	def set_hbond_geom(r_oh, theta, rho):
		
		raise NotImplementedError





class nmr_file:


	def __init__(self, filename):
		self._filename = filename

		self._atom_type = ""

	        # Casper Steinmann's boilerplate code   
	        self._obmol = openbabel.OBMol()
	        self._obpat = openbabel.OBSmartsPattern()
	        self._obconv = openbabel.OBConversion()
	        self._obconv.SetInFormat("g09")
	        self._obconv.ReadFile(self._obmol, self._filename)


	def get_matches(self, atom_type):
		atom = atom_type[0]
		pattern = smiles[atom_type]
		self._obpat.Init(pattern)
		self._obpat.Match(self._obmol)
		print pattern
		matches = [m[0] for m in self._obpat.GetUMapList()]
		print matches
		return matches


	def get_chem_shifts(self, atom_type):
		matches = self.get_matches(atom_type)

		matches_shieldings = []		

		for i in matches:
			command = "grep \"" + str(i) + "  " + atom_type[0] + "    Isotropic \" " + self._filename
			output = os.popen(command)
			shielding = string.split(output.readlines()[0][0:-1])[4]
			data = i, float(shielding)

			matches_shieldings.append(data)

		return matches_shieldings


class load_opt_file:

        def __init__(self, filename):

                self._filename  = filename
		self._converged = self.check_convergence
		self._molecule  = pybel.readfile("g09", self._filename).next()

		self._method_and_basis 	= "b3lyp/6-311+g(2d,p)"
		# Don't use solvent as default
		# self._solvation	= "scrf=(solvent=water)"
		self._solvation		= ""
		self._nmr_method	= "giao"
		self._memory		= "1900mb"
		self._nprocs		= "1"


	def check_convergence(self):

		tests = []

		command = "grep -c \"Normal termination of Gaussian 09 at\" " + self._filename
		output = os.popen(command).readlines()

		if int(output[0]) == 1:
			tests.append(True)
		else:
			tests.append(False)

		command = "grep -c \"Stationary point found\" " + self._filename
		output = os.popen(command).readlines()

		if int(output[0]) == 1:
			tests.append(True)
		else:
			tests.append(False)

		if False in tests:
			return False
		else:
			return True


	def set_solvent(self, solvent_type):
		if solvent_type == "":
			self._solvation = ""
		else:
			self._solvation = "scrf=(solvent=" + solvent_type + ")"


	def set_method_and_basis_set(self, command):
		self._method_and_basis = command


	def set_memory(self, amount):
		self._memory = amount


	def set_nprocs(self, nprocs):
		self._nprocs = str(nprocs)


	def set_nmr_method(self, method):
		self._nmr_method = method


	def write_g09_nmr_file(self, filename):
		for i in filename:
			if i not in regular_chars:
				print "ERROR: Writing filename, but filename illegal character:", i
				exit(1)
		output = pybel.Outputfile("com", filename, overwrite=True)
		output.write(self._molecule)
		output.close()
		
		options     = "%mem=" + self._memory + "\\n"
		options    += "%nprocs=" + self._nprocs + "\\n"
		options    += "#p nmr=" + self._nmr_method + " " + self._method_and_basis + " " + self._solvation
		old_options = "#Put Keywords Here, check Charge and Multiplicity."



#		print "sed -i \"s@" + old_options + "@" + options + "@\" "+ filename
		os.system("sed -i \"s@" + old_options + "@" + options + "@\" "+ filename)

	def write_output(self, filename):
		write_g09_nmr_file(filename)
		run_g09(filename)


def get_chi(residue):
        if residue.get_resname() == 'ALA':
                return []
        if residue.get_resname() == 'GLY':
                return []
        if residue.get_resname() == 'ARG':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['NE'].get_vector()
                sc_atom7 = residue['CZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'ASN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ASP':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'CYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['SG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'GLU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'GLN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'HIS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD2'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ILE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LEU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                sc_atom7 = residue['NZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'MET':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['SD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'PHE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'PRO':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom1)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom1, sc_atom2)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'SER':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'THR':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'TRP':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'TYR':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'VAL':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        else:
                return "FAILLLL"


def read_angles(filename):

	parser=PDBParser(QUIET=True)
	structure=parser.get_structure("LOL", filename)
	dbn = TorusDBN(pdb_filename=filename, input_aa=True, input_angles=True)
	phipsi_angles = dbn.get_angles()
	
	os.system("cp " + filename + " tmp.pdb")
	os.system("sh /home/andersx/programs/fragbuilder/dssp_torsome.sh 2>/dev/null")
	dssp_file = os.popen("cat tmp.sec")
	dssp = dssp_file.readlines()[0]
	dssp_file.close()
	
	first_residue_id = 0
	for model in structure[0]:
		for residue in model:
			if is_aa(residue):
				first_residue_id = residue.get_id()[1]
				break
		break
	
	angles_dictionary = dict()
	last_phi = 0
	last_psi = 0
	for model in structure[0]:
		for residue in model:
			if is_aa(residue):
				angles_dictionary[residue.get_id()[1]] = dict()

				angles_dictionary[residue.get_id()[1]]['res_name'] = residue.get_resname()
				angles_dictionary[residue.get_id()[1]]['res_letter'] = three_to_one(residue.get_resname())
				angles_dictionary[residue.get_id()[1]]['DSSP'] = dssp[residue.get_id()[1]-1-first_residue_id]

				phi = phipsi_angles[residue.get_id()[1]-first_residue_id][0]*180/math.pi
				psi = phipsi_angles[residue.get_id()[1]-first_residue_id][1]*180/math.pi

                                angles_dictionary[residue.get_id()[1]]['bb_angles'] = [phi, psi]  

				chi_angles = []
				for chi_angle in get_chi(residue):
					chi_angles.append(chi_angle/math.pi*180.0)

				angles_dictionary[residue.get_id()[1]]['chi_angles'] = chi_angles

	return angles_dictionary


	
