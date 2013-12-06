import fragbuilder
import string
import math

sequence = "GLG"

my_peptide = fragbuilder.peptide(sequence, nterm="methyl", cterm="methyl")

print my_peptide.get_smiles()
print "Charge:  ", my_peptide.get_charge()
print "Length:  ", my_peptide.get_length()
print "Sequence:", my_peptide.get_sequence()

my_peptide.set_residue_bb_angles(1, [ -156.3,  142.3])
my_peptide.set_residue_bb_angles(2, [   48.8,   42.3])
my_peptide.set_residue_bb_angles(3, [   81.3,    0.2])

print my_peptide.sample_residue_bb_angles(2)
print my_peptide.get_residue_bb_angles(2)

print my_peptide.sample_residue_chi_angles(2)
print my_peptide.get_residue_chi_angles(2)

my_peptide.set_residue_chi_angles(2, [  -50.6, -75.1])
print my_peptide.get_residue_chi_angles(2)

print my_peptide.get_bb_angles()
my_peptide.write_xyz("D101.xyz")



