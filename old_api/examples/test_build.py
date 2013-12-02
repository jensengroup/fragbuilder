import fragbuild
import string
import math
#my_peptide = fragbuild.peptide("ACDEFGHIKLMNPQRSTVWY")
sequence = "ACDEFGHIKLMNPQRSTVWY"

#sequence = "AR"

my_peptide = fragbuild.peptide(sequence)

print my_peptide.get_smiles()

print "Charge:  ",  my_peptide.get_charge()

print "Length:  ", my_peptide.get_length()

print "Sequence:", my_peptide.get_sequence()

my_peptide.set_residue_bb_angles(2, [math.pi, math.pi])

my_peptide.write_xyz("my_peptide.xyz")

my_peptide.write_g09_opt_file("my_opt.com")


#for i, angles in enumerate(fragbuild.Arginine.Rotamer):

#	angles_in_rad = [fragbuild.to_rad(angle) for angle in angles]

#	print angles_in_rad
#	my_peptide.set_residue_chi_angles(2, angles_in_rad)
#	my_peptide.write_xyz("arg_test/arg_" + str(i) + ".xyz")


