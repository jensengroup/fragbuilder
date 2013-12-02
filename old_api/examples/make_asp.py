import fragbuild
import string
import math

sequence = "ADG"

my_peptide = fragbuild.peptide(sequence)

print my_peptide.get_smiles()
print "Charge:  ",  my_peptide.get_charge()
print "Length:  ", my_peptide.get_length()
print "Sequence:", my_peptide.get_sequence()


#ASP119 angles

phi_psi1    = [fragbuild.to_rad( -58.8),
               fragbuild.to_rad( -33.1)]

phi_psi2    = [fragbuild.to_rad( -79.6),
               fragbuild.to_rad(  -0.6)]

phi_psi3    = [fragbuild.to_rad(  87.0),
               fragbuild.to_rad( 177.6)]



my_peptide.set_residue_bb_angles(1, phi_psi1)
my_peptide.set_residue_bb_angles(2, phi_psi2)
my_peptide.set_residue_bb_angles(3, phi_psi3)


chi_angles = [fragbuild.to_rad(-70.9),
              fragbuild.to_rad(-30.9)]
my_peptide.set_residue_chi_angles(2, chi_angles)


my_peptide.write_xyz("my_peptide.xyz")
my_peptide.write_g09_opt_file("my_opt.com")





