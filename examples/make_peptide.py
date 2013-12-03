import fragbuilder
import string
import math

sequence = "YDG"

my_peptide = fragbuilder.peptide(sequence)

print my_peptide.get_smiles()
print "Charge:  ",  my_peptide.get_charge()
print "Length:  ", my_peptide.get_length()
print "Sequence:", my_peptide.get_sequence()


#ASPD11 angles

phi_psi1    = [fragbuilder.to_rad( -156.3),
               fragbuilder.to_rad( 142.3)]

phi_psi2    = [fragbuilder.to_rad(48.8),
               fragbuilder.to_rad(42.3)]

phi_psi3    = [fragbuilder.to_rad(81.3),
               fragbuilder.to_rad(0.2)]



my_peptide.set_residue_bb_angles(1, phi_psi1)
my_peptide.set_residue_bb_angles(2, phi_psi2)
my_peptide.set_residue_bb_angles(3, phi_psi3)


chi_angles = [fragbuilder.to_rad(-50.6),
              fragbuilder.to_rad(-75.1)]  #OBS TJEK
my_peptide.set_residue_chi_angles(2, chi_angles)

my_peptide.regularize()

my_peptide.write_xyz("D101.xyz")
my_peptide.write_g09_opt_file("my_opt_D101.com")





