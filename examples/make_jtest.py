import fragbuild
import string
import math

sequence = "A"

my_peptide = fragbuild.peptide(sequence)

print my_peptide.get_smiles()
print "Charge:  ",  my_peptide.get_charge()
print "Length:  ", my_peptide.get_length()
print "Sequence:", my_peptide.get_sequence()


psis = [170.0, 120.0, -60.0]
phi = -60.0

for i, psi in enumerate(psis):
    phi_psi    = [fragbuild.to_rad(phi),
                   fragbuild.to_rad(psi)]
    my_peptide.set_residue_bb_angles(1, phi_psi)
    my_peptide.write_g09_opt_file("jtests/ala_simple%i.com" % (i))





