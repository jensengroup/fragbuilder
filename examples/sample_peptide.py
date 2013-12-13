from fragbuilder import Peptide

# Make a GLAG peptide object
sequence = "GLAG"
pep = Peptide(sequence)

# Loop over residues:
for i in pep.get_residue_numbers():

    pep.sample_bb_angles(i)
    pep.sample_chi_angles(i)

# Regularize to remove any steric clashes
pep.regularize()

# Write files
pep.write_xyz("pep.xyz")



