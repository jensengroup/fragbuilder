from fragbuilder import Peptide

# Make a GLG peptide object
sequence = "GDdG"

pep = Peptide(sequence)


# Regularize
pep.regularize()

# Write files
pep.write_xyz("pep.xyz")



