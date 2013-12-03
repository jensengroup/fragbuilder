import fragbuild

# Peptide sequence
sequence = "DAAAAAAAK"

# Initialize peptide
my_peptide = fragbuild.peptide(sequence)

# N- and C-terminii are methyl caps by default.
# Possible options are "neutral", "charged" and "methyl" (default).
#my_peptide.set_nterm("charged") 
#my_peptide.set_cterm("charged") 

# List with [phi, psi] backbone angles
angles_in_deg = [-37.0, -63.0]
angles_in_rad = [fragbuild.to_rad(angle) for angle in angles_in_deg]

# Set bb angles for all residues
for i in range(len(sequence)):
	my_peptide.set_residue_bb_angles(i+1, angles_in_rad)

# Write XYZ file
my_peptide.write_xyz("ahelix.xyz")

