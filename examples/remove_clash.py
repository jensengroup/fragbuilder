from fragbuilder import Peptide

# Make a peptide object
sequence = "YDG"
pep = Peptide(sequence)

# Set some angles that will clash
pep.set_bb_angles(1, [ -156.3, 142.3])
pep.set_bb_angles(2, [48.8, 42.3])
pep.set_bb_angles(3, [81.3, 0.2])

# Write an XYZ file
pep.write_xyz("clash.xyz")

# Remove clashes (keep dihedral angles)
pep.regularize()

# Write another XYZ file
pep.write_xyz("noclash.xyz")



