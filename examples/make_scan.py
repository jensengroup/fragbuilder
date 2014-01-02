from fragbuilder import Peptide

# Create a peptide object with the sequence
# glycine-leucine-glycine.
sequence = "GLG"
pep = Peptide(sequence)

# Define a list of angles.
angles = range(-180, 180, 60)

# Double loop over phi/psi angles in the desired range.
for psi in angles:
    for phi in angles:

        # Set the second (leucine) residue to (phi, psi).
        pep.set_bb_angles(2, [phi, psi])

        # Save each conformation.
        pep.write_xyz("pep_%04i_%04i.xyz" % (phi, psi))




