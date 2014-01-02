from fragbuilder import PDB
from fragbuilder import Peptide

# Make a PDB object
pdb = PDB("1xnb.pdb")

# Extract length and sequence
n = pdb.get_length()
seq = pdb.get_sequence()

# Loop over residue numbers
for i in pdb.get_residue_numbers():

    # Ignore first and last few residues as omega,
    # phi, or psi may be undefined.
    if (i < 3) or (i > n - 2):
        continue

    # Determine sequence of the three-residues
    # peptide model and make a Peptide object
    pep_seq = seq[i-2:i+1]
    pep = Peptide(pep_seq)

    print i, pep_seq

    # Loop over preceding, central, and following residue.
    for j in [-1, 0, 1]:

        # Read backbone and chi angles from pdb object
        bb = pdb.get_residue_bb_angles(i + j)
        chi = pdb.get_residue_chi_angles(i + j)

        # Set matching angles in the Peptide object
        pep.set_bb_angles(2 + j, bb)
        pep.set_chi_angles(2 + j, chi)

    # Short optimization and save xyz file
    pep.regularize()
    pep.write_xyz("pep_%04i.xyz" % (i))
