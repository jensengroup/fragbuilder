from fragbuilder import Peptide

# Make a GLG peptide object
sequence = "GLG"
pep = Peptide(sequence, nterm="neutral", cterm="charged")

print "SMILES:  ", pep.get_smiles()
print "Charge:  ", pep.get_charge()
print "Length:  ", pep.get_length()
print "Sequence:", pep.get_sequence()

# Set [PHI, PSI] angles for residues 1, 2 and 3
pep.set_bb_angles(1, [ -156.3,  142.3])
pep.set_bb_angles(2, [   48.8,   42.3])
pep.set_bb_angles(3, [   81.3,    0.2])

# Set and read chi angles for residue 2
pep.set_chi_angles(2, [  -50.6, -75.1])
print pep.get_bb_angles(2)

# Sample (and set) backbone angles for residue 2
print pep.sample_bb_angles(2)
print pep.get_bb_angles(2)

# Sample (and set) chi angles for residue 2
print pep.sample_chi_angles(2)
print pep.get_chi_angles(2)

# Write files
pep.write_xyz("pep.xyz")
pep.write_pdb("pep.pdb", QUIET=False)



