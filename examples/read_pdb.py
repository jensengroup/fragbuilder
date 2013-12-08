import fragbuilder

pdbfile = fragbuilder.PDB("1xnb.pqr")

n = pdbfile.get_length()

print "Protein contains %i residues!" % (n)

print pdbfile.get_residue_bb_angles(3)
print pdbfile.get_residue_chi_angles(3)

print "List of backbone phi/psi/omega angles:"

for i in pdbfile.get_residue_numbers():
    print i, pdbfile.get_residue_bb_angles(i)

