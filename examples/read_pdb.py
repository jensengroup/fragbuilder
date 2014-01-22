import fragbuilder
import sys

filename = sys.argv[1]

pdbfile = fragbuilder.PDB(filename)

n = pdbfile.get_length()

print "Protein contains %i residues!" % (n)

print pdbfile.get_bb_angles(3)
print pdbfile.get_chi_angles(3)


seq = pdbfile.get_sequence()

print "Sequence is:", seq

print "List of backbone phi/psi/omega angles:"

for i in pdbfile.get_residue_numbers():
    print i, pdbfile.get_resname(i),
    print pdbfile.get_bb_angles(i), pdbfile.get_chi_angles(i)

