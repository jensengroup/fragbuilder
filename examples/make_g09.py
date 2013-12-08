import fragbuilder

sequence = "GLG"

pep = fragbuilder.peptide(sequence, nterm="methyl", cterm="methyl")

print "Charge:", pep.get_charge()

for i in pep.get_residue_numbers():
    if i != 1: continue
    pep.sample_residue_bb_angles(i)
    pep.sample_residue_chi_angles(i)


com = fragbuilder.G09_opt(pep)

com.set_method_and_basis_set("hf/sto-3g")
com.write_com("%s.com" % (sequence))


com.write_com("%s_constraint.com" % (sequence), constraint_dihedrals=True)
