from fragbuilder import Peptide
from fragbuilder import G09_opt
from fragbuilder import G09_energy
from fragbuilder import G09_NMR

sequence = "GLG"
pep = Peptide(sequence, nterm="methyl", cterm="methyl")

print "Charge:", pep.get_charge()

for i in pep.get_residue_numbers():
    pep.sample_bb_angles(i)
    pep.sample_chi_angles(i)

pep.optimize()

opt = G09_opt(pep)
opt.set_method("B3LYP/6-31G(d)")
opt.write_com("g09_opt_%s.com" % (sequence))
opt.write_com("g09_opt_%s_constraint.com" % (sequence), constraint_dihedrals=True)

energy = G09_energy(pep)
energy.set_method("PM3")
energy.write_com("g09_energy_%s.com" % (sequence))

nmr = G09_NMR(pep)
nmr.set_method("B3LYP/6-31+G(d,p)")
nmr.write_com("g09_nmr_%s.com" % (sequence))

