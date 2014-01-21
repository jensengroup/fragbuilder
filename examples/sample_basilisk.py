from fragbuilder import Basilisk_DBN
from fragbuilder import set_seed

set_seed(12)

dbn = Basilisk_DBN()

chi, bb, ll = dbn.get_sample("K")

print "Chi angles:      ", chi
print "Phi/Psi angles:  ", bb
print "Log likelihood:  ", ll
