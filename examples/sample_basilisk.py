from fragbuilder import Basilisk_DBN

dbn = Basilisk_DBN()

chi, bb, ll = dbn.get_sample("K")

print "Chi angles:      ", chi
print "Phi/Psi angles:  ", bb
print "Log likelihood:  ", ll
