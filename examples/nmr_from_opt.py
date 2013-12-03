## FragBuild 0.1 Copyright (c) 2012 Anders S. Christensen
## Examples
import fragbuild
import string
import math
import sys

opt_filename = sys.argv[1]

# Load optimized structure and prepare nmr_file
nmr_file = fragbuild.load_opt_file(opt_filename)

# Check convergence
convergence = nmr_file.check_convergence()
print "Optimization is converged:", convergence
if not convergence:
	print "Bailing out .. optimization not converged."
	exit(0)

# Setup NMR calculation
nmr_file.set_solvent("water")
nmr_file.set_memory("1900mb")
nmr_file.set_nprocs(2)
nmr_file.set_nmr_method("giao")
nmr_file.set_method_and_basis_set("b3lyp/6-31g(d)")

nmr_filename = opt_filename.rstrip(".log") + "_nmr.com"

# Write and execute optimization
nmr_file.write_g09_nmr_file(nmr_filename)
#fragbuild.run_g09(nmr_filename)

# Load results
#nmr_data = fragbuild.nmr_file(nmr_filename.rstrip(".com") + ".log")

# Print all CA shielding values
#print nmr_data.get_chem_shifts("CA")

