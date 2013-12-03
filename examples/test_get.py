## FragBuild 0.1 Copyright (c) 2012 Anders S. Christensen
## Examples
import fragbuild
import string
import math

# Generate a alanine-valine peptide with methyl caps
my_peptide = fragbuild.peptide("V")

# Print some about the molecule
# print my_peptide.get_smiles() # Does not really work
print "Charge:  ",  my_peptide.get_charge()
print "Length:  ", my_peptide.get_length()
print "Sequence:", my_peptide.get_sequence()

# Set [phi, psi] angles of residue #2 to 180 degrees
# Note: Takes list in the format of either [phi, psi] or [phi, psi, omega]
my_peptide.set_residue_bb_angles(1, [math.pi, math.pi])

# Set [chi] angles of residue #2 to 63 degrees. 
# Note: Takes a list of angles as argument
chi1 = fragbuild.to_rad(63) # Convert 63 degrees to radians
angles = [chi1]
my_peptide.set_residue_chi_angles(1, angles)

# Setup optimization procedure
my_peptide.set_solvent("")
my_peptide.set_memory("3900mb")
my_peptide.set_nprocs(1)
my_peptide.set_method_and_basis_set("pm6")

# Write and execute optimization
my_peptide.write_g09_opt_file("test_opt.com")
fragbuild.run_g09("test_opt.com")

# Load optimized structure and prepare nmr_file
nmr_file = fragbuild.load_opt_file("test_opt.log")

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
nmr_file.set_method_and_basis_set("hf/sto-3g")

# Write and execute optimization
nmr_file.write_g09_nmr_file("test_nmr.com")
fragbuild.run_g09("test_nmr.com")

# Load results
nmr_data = fragbuild.nmr_file("test_nmr.log")

# Print all CA shielding values
print nmr_data.get_chem_shifts("CA")

