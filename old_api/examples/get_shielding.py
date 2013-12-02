import openbabel
import os
import string

filenames = os.popen("ls minimized_pcm/*.log")
#print filenames.readlines()

alpha_carbon = "[$(CC(=O)[N,O])]"
alpha_hydrogen = "[$([H]C(NC)C(=O))]"

atom_type = "C"

pattern = alpha_carbon

data = []	
for filename_with_newline in filenames.readlines():
	# Remove trailing newline character
	filename = filename_with_newline[0:-1]

	# Casper Steinmann's boilerplate code	
	obmol = openbabel.OBMol()
	obpat = openbabel.OBSmartsPattern()
	obconv = openbabel.OBConversion()
	obconv.SetInFormat("g09")
	obconv.ReadFile(obmol, filename)
	obpat.Init(pattern)
	obpat.Match(obmol)
	matches = [m[0] for m in obpat.GetUMapList()]

	# Use second alpha carbon	
	i = matches[1]

	# Grep for the match and the shielding constant	
	command = "grep '" + str(i) + "  " + atom_type + "    Isotropic ' " + filename
	output = os.popen(command)
	res_nr = int(string.split(filename[0:-4], "X")[2])
	shielding = float(string.split(output.readlines()[0])[4])	

	data.append([res_nr, shielding])

data = sorted(data)
for i in data:
	print i[0]+17, i[1]
















