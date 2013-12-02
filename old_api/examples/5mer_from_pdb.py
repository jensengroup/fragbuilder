import fragbuild
import openbabel

# filename = "1CEX_babel.pdb"
filename = "5PTI.pdb"
output_dir = "5PTI_5mers"

residues_dictionary = fragbuild.read_angles(filename)

for res_num, data in residues_dictionary.items():

	if res_num+1 in residues_dictionary and res_num+2 in residues_dictionary and res_num+3 in residues_dictionary and res_num+4 in residues_dictionary:

		print res_num, data
		print res_num+1, residues_dictionary[res_num+1]
		print res_num+2, residues_dictionary[res_num+2]
		print res_num+3, residues_dictionary[res_num+3]
		print res_num+4, residues_dictionary[res_num+4]

		sequence = fragbuild.three_to_one(residues_dictionary[res_num+0]['res_name'])  + fragbuild.three_to_one(residues_dictionary[res_num+1]['res_name']) + fragbuild.three_to_one(residues_dictionary[res_num+2]['res_name']) + fragbuild.three_to_one(residues_dictionary[res_num+3]['res_name']) + fragbuild.three_to_one(residues_dictionary[res_num+4]['res_name'])


		my_peptide = fragbuild.peptide(sequence)


		for i, letter in enumerate(sequence):
	                angles = residues_dictionary[res_num+i]['bb_angles']
			angles_in_rad = [fragbuild.to_rad(angle) for angle in angles]
			my_peptide.set_residue_bb_angles(i+1, angles_in_rad)

		for i, letter in enumerate(sequence):
                        angles = residues_dictionary[res_num+i]['chi_angles']
                        angles_in_rad = [fragbuild.to_rad(angle) for angle in angles]
			if letter == 'P':
				continue
                        my_peptide.set_residue_chi_angles(i+1, angles_in_rad)

		my_peptide.write_g09_opt_file(output_dir + "/%03i" % (res_num + 1) + sequence + ".com")


print "Done!"
