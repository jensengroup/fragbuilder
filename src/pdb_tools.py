def get_chi(residue):
        if residue.get_resname() == 'ALA':
                return []
        if residue.get_resname() == 'GLY':
                return []
        if residue.get_resname() == 'ARG':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['NE'].get_vector()
                sc_atom7 = residue['CZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'ASN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ASP':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['OD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'CYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['SG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'GLU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'GLN':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['OE1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'HIS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD2'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'ILE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LEU':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'LYS':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                sc_atom7 = residue['NZ'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom6, sc_atom7)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'MET':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['SD'].get_vector()
                sc_atom6 = residue['CE'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom6)
                return [chi1, chi2, chi3]
        if residue.get_resname() == 'PHE':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'PRO':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom1)
                chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom1, sc_atom2)
                return [chi1, chi2, chi3, chi4]
        if residue.get_resname() == 'SER':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'THR':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['OG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        if residue.get_resname() == 'TRP':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'TYR':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG'].get_vector()
                sc_atom5 = residue['CD1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                return [chi1, chi2]
        if residue.get_resname() == 'VAL':
                sc_atom1 = residue['N'].get_vector()
                sc_atom2 = residue['CA'].get_vector()
                sc_atom3 = residue['CB'].get_vector()
                sc_atom4 = residue['CG1'].get_vector()
                chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                return [chi1]
        else:
                return "FAILLLL"


def read_angles(filename):

	parser=PDBParser(QUIET=True)
	structure=parser.get_structure("LOL", filename)
	dbn = TorusDBN(pdb_filename=filename, input_aa=True, input_angles=True)
	phipsi_angles = dbn.get_angles()
	
	#os.system("cp " + filename + " tmp.pdb")
	#os.system("sh /home/andersx/programs/fragbuilder/dssp_torsome.sh 2>/dev/null")
	#dssp_file = os.popen("cat tmp.sec")
	#dssp = dssp_file.readlines()[0]
	#dssp_file.close()
	
	first_residue_id = 0
	for model in structure[0]:
		for residue in model:
			if is_aa(residue):
				first_residue_id = residue.get_id()[1]
				break
		break
	
	angles_dictionary = dict()
	last_phi = 0
	last_psi = 0
	for model in structure[0]:
		for residue in model:
			if is_aa(residue):
				angles_dictionary[residue.get_id()[1]] = dict()

				angles_dictionary[residue.get_id()[1]]['res_name'] = residue.get_resname()
				angles_dictionary[residue.get_id()[1]]['res_letter'] = three_to_one(residue.get_resname())
				#angles_dictionary[residue.get_id()[1]]['DSSP'] = dssp[residue.get_id()[1]-1-first_residue_id]

				phi = phipsi_angles[residue.get_id()[1]-first_residue_id][0]*180/math.pi
				psi = phipsi_angles[residue.get_id()[1]-first_residue_id][1]*180/math.pi

                                angles_dictionary[residue.get_id()[1]]['bb_angles'] = [phi, psi]  

				chi_angles = []
				for chi_angle in get_chi(residue):
					chi_angles.append(chi_angle/math.pi*180.0)

				angles_dictionary[residue.get_id()[1]]['chi_angles'] = chi_angles

	return angles_dictionary


	
