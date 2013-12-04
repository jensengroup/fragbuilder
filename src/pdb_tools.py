import math

import Bio.PDB
from Bio.PDB import calc_dihedral

from names import *



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


def get_first_chain(filename):
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", filename)

    for model in structure:
        for chain in model:
            return chain

def read_phi_psi(chain):
    polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)

    for poly_index, poly in enumerate(polypeptides):
        return poly.get_phi_psi_list()




def read_angles(filename):
    """ Returns a dictionary with the following format:

        resnum is the number of a residue

        d[resnum]['resname']    : Returns the residue name
        d[resnum]['res_letter'] : Returns the one letter residue code
        d[resnum]['bb_angles']  : Returns a list containing [phi, psi]
        d[resnum]['chi_angles'] : Returns a list containing chi angles
    """

    chain = get_first_chain(filename)

    phipsi_angles = read_phi_psi(chain)

    first_residue_id = 0
    for residue in chain:
        if Bio.PDB.is_aa(residue):
            first_residue_id = residue.get_id()[1]
            break

    angles_dictionary = dict()

    for residue in chain:
        if not Bio.PDB.is_aa(residue):
            continue
        angles_dictionary[residue.get_id()[1]] = dict()

        angles_dictionary[residue.get_id()[1]]['res_name'] = residue.get_resname()
        angles_dictionary[residue.get_id()[1]]['res_letter'] = three_to_one(residue.get_resname())

        phi = phipsi_angles[residue.get_id()[1]-first_residue_id][0]
        if phi is not None:
            phi *= 180.0/math.pi
        psi = phipsi_angles[residue.get_id()[1]-first_residue_id][1]
        if psi is not None:
            psi *= 180.0/math.pi

        angles_dictionary[residue.get_id()[1]]['bb_angles'] = [phi, psi]  

        chi_angles = []
        for chi_angle in get_chi(residue):
            chi_angles.append(chi_angle/math.pi*180.0)

        angles_dictionary[residue.get_id()[1]]['chi_angles'] = chi_angles

    return angles_dictionary


