import math
import sys

from math_utils import *
from names import *



class PDB:
    """ Usually instantiated from something like:

        pdbfile = fragbuilder.PDB("1UBQ.pdb")

        NOTE: This module requires Biopython to function. Will throw an error if
        Biopython is not installed. See http://www.biopython.org for download and help.
    """


    def __init__(self, pdbfile):
        """ Wrapper class for Bio.PDB which makes it convenient to
            read phi/psi/omega/chi torsion angles from a PDB-file.

            Arguments:
            pdbfile -- The PDB file you wish to read.

            WARNING: Will throw an error if Bio.PDB cannot be imported.
        """
        # Workaround for not having the entire fragbuilder module require
        # Biopython.
        try:
            import Bio.PDB
        except ImportError :
            sys.stderr.write("The PDB class requires Biopython to run. Please make sure you have biopython installed.\n See http://www.biopython.org for download and help.\n\n")
            sys.exit()

        from Bio.PDB import PDBParser
        #from Bio.PDB import PPBuilder

        try:
            self._parser = PDBParser(QUIET=True)
        except:
            # Workaround for missing QUIET keyword
            # in certain versions of Biopython.
            self._parser = PDBParser()

        self._pdbfile = pdbfile
        self._structure = self._parser.get_structure("pdb", self._pdbfile)
        self._chain = self._get_first_chain(self._structure)
        #self._polypeptide = PPBuilder().build_peptides(self._chain)

        #self._phipsi = None


    def get_length(self):
        """ Returns the length of the protein.
        """

        from Bio.PDB import is_aa

        length = 0
        for residue in self._chain:
            if is_aa(residue):
                length += 1

        return length

    def get_residue_numbers(self):
        """ Returns a list with indexes of all amino acids in the chain.
        
            Can be used for iterating over residues, e.g.:

            >>> for i in pdbfile.get_residue_numbers():
            ...     print i, pdbfile.get_residue_bb_angles(i)
        """
        length = self.get_length()

        return range(1, length + 1)

    def get_residue_chi_angles(self, resnum):
        """ Returns a list of chi angles for residue #resnum.
        """
        angles_rad = self._get_chi(self._chain[resnum])

        angles_deg = [angle * RAD_TO_DEG for angle in angles_rad]

        return angles_deg


    def get_residue_bb_angles(self, resnum):
        """ Returns a list of [phi, psi, omega] angles for residue #resnum.
        """
        from Bio.PDB import calc_dihedral

        length = self.get_length()

        angles_deg = []

        if resnum == 1:

            res_1 = self._chain[resnum]
            res_2 = self._chain[resnum + 1]

            N1  = res_1['N' ].get_vector()
            CA1 = res_1['CA'].get_vector()
            C1  = res_1['C' ].get_vector()
            N2  = res_2['N' ].get_vector()

            phi   = None
            psi   = calc_dihedral(N1, CA1, C1, N2) * RAD_TO_DEG
            omega = None

            angles_deg = [phi, psi, omega]

        elif resnum == length:

            res_0 = self._chain[resnum - 1]
            res_1 = self._chain[resnum]

            CA0 = res_0['CA'].get_vector()
            C0  = res_0['C' ].get_vector()
            N1  = res_1['N' ].get_vector()
            CA1 = res_1['CA'].get_vector()
            C1  = res_1['C' ].get_vector()

            phi   = calc_dihedral(C0, N1, CA1, C1) * RAD_TO_DEG
            psi   = None
            omega = calc_dihedral(CA0, C0, N1, CA1) * RAD_TO_DEG

            angles_deg = [phi, psi, omega]

        else:

            res_0 = self._chain[resnum - 1]
            res_1 = self._chain[resnum]
            res_2 = self._chain[resnum + 1]

            CA0 = res_0['CA'].get_vector()
            C0  = res_0['C' ].get_vector()
            N1  = res_1['N' ].get_vector()
            CA1 = res_1['CA'].get_vector()
            C1  = res_1['C' ].get_vector()
            N2  = res_2['N' ].get_vector()

            phi   = calc_dihedral(C0, N1, CA1, C1) * RAD_TO_DEG
            psi   = calc_dihedral(N1, CA1, C1, N2) * RAD_TO_DEG
            omega = calc_dihedral(CA0, C0, N1, CA1) * RAD_TO_DEG

            angles_deg = [phi, psi, omega]

        return angles_deg


    def _get_chi(self, residue):
            from Bio.PDB import calc_dihedral
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
                    # sc_atom1 = residue['N'].get_vector()
                    # sc_atom2 = residue['CA'].get_vector()
                    # sc_atom3 = residue['CB'].get_vector()
                    # sc_atom4 = residue['CG'].get_vector()
                    # sc_atom5 = residue['CD'].get_vector()
                    # chi1 = calc_dihedral(sc_atom1, sc_atom2, sc_atom3, sc_atom4)
                    # chi2 = calc_dihedral(sc_atom2, sc_atom3, sc_atom4, sc_atom5)
                    # chi3 = calc_dihedral(sc_atom3, sc_atom4, sc_atom5, sc_atom1)
                    # chi4 = calc_dihedral(sc_atom4, sc_atom5, sc_atom1, sc_atom2)
                    # return [chi1, chi2, chi3, chi4]
                    return []
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


    def _get_first_chain(self, structure):
        for model in structure:
            for chain in model:
                return chain

    def get_sequence(self):

        return self._get_sequence_from_chain(self._chain)



    def _get_sequence_from_chain(self, chain):

        sequence = ""

        for residue in chain:
            sequence += three_to_one(residue.get_resname())

        return sequence

    # def _calc_phi_psi(self, chain):
    #     for poly_index, poly in enumerate(self._polypeptide):
    #         return poly.get_phi_psi_list()




    # def _read_angles(self, filename):
    #     """ Returns a dictionary with the following format:

    #         resnum is the number of a residue

    #         d[resnum]['resname']    : Returns the residue name
    #         d[resnum]['res_letter'] : Returns the one letter residue code
    #         d[resnum]['bb_angles']  : Returns a list containing [phi, psi]
    #         d[resnum]['chi_angles'] : Returns a list containing chi angles
    #     """

    #     from Bio.PDB import is_aa
    #     chain = get_first_chain(filename)

    #     phipsi_angles = read_phi_psi(chain)

    #     first_residue_id = 0
    #     for residue in chain:
    #         if is_aa(residue):
    #             first_residue_id = residue.get_id()[1]
    #             break

    #     angles_dictionary = dict()

    #     for residue in chain:
    #         if not is_aa(residue):
    #             continue
    #         angles_dictionary[residue.get_id()[1]] = dict()

    #         angles_dictionary[residue.get_id()[1]]['res_name'] = residue.get_resname()
    #         angles_dictionary[residue.get_id()[1]]['res_letter'] = three_to_one(residue.get_resname())

    #         phi = phipsi_angles[residue.get_id()[1]-first_residue_id][0]
    #         if phi is not None:
    #             phi *= RAD_TO_DEG
    #         psi = phipsi_angles[residue.get_id()[1]-first_residue_id][1]
    #         if psi is not None:
    #             psi *= RAD_TO_DEG

    #         angles_dictionary[residue.get_id()[1]]['bb_angles'] = [phi, psi]  

    #         chi_angles = []
    #         for chi_angle in get_chi(residue):
    #             chi_angles.append(chi_angle * RAD_TO_DEG)

    #         angles_dictionary[residue.get_id()[1]]['chi_angles'] = chi_angles

    #     return angles_dictionary


