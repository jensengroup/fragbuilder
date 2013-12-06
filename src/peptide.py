import os
import openbabel
import uuid

import basilisk_lib

from residues import *
from math_utils import *

class peptide:

    def _init_residues(self, sequence):

        residues = [LeftMethylCap()]

        for i in sequence:
            residues.append(aa_dictionary[i]())

        residues.append(RightMethylCap())

        return residues

    def _get_residue_length(self, residue):
        i = 0
        for atom in residue.Mol:
            i += 1
        return i

    def _rotate_xyz(self, mol, axis, center, angle):
        for atom in mol:

                #Translate Atom, so Center corresponds to (0,0,0).
                atom[1] = atom[1] - center
                #Rotate!
                atom[1] = rotate(atom[1], axis - center, angle)
                atom[1] = atom[1] + center

        return mol


    def _get_all_bb_angles(self):

        angles = []
        for res_nr, residue in enumerate(self._residues):

            if res_nr == 0 or res_nr == len(self._residues) - 1:
                continue

            offset_prev = 0
            offset_this = 0
            offset_next = 0
            for j, residue2 in enumerate(self._residues):
                atoms_in_residue = self._get_residue_length(residue2)
                if j < res_nr - 1:
                    offset_prev += atoms_in_residue
                if j < res_nr:
                    offset_this += atoms_in_residue
                if j < res_nr + 1:
                    offset_next += atoms_in_residue

            omega = False
            phi = False
            psi = False

            # Special case for first residue with 
            if res_nr == 1 and self._state_nterm != "methyl":
                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

                omega = False
                phi = False
                psi = self._molecule.OBMol.GetTorsion(NH1, CA1, CO1, NH2) #PSI

            # For standard residues
            else:
                CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
                CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
                NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
                CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
                CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
                NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

                omega = self._molecule.OBMol.GetTorsion(CA0, CO0, NH1, CA1) #Omega
                phi = self._molecule.OBMol.GetTorsion(CO0, NH1, CA1, CO1) #PHI
                psi = self._molecule.OBMol.GetTorsion(NH1, CA1, CO1, NH2) #PSI

            angles.append([res_nr, [phi, psi, omega]])

        return angles

    def print_angles(self):
        print self._get_all_bb_angles()

    def _get_all_sc_angles(self):
        return 0



    def __init__(self, init_sequence, nterm="methyl", cterm="methyl"):
        """ Create a fragbuilder peptide with a specific sequence
            Default backbone angles are (phi, psi, omega) = (-120, 140, 180).

            Arguments:
            init_sequence -- sequence in one letter format. E.g. ""GGG""
            for a triple glycine peptide.

            Keyword arguments:
            nterm -- Type of cap at the n-terminal. (Default "methyl")
            cterm -- Type of cap at the c-terminal. (Default "methyl")

            NOTE: Options for cap types are "methyl", "neutral" "charged"

        """


        self._sequence = init_sequence.upper()

        self._state_nterm = nterm
        self._state_cterm = cterm

        self._molecule = self._assemble_peptide(self._sequence)
        self._charge   = self._calc_charge()
        self._length   = len(init_sequence)
        self._residues = self._init_residues(init_sequence)

        self._dbn = None




    def _reassemble_peptide(self, sequence):

        all_bb_angles = self._get_all_bb_angles()
        self._molecule = self._assemble_peptide(self._sequence)


    def get_smiles(self):
        return self._molecule

    def write_xyz(self, filename):
        self._molecule.write("xyz", filename, overwrite=True)


    def get_length(self):
        return self._length

    def set_residue_bb_angles(self, res_nr, angles):


        if len(angles) == 3:
            phi   = angles[0] * DEG_TO_RAD
            psi   = angles[1] * DEG_TO_RAD
            omega = angles[2] * DEG_TO_RAD
        elif len(angles) == 2:
            phi   = angles[0]  * DEG_TO_RAD
            psi   = angles[1]  * DEG_TO_RAD
            omega = math.pi
        else: 
            print "ERROR: Angles must be of type: [phi, psi, omega] or [phi, psi]"
            sys.exit(1)

        if res_nr < 1 or res_nr > len(self._residues) - 2:
                        print "ERROR: Error in set_residue_bb_angles. User supplied index:", res_nr, "Allowed range: 1 -", len(self._residues) - 2
                        sys.exit(1)

        offset_prev = 0
        offset_this = 0
        offset_next = 0
#       print self._residues

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)
            #print i, atoms_in_residue, self._residues[i]
            if i < res_nr - 1:
                offset_prev += atoms_in_residue
            if i < res_nr:
                offset_this += atoms_in_residue
            if i < res_nr + 1:
                offset_next += atoms_in_residue
#       print offset_prev, offset_this, offset_next
#       print "lol", self._residues[res_nr].Filename
        if res_nr == 1 and self._state_nterm != "methyl":
            NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
            CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
            CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
            NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

            self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI
        else:
            CA0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-2] + offset_prev)
            CO0 = self._molecule.OBMol.GetAtom(self._residues[res_nr - 1].BB[-1] + offset_prev)
            NH1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[0] + offset_this)
            CA1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[1] + offset_this)
            CO1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].BB[2] + offset_this)
            NH2 = self._molecule.OBMol.GetAtom(self._residues[res_nr + 1].BB[0] + offset_next)

            self._molecule.OBMol.SetTorsion(CA0, CO0, NH1, CA1, omega) #Omega
            self._molecule.OBMol.SetTorsion(CO0, NH1, CA1, CO1, phi) #PHI
            # If c-term cap is not methyl, there is no NH2 atom
            if NH2 != None:
                self._molecule.OBMol.SetTorsion(NH1, CA1, CO1, NH2, psi) #PSI


    def set_residue_chi_angles(self, res_nr, angles_deg):
        """ Sets the side chain torsion angles of a residue.

            Arguments:
            res_nr -- The index for the residue.
            angles_deg -- List of chi angles (in degrees).
        """

        angles = np.array(angles_deg) * DEG_TO_RAD

        offset = 0

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)
            if i < res_nr:
                offset += atoms_in_residue

        for i, angle in enumerate(angles):
            chi_atom1 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][0] + offset)
            chi_atom2 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][1] + offset)
            chi_atom3 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][2] + offset)
            chi_atom4 = self._molecule.OBMol.GetAtom(self._residues[res_nr].SC[i][3] + offset)

            self._molecule.OBMol.SetTorsion(chi_atom1, chi_atom2, chi_atom3, chi_atom4, angle)
        return

    _nprocs = 1
    _mem_in_mb = "400mb"
    _method_and_basis = "pm6"
    _comment = "No Title"
    _solvation = ""

    def set_solvent(self, solvent_type):
        if solvent_type == "":
            self._solvation = ""
        else:
            self._solvation = "scrf=(solvent=" + solvent_type + ")"

    def set_method_and_basis_set(self, command):
        self._method_and_basis = command


    def set_memory(self, amount):
        self._mem_in_mb = amount


    def set_nprocs(self, nprocs):
        self._nprocs = str(nprocs)


    def set_comment(self, comment):
        self._comment = comment


    def write_g09_opt_file(self, filename, constraint_dihedral=True):
        
        output_stream  = "%mem = "
        output_stream += self._mem_in_mb + "\n"
        output_stream += "%nprocs= " + str(self._nprocs)
        output_stream += "\n#t opt=modredundant " 
        output_stream += self._method_and_basis
        output_stream += " " + self._solvation
        output_stream += "\n\n" + self._comment+ "\n\n " 
        output_stream += str(self.get_charge()) + " 1\n" 

        for atom in self._molecule: 
            output_stream += atom.type[0] + " " + str(atom.coords[0]) + " " + str(atom.coords[1]) + " " + str(atom.coords[2]) + "\n"

        if constraint_dihedral:
            output_stream += self._get_freeze_string()
        output_stream += "\n\n"

        g09_file = open(filename, "w")
        g09_file.write(output_stream)
        g09_file.close()


    def get_sequence(self):
        return self._sequence


    def _calc_charge(self):

        charge = 0
        for i in self._sequence:
            charge += aa_dictionary[i].Charge
        return charge




    def get_charge(self):
        return self._charge


    # A list of dihedral angles to be frozen in G09 optimization.
    def _get_freeze_string(self):

        backbone_chain = []

        offset = 0

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)
            for atom in residue.BB:
                backbone_chain.append(atom + offset)
            offset += atoms_in_residue


        backbone_chain.sort()

        freeze_string = "\n"

        for i in range(len(backbone_chain)-3):
            freeze_string += "D " 
            for j in range(4):
                freeze_string += str(backbone_chain[i+j]) + " "
            freeze_string +="F\n"

        offset = 0

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)
            for dihedral in residue.SC:
                freeze_string += "D "
                for atom_index in dihedral:
                    freeze_string += str(atom_index + offset) + " "
                freeze_string +="F\n"
            offset += atoms_in_residue

        return freeze_string

    # Returns the peptide as an OBMol with phi/psi angles -120/140 degrees
    # which corresponds to nicely to a straight B-sheet strand.
    def _assemble_peptide(self, sequence):

        residues = [LeftMethylCap()]

        for i in sequence:
            residues.append(aa_dictionary[i]())

        residues.append(RightMethylCap())

        Fragment = residues

        TotalAtoms = 0
        FragmentCharge = 0
        for Residue in Fragment:
            TotalAtoms = TotalAtoms + len(Residue.Mol.atoms)
            FragmentCharge = FragmentCharge + Residue.Charge

        for i in range(len(Fragment)-1):
            #1.1 Get coordinate of connecting point.
            V3 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 1]-1][1]
            V2 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 2]-1][1]
            V1 = Fragment[i].AwesomeMol[Fragment[i].BB[len(Fragment[i].BB) - 3]-1][1]
            ConnectionPoint = V3 + (V2 - V1)/length(V2 - V1)*1.4
            ConnectionVector = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[0]-1][1] - ConnectionPoint

            #1.2 Translocate 
            for Atom in Fragment[i+1].AwesomeMol:
                Atom[1] = Atom[1] - ConnectionVector

            #2.1 Get rotation
            V4 = V3 - V2 + ConnectionPoint
            Axis1 = perp_vector(Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1], ConnectionPoint, V4 )
            Angle1 = - bondangle(Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1], ConnectionPoint, V4 )
            Center1 = ConnectionPoint

            #2.2 Rotate all coordinates in the AwesomeMol
            Fragment[i+1].AwesomeMol = self._rotate_xyz(Fragment[i+1].AwesomeMol, Axis1, Center1, Angle1)

            #2.3 Get other rotation around the dihedral
            Axis2 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1] - ConnectionPoint
            Axis2 = Axis2/length(Axis2) + ConnectionPoint
            D1 = V3
            D2 = ConnectionPoint
            D3 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[1]-1][1]
            D4 = Fragment[i+1].AwesomeMol[Fragment[i+1].BB[2]-1][1]
            # If next residue is proline, a little tweak is necessary.
            # using try/except, since there may be no next residue, lawl.
            try:
                if Fragment[i+1].ResName == "P":
                    Angle2 = math.pi - dihedral(D1, D2, D3, D4) + 90.0/180.0*math.pi
                else:
                    Angle2 = math.pi - dihedral(D1, D2, D3, D4)
            # If there is no next residue
            except:
                Angle2 = math.pi - dihedral(D1, D2, D3, D4) 

            Center2 = ConnectionPoint

            #2.4 Rotate all coordinates in the AwesomeMol
            Fragment[i+1].AwesomeMol = self._rotate_xyz(Fragment[i+1].AwesomeMol, Axis2, Center2, Angle2)

        # So far the peptide has been assembled with methyl-caps. Now I replace with charged/neutral caps
        # if the n-term and c-term states are not set to "methyl"

        # A little houskeeping here. Adjust total number of atoms and 
        if self._state_nterm == "methyl":
            temp = 0 # do nothing
        elif self._state_nterm == "neutral":
            TotalAtoms = TotalAtoms - len(Fragment[0].Mol.atoms) + 1
            Fragment[0] = LeftNeutralCap()

            NH = Fragment[1].AwesomeMol[Fragment[1].BB[0]-1][1]
            CA = Fragment[1].AwesomeMol[Fragment[1].BB[1]-1][1]
            CO = Fragment[1].AwesomeMol[Fragment[1].BB[2]-1][1]

            Fragment[0].AwesomeMol[0][1] = NH - ( (CO - CA)* 0.6  + perp_vector(NH, CA, CO) - CA ) * 0.8

        elif self._state_nterm == "charged":
            TotalAtoms = TotalAtoms - len(Fragment[0].Mol.atoms) + 2
            Fragment[0] = LeftChargedCap()

            NH = Fragment[1].AwesomeMol[Fragment[1].BB[0]-1][1]
            CA = Fragment[1].AwesomeMol[Fragment[1].BB[1]-1][1]
            CO = Fragment[1].AwesomeMol[Fragment[1].BB[2]-1][1]

            Fragment[0].AwesomeMol[0][1] = NH - ( (CO - CA)* 0.6  + perp_vector(NH, CA, CO) - CA ) * 0.8
            Fragment[0].AwesomeMol[1][1] = NH - ( (CO - CA)* 0.6  - perp_vector(NH, CA, CO) + CA ) * 0.8


        if self._state_cterm == "methyl":
            temp = 0 # do nothing

        elif self._state_cterm == "neutral":
            TotalAtoms = TotalAtoms - len(Fragment[-1].Mol.atoms) + 2
            Fragment[-1] = RightNeutralCap()

            NH = Fragment[-2].AwesomeMol[Fragment[-2].BB[0]-1][1]
            CA = Fragment[-2].AwesomeMol[Fragment[-2].BB[1]-1][1]
            CO = Fragment[-2].AwesomeMol[Fragment[-2].BB[2]-1][1]

            Fragment[-1].AwesomeMol[0][1] = CO - (NH - CA) * 0.8
            Fragment[-1].AwesomeMol[1][1] = CO + (CO - NH) * 0.8



        elif self._state_cterm == "charged":
            TotalAtoms = TotalAtoms - len(Fragment[-1].Mol.atoms) + 1
            Fragment[-1] = RightChargedCap()

            NH = Fragment[-2].AwesomeMol[Fragment[-2].BB[0]-1][1]
            CA = Fragment[-2].AwesomeMol[Fragment[-2].BB[1]-1][1]
            CO = Fragment[-2].AwesomeMol[Fragment[-2].BB[2]-1][1]

            Fragment[-1].AwesomeMol[0][1] = CO - (NH - CA) * 0.8

    
        uid = uuid.uuid4()
        temp_xyz = uid.hex  + ".xyz"
        file_out = open(temp_xyz, "w")
        file_out.write(str(TotalAtoms) + "\n")
        file_out.write("\n")
    
#       print Fragment
        for Residue in Fragment:
                for Atom in Residue.AwesomeMol:
                        file_out.write(Atom[0] + "  " +  str(Atom[1][0]) + "  " + str(Atom[1][1]) + "  " + str(Atom[1][2]) + "\n")
        file_out.write("\n")
        file_out.close()

        mol = pybel.readfile("xyz", temp_xyz).next()
        for k in range(3):
            for i in range(len(Fragment)-2):

                ThisResidue = i + 1
                NextOffset = 0
                for j in range(ThisResidue+1):
                    NextOffset = NextOffset + len(Fragment[j].Mol.atoms)
                Offset = 0
                for j in range(ThisResidue):
                    Offset = Offset + len(Fragment[j].Mol.atoms)
                PrevOffset = 0
                for j in range(ThisResidue-1):
                    PrevOffset = PrevOffset + len(Fragment[j].Mol.atoms)
                #try:
                CA0 = mol.OBMol.GetAtom(Fragment[ThisResidue - 1].BB[-2] + PrevOffset)
                CO0 = mol.OBMol.GetAtom(Fragment[ThisResidue - 1].BB[-1] + PrevOffset)
                NH1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[0] + Offset)
                CA1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[1] + Offset)
                CO1 = mol.OBMol.GetAtom(Fragment[ThisResidue].BB[2] + Offset)
                NH2 = mol.OBMol.GetAtom(Fragment[ThisResidue + 1].BB[0] + NextOffset )

                mol.OBMol.SetTorsion(CA0, CO0, NH1, CA1, math.pi) #Omega
                mol.OBMol.SetTorsion(CO0, NH1, CA1, CO1, -120.0/180.0*math.pi)
                mol.OBMol.SetTorsion(NH1, CA1, CO1, NH2, 140.0/180.0*math.pi)
                #print mol.OBMol.GetTorsion(CA0, CO0, NH1, CA1)
                #print mol.OBMol.GetTorsion(CO0, NH1, CA1, CO1)
                #print mol.OBMol.GetTorsion(NH1, CA1, CO1, NH2)
                #except:
                #    temp = 0

        os.remove(temp_xyz)

        return mol



    def get_residue_nr(self, atom_nr):
        """ Returns the residue number which has the atom of index atom_nr.

            Arguments:
            atom_nr -- Index of the atom in question.

        """

        correct_residue_nr = 0
        nr_of_atoms_up_to_this_residues = 0

        for i, residue in enumerate(self._residues):
            nr_of_atoms_up_to_this_residues += self._get_residue_length(residue)
            if atom_nr <= nr_of_atoms_up_to_this_residues:
                return i
        print "ERROR: Atom nr", atom_nr, "not found.", nr_of_atoms_up_to_this_residues," in molecule."
        exit(1)

    def _get_dihedral_constraints(self):

        restraints = []

        backbone_chain = []

        offset = 0

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)
            for atom in residue.BB:
                backbone_chain.append(atom + offset)
            offset += atoms_in_residue


        for i in range(len(backbone_chain)-3):

            local_restraint = []

            for j in range(4):
                 local_restraint.append(backbone_chain[i+j])
            restraints.append(local_restraint)

        offset = 0

        for i, residue in enumerate(self._residues):
            atoms_in_residue = self._get_residue_length(residue)

            for dihedral in residue.SC:
                local_restraint = []

                for atom_index in dihedral:
                    local_restraint.append(atom_index + offset)
                restraints.append(local_restraint)
            offset += atoms_in_residue


        return restraints

    def regularize(self, iterations=10, opt_steps = 20):
        """ Optimize the structure with fixed torsion angles.

            Keyword arguments:
            iterations -- Number of opimization+angle reset cycles. (default 10)
            opt_steps -- Number of optimization steps each cycle. (default 20)

        """
        restraints = self._get_dihedral_constraints()
        a = []

        for r in restraints:
            angle = self._molecule.OBMol.GetTorsion(r[0], r[1], r[2], r[3])
            #print "A", angle, r
            a.append(angle)

        for i in range(iterations):
            self.optimize(constraint=True, steps= opt_steps)
            for j, r in enumerate(restraints):
                self._molecule.OBMol.SetTorsion(r[0], r[1], r[2], r[3], a[j]/180.0*numpy.pi)
                #print self._molecule.OBMol.GetTorsion(r[0], r[1], r[2], r[3]), a[j]
        return




    def optimize(self, constraint=True, steps=100):
        """ Perform conjugate gradient optimization with the MMFF94 force field.

            Keyword arguments:
            constraint -- Apply harmonic constraint to dihedral angles. (default True)
            steps -- Max number of CG steps. (default 100)

        """

        mol = self._molecule.OBMol

        restraints = self.get_dihedral_constraints()

        cnstr = openbabel.OBFFConstraints()

        if constraint:
            for r in restraints:
                a = mol.GetTorsion(r[0], r[1], r[2], r[3])
                cnstr.AddTorsionConstraint(r[0], r[1], r[2], r[3], a)

        FF = openbabel.OBForceField.FindForceField("MMFF94")
        FF.Setup(mol, cnstr)
        FF.SetConstraints(cnstr)
        FF.ConjugateGradients(steps)
        FF.GetCoordinates(mol)

        return


    def sample_chi_angles(self, resnum):
        """ Resamples chi angles for a residue. 
            Returns the new chi angles.

            Arguments:
            resnum -- number of the residue to resample

        """

        if self._dbn is None:
            dbn = basilisk_lib.basilisk_dbn()


        angles = self._get_all_bb_angles()[resnum - 1]
        phi = angles[1][0]/180.0*numpy.pi
        psi = angles[1][1]/180.0*numpy.pi

        aa_letter = self._sequence[resnum - 1]

        aa = int(basilisk_lib.basilisk_utils.one_to_index(aa_letter))
        chis, bb, ll = dbn.get_sample(aa, phi, psi)

        return chis

    def energy(self):
        """ Returns the MMFF94 energy of the peptide in kcal/mol.

        """
        mol = self._molecule.OBMol
        FF = openbabel.OBForceField.FindForceField("MMFF94")
        FF.Setup(mol)
        return FF.Energy()


