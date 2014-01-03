import sys


class G09_base:


    def __init__(self, pep):
        """ Class to generate input file for a Gaussian optimization.

            Arguments:
            peptide -- A an object of fragbuilder.Peptide type.
        """
        self._nprocs = 1
        self._mem_in_mb = "1gb"
        self._method_and_basis = "pm6"
        self._comment = "No Title"
        self._solvation = ""
        self._header = ""
        self._peptide = pep

    def set_solvent(self, solvent_type):
        """ Set the header of the input Gaussian 09 .com file.

            E.g:
            com.set_solvent("water")

            This inserts the keyword "scrf=(solvent=water)".
        """
        if solvent_type == "":
            self._solvation = ""
        else:
            self._solvation = "scrf=(solvent=" + solvent_type + ")"

    def set_method(self, command):
        """ Set the method and basis set in the input Gaussian 09 .com 
            file using standard Gaussian 09 notation.

            E.g.: 
            com.set_method("B3LYP/6-31G(d)")
        """
        self._method_and_basis = command


    def set_memory(self, amount):
        """ Set amount of memory of the input Gaussian 09 .com file.

            E.g:
            com.set_memory("2gb")
        """
        self._mem_in_mb = amount


    def set_nprocs(self, nprocs):
        """ Set the nprocs keyword of the input Gaussian 09 .com file.

            E.g:
            com.set_nprocs(8)
        """
        self._nprocs = str(nprocs)


    def set_comment(self, comment):
        """ Set the comment section of the input Gaussian 09 .com file.
        """
        self._comment = comment

    def set_header(self, header):
        """ Set the header of the put Gaussian 09 .com file.

            E.g:
            com.set_header("%mem=2gb\\n# opt B3LYP/6-31G(d)")

            WARNING: Using this option overrides all other header options!
        """
        self._header = header


class G09_opt(G09_base):
    """ Class to setup a Gaussian 09 optimization from a
        a peptide.

        Requires a peptide objec. A simple use case could be:

        from fragbuilder import Peptide, G09_opt
        pep = fragbuilder.Peptide("GGGG")
        opt = G09_opt(pep)
        opt.write_com()
    """
    def __init__(self, pep):
        """ Class to generate input file for a Gaussian optimization.

            Arguments:
            peptide -- A an object of fragbuilder.Peptide type.
        """
        self._nprocs = 1
        self._mem_in_mb = "1gb"
        self._method_and_basis = "pm6"
        self._comment = "No Title"
        self._solvation = ""
        self._header = ""
        self._peptide = pep


    def write_com(self, filename, constraint_dihedrals=False):
        """ Writes a Gaussian 09 optimization optimization .com file.

            Arguments:
            filename -- Name of the file  to be written.

            Keyword arguments:
            constraint_dihedrals -- Whether to constraint dihedral angles
            during geometry optimization. (Default False)


            NOTES: The default header is:

            %mem=1gb
            %nprocs=1
            #t opt pm6

            WARNING: Will overwrite an existing file with the same name.
        """

        if self._header == "":
            output_stream  = "%mem = "
            output_stream += self._mem_in_mb + "\n"
            output_stream += "%nprocs= " + str(self._nprocs)
            if constraint_dihedrals:
                output_stream += "\n#t opt=modredundant "
            else:
                output_stream += "\n#t opt "
            output_stream += self._method_and_basis
            output_stream += " " + self._solvation
        else:
            outpt_stream += self._header

        output_stream += "\n\n" + self._comment+ "\n\n "
        output_stream += str(self._peptide.get_charge()) + " 1\n"

        for atom in self._peptide._molecule:
            #output_stream += atom.type[0] + " " + str(atom.coords[0]) + " " + str(atom.coords[1]) + " " + str(atom.coords[2]) + "\n"
            output_stream += " %2s           %16.10f %16.10f %16.10f\n" % (atom.type[0],
                    atom.coords[0], atom.coords[1], atom.coords[2])

        if constraint_dihedrals:
            output_stream += self._get_freeze_string()
        output_stream += "\n\n"

        g09_file = open(filename, "w")
        g09_file.write(output_stream)
        g09_file.close()


    # A list of dihedral angles to be frozen in G09 optimization.
    def _get_freeze_string(self):
        """ Returns a list of dihedral constraints in G09 format.
            Reads from a peptide object.
        """

        backbone_chain = []

        offset = 0

        for i, residue in enumerate(self._peptide._residues):
            atoms_in_residue = self._peptide._get_residue_length(residue)
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

        for i, residue in enumerate(self._peptide._residues):
            atoms_in_residue = self._peptide._get_residue_length(residue)
            for dihedral in residue.SC:
                freeze_string += "D "
                for atom_index in dihedral:
                    freeze_string += str(atom_index + offset) + " "
                freeze_string +="F\n"
            offset += atoms_in_residue

        return freeze_string



class G09_energy(G09_base):
    """ Class to setup a Gaussian 09 single point energy 
        calculation from a peptide.

        Requires a peptide objec. A simple use case could be:

        from fragbuilder import Peptide, G09_energy
        pep = fragbuilder.Peptide("GGGG")
        opt = G09_energy(pep)
        opt.write_com()
    """

    def __init__(self, pep):
        """ Class to generate input file for a Gaussian optimization.

            Arguments:
            peptide -- A an object of fragbuilder.Peptide type.
        """
        self._nprocs = 1
        self._mem_in_mb = "1gb"
        self._method_and_basis = "pm6"
        self._comment = "No Title"
        self._solvation = ""
        self._header = ""
        self._peptide = pep


    def write_com(self, filename):
        """ Writes a Gaussian 09 optimization optimization .com file.

            Arguments:
            filename -- Name of the file  to be written.

            NOTES: The default header is:

            %mem=1gb
            %nprocs=1
            # pm6

            WARNING: Will overwrite an existing file with the same name.
        """

        if self._header == "":
            output_stream  = "%mem = "
            output_stream += self._mem_in_mb + "\n"
            output_stream += "%nprocs= " + str(self._nprocs) 
            output_stream += "\n# "
            output_stream += self._method_and_basis
            output_stream += " " + self._solvation
        else:
            outpt_stream += self._header

        output_stream += "\n\n" + self._comment+ "\n\n "
        output_stream += str(self._peptide.get_charge()) + " 1\n"

        for atom in self._peptide._molecule:
            output_stream += " %2s           %16.10f %16.10f %16.10f\n" % (atom.type[0],
                    atom.coords[0], atom.coords[1], atom.coords[2])
        output_stream += "\n\n"

        g09_file = open(filename, "w")
        g09_file.write(output_stream)
        g09_file.close()


class G09_NMR(G09_base):
    """ Class to setup a Gaussian 09 NMR chemical shielding 
        calculation from a peptide.

        Requires a peptide objec. A simple use case could be:

        from fragbuilder import Peptide, G09_NMR
        pep = fragbuilder.Peptide("GGGG")
        opt = G09_NMR(pep)
        opt.write_com()
    """

    def __init__(self, pep):
        """ Class to generate input file for a Gaussian NMR shielding calculation.

            Arguments:
            peptide -- A an object of fragbuilder.Peptide type.
        """
        self._nprocs = 1
        self._mem_in_mb = "1gb"
        self._method_and_basis = "b3lyp/6-311G+(2d,p)"
        self._comment = "No Title"
        self._solvation = ""
        self._header = ""
        self._peptide = pep


    def write_com(self, filename):
        """ Writes a Gaussian 09 optimization optimization .com file.

            Arguments:
            filename -- Name of the file  to be written.

            NOTES: The default header is:

            %mem=1gb
            %nprocs=1
            #p nmr=giao b3lyp/6-311G+(2d,p)

            WARNING: Will overwrite an existing file with the same name.
        """

        if self._header == "":
            output_stream  = "%mem = "
            output_stream += self._mem_in_mb + "\n"
            output_stream += "%nprocs= " + str(self._nprocs) 
            output_stream += "\n#p nmr=giao "
            output_stream += self._method_and_basis
            output_stream += " " + self._solvation
        else:
            outpt_stream += self._header

        output_stream += "\n\n" + self._comment+ "\n\n "
        output_stream += str(self._peptide.get_charge()) + " 1\n"

        for atom in self._peptide._molecule:
            output_stream += " %2s           %16.10f %16.10f %16.10f\n" % (atom.type[0],
                    atom.coords[0], atom.coords[1], atom.coords[2])
        output_stream += "\n\n"

        g09_file = open(filename, "w")
        g09_file.write(output_stream)
        g09_file.close()


