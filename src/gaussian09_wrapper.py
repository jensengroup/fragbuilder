class nmr_file:


	def __init__(self, filename):
		self._filename = filename

		self._atom_type = ""

	        # Casper Steinmann's boilerplate code   
	        self._obmol = openbabel.OBMol()
	        self._obpat = openbabel.OBSmartsPattern()
	        self._obconv = openbabel.OBConversion()
	        self._obconv.SetInFormat("g09")
	        self._obconv.ReadFile(self._obmol, self._filename)


	def get_matches(self, atom_type):
		atom = atom_type[0]
		pattern = smiles[atom_type]
		self._obpat.Init(pattern)
		self._obpat.Match(self._obmol)
		print pattern
		matches = [m[0] for m in self._obpat.GetUMapList()]
		print matches
		return matches


	def get_chem_shifts(self, atom_type):
		matches = self.get_matches(atom_type)

		matches_shieldings = []		

		for i in matches:
			command = "grep \"" + str(i) + "  " + atom_type[0] + "    Isotropic \" " + self._filename
			output = os.popen(command)
			shielding = string.split(output.readlines()[0][0:-1])[4]
			data = i, float(shielding)

			matches_shieldings.append(data)

		return matches_shieldings


class load_opt_file:

        def __init__(self, filename):

                self._filename  = filename
		self._converged = self.check_convergence
		self._molecule  = pybel.readfile("g09", self._filename).next()

		self._method_and_basis 	= "b3lyp/6-311+g(2d,p)"
		# Don't use solvent as default
		# self._solvation	= "scrf=(solvent=water)"
		self._solvation		= ""
		self._nmr_method	= "giao"
		self._memory		= "1900mb"
		self._nprocs		= "1"


	def check_convergence(self):

		tests = []

		command = "grep -c \"Normal termination of Gaussian 09 at\" " + self._filename
		output = os.popen(command).readlines()

		if int(output[0]) == 1:
			tests.append(True)
		else:
			tests.append(False)

		command = "grep -c \"Stationary point found\" " + self._filename
		output = os.popen(command).readlines()

		if int(output[0]) == 1:
			tests.append(True)
		else:
			tests.append(False)

		if False in tests:
			return False
		else:
			return True


	def set_solvent(self, solvent_type):
		if solvent_type == "":
			self._solvation = ""
		else:
			self._solvation = "scrf=(solvent=" + solvent_type + ")"


	def set_method_and_basis_set(self, command):
		self._method_and_basis = command


	def set_memory(self, amount):
		self._memory = amount


	def set_nprocs(self, nprocs):
		self._nprocs = str(nprocs)


	def set_nmr_method(self, method):
		self._nmr_method = method


	def write_g09_nmr_file(self, filename):
		for i in filename:
			if i not in regular_chars:
				print "ERROR: Writing filename, but filename illegal character:", i
				exit(1)
		output = pybel.Outputfile("com", filename, overwrite=True)
		output.write(self._molecule)
		output.close()
		
		options     = "%mem=" + self._memory + "\\n"
		options    += "%nprocs=" + self._nprocs + "\\n"
		options    += "#p nmr=" + self._nmr_method + " " + self._method_and_basis + " " + self._solvation
		old_options = "#Put Keywords Here, check Charge and Multiplicity."



#		print "sed -i \"s@" + old_options + "@" + options + "@\" "+ filename
		os.system("sed -i \"s@" + old_options + "@" + options + "@\" "+ filename)

	def write_output(self, filename):
		write_g09_nmr_file(filename)
		run_g09(filename)


###  _nprocs = 1
###  _mem_in_mb = "400mb"
###  _method_and_basis = "pm6"
###  _comment = "No Title"
###  _solvation = ""

###  def set_solvent(self, solvent_type):
###      if solvent_type == "":
###          self._solvation = ""
###      else:
###          self._solvation = "scrf=(solvent=" + solvent_type + ")"

###  def set_method_and_basis_set(self, command):
###      self._method_and_basis = command


###  def set_memory(self, amount):
###      self._mem_in_mb = amount


###  def set_nprocs(self, nprocs):
###      self._nprocs = str(nprocs)


###  def set_comment(self, comment):
###      self._comment = comment


###  def write_g09_opt_file(self, filename, constraint_dihedral=True):
###
###      output_stream  = "%mem = "
###      output_stream += self._mem_in_mb + "\n"
###      output_stream += "%nprocs= " + str(self._nprocs)
###      output_stream += "\n#t opt=modredundant "
###      output_stream += self._method_and_basis
###      output_stream += " " + self._solvation
###      output_stream += "\n\n" + self._comment+ "\n\n "
###      output_stream += str(self.get_charge()) + " 1\n"

###      for atom in self._molecule:
###          output_stream += atom.type[0] + " " + str(atom.coords[0]) + " " + str(atom.coords[1]) + " " + str(atom.coords[2]) + "\n"

###      if constraint_dihedral:
###          output_stream += self._get_freeze_string()
###      output_stream += "\n\n"

###      g09_file = open(filename, "w")
###      g09_file.write(output_stream)
###      g09_file.close()






###    # A list of dihedral angles to be frozen in G09 optimization.
###    def _get_freeze_string(self):
###
###        backbone_chain = []
###
###        offset = 0
###
###        for i, residue in enumerate(self._residues):
###            atoms_in_residue = self._get_residue_length(residue)
###            for atom in residue.BB:
###                backbone_chain.append(atom + offset)
###            offset += atoms_in_residue
###
###
###        backbone_chain.sort()
###
###        freeze_string = "\n"
###
###        for i in range(len(backbone_chain)-3):
###            freeze_string += "D "
###            for j in range(4):
###                freeze_string += str(backbone_chain[i+j]) + " "
###            freeze_string +="F\n"
###
###        offset = 0
###
###        for i, residue in enumerate(self._residues):
###            atoms_in_residue = self._get_residue_length(residue)
###            for dihedral in residue.SC:
###                freeze_string += "D "
###                for atom_index in dihedral:
###                    freeze_string += str(atom_index + offset) + " "
###                freeze_string +="F\n"
###            offset += atoms_in_residue
###
###        return freeze_string



