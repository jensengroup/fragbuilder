import openbabel
filename = "my_peptide.xyz"

alpha_carbon = "[$(CC(=O)[N,O])]"

alpha_hydrogen = "[$([H]C(NC)C(=O))]"

pattern = alpha_hydrogen



obmol = openbabel.OBMol()
obpat = openbabel.OBSmartsPattern()
obconv = openbabel.OBConversion()
obconv.SetInFormat(filename[-3:])
obconv.ReadFile(obmol, filename)
obpat.Init(pattern)
obpat.Match(obmol)
matches = [m[0] for m in obpat.GetUMapList()]

print matches
