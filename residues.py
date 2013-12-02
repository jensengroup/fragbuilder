class LeftMethylCap:
    Charge      = 0
    Index       = 23
    Filename    = fragbuild_dir + "residues/left.xyz"
    ResName     = "X"
    BB      = [6, 3, 1]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightMethylCap:
    Charge      = 0
    Index       = 24
    Filename    = fragbuild_dir + "residues/right.xyz"
    ResName     = "X"
    BB      = [1, 3, 4]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class LeftNeutralCap:
    Charge      = 0
    Index       = 25
    Filename    = fragbuild_dir + "residues/leftneutral.xyz"
    ResName     = "N"
    BB      = []
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightNeutralCap:
    Charge      = 0
    Index       = 26
    Filename    = fragbuild_dir + "residues/rightneutral.xyz"
    ResName     = "N"
    BB      = [1]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class LeftChargedCap:
    Charge      = 1
    Index       = 27
    Filename    = fragbuild_dir + "residues/leftcharged.xyz"
    ResName     = "Z"
    BB      = []
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class RightChargedCap:
    Charge      = -1
    Index       = 28
    Filename    = fragbuild_dir + "residues/rightcharged.xyz"
    ResName     = "Z"
    BB      = [1]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)



class Alanine:
    Charge      = 0
    Index       = 1
    Filename    = fragbuild_dir + "residues/ala.xyz"
    ResName     = "A"
    BB      = [1, 3, 6]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Arginine:
    Charge      = 1
    Index       = 20
    Filename    = fragbuild_dir + "residues/arg.xyz"
    ResName     = "R"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13],[5,9,13,16],[9,13,16,18]]
    Rotamer     = [[62, 180, 65, 85],
              [62, 180, 65, -175],
              [62, 180, 180, 85],
              [62, 180, 180, 180],
              [62, 180, 180, -85],
              [62, 180, -65, 175],
              [62, 180, -65, -85],
              [-177, 65, 65, 85],
              [-177, 65, 65, -175],
              [-177, 65, 180, 85],
              [-177, 65, 180, 180],
              [-177, 180, 65, 85],
              [-177, 180, 65, -175],
              [-177, 180, 65, -105],
              [-177, 180, 180, 85],
              [-177, 180, 180, 180],
              [-177, 180, 180, -85],
              [-177, 180, -65, 105],
              [-177, 180, -65, 175],
              [-177, 180, -65, -85],
              [-67, 180, 65, 85],
              [-67, 180, 65, -175],
              [-67, 180, 65, -105],
              [-67, 180, 180, 85],
              [-67, 180, 180, 180],
              [-67, 180, 180, -85],
              [-67, 180, -65, 105],
              [-67, 180, -65, 175],
              [-67, -167, -65, -85],
              [-62, -68, 180, 85],
              [-62, -68, 180, 180],
              [-62, -68, 180, -85],
              [-62, -68, -65, 175],
              [-62, -68, -65, -85]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Aspargine:
    Charge      = 0
    Index       = 16
    Filename    = fragbuild_dir + "residues/asn.xyz"
    ResName     = "N"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,11]]
    Rotamer     = [[62, -10],
              [62,  30],
              [-174, -20],
              [-177, 30],
              [-65, -20],
              [-65, -75],
              [-65, 120]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Aspartate:
    Charge      = -1
    Index       = 14
    Filename    = fragbuild_dir + "residues/asp.xyz"
    ResName     = "D"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,11]]
    Rotamer     = [[62, -10],
              [62, 30],
              [-177, 0],
              [-177, 65]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)



class Cysteine:
    Charge      = 0
    Index       = 8
    Filename    = fragbuild_dir + "residues/cys.xyz"
    ResName     = "C"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9]]
    Rotamer     = [[62],[-177],[-65]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glutamate:
    Charge      = -1
    Index       = 15
    Filename    = fragbuild_dir + "residues/glu.xyz"
    ResName     = "E"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13],[5,9,13,15]]
    Rotamer     = [[62,  180, -20],
              [70,   -80, 0],
              [-177,  65, 10],
              [-177, 180, 0],
              [-80,  -50, -25],
              [-65,   85, 0],
              [-67,  180, -10],
              [-65,  -65, -40]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glutamine:
    Charge      = 0
    Index       = 17
    Filename    = fragbuild_dir + "residues/gln.xyz"
    ResName     = "Q"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13],[5,9,13,14]]
    Rotamer     = [[62, 180, 20],
              [70,   -75, 0],
              [-177,  65, -100],
              [-177,  65, 60],
              [-177, 180, 0],
              [-65,   85, 0],
              [-67,  180, -25],
              [-65,  -65, -40],
              [-65,  -65, 100]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Glycine:
    Charge      = 0
    Index       = 5
    Filename    = fragbuild_dir + "residues/gly.xyz"
    ResName     = "G"
    BB      = [1, 3, 4]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Histidine:
    Charge      = 1
    Index       = 18
    Filename    = fragbuild_dir + "residues/hip.xyz"
    ResName     = "H"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,12]]
    Rotamer     = [[62, -75],
              [62, 80],
              [-177, -165],
              [-177, -80],
              [-177, 60],
              [-65, -70],
              [-65, 165]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Isoleucine:
    Charge      = 0
    Index       = 3
    Filename    = fragbuild_dir + "residues/ile.xyz"
    ResName     = "I"
    BB      = [1, 3, 6]
    SC      = [[1, 3, 5, 9],[3, 5, 9, 16]]
    Rotamer     = [[ 62, 100],
              [  62, 170],
              [-177,  66],
              [-177, 165],
              [ -65, 100],
              [ -65, 170],
              [ -57, -60]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)



class Leucine:
    Charge      = 0
    Index       = 4
    Filename    = fragbuild_dir + "residues/leu.xyz"
    ResName     = "L"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13]]
    Rotamer     = [[ 62,  80],
              [-177,  65],
              [-172, 145],
              [ -85,  65],
              [ -65, 175]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Lysine:
    Charge      = 1
    Index       = 19
    Filename    = fragbuild_dir + "residues/lys.xyz"
    ResName     = "K"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13],[5,9,13,16],[9,13,16,19]]
    Rotamer     = [[62, 180, 68, 180],
              [62, 180, 180, 65],
              [62, 180, 180, 180],
              [62, 180, 180, -65],
              [62, 180, -68, 180],
              [-177, 68, 180, 65],
              [-177, 68, 180, 180],
              [68, -65, 9, 10],
              [-177, 180, 68, 65],
              [-177, 180, 68, 180],
              [-177, 180, 180, 65],
              [-177, 180, 180, 180],
              [-177, 180, 180, -65],
              [-177, 180, -68, 180],
              [-177, 180, -68, -65],
              [-90, 68, 180, 180],
              [-67, 180, 68, 65],
              [-67, 180, 68, 180],
              [-67, 180, 180, 65],
              [-67, 180, 180, 180],
              [-67, 180, 180, -65],
              [-67, 180, -68, 180],
              [-67, 180, -68, -65],
              [-62, -68, 180, 65],
              [-62, -68, 180, 180],
              [-62, -68, 180, -65],
              [-62, -68, -68, 180]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Methionine:
    Charge      = 0
    Index       = 9
    Filename    = fragbuild_dir + "residues/met.xyz"
    ResName     = "M"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,13],[5,9,13,14]]
    Rotamer     = [[62, 180, 75],
               [62, 180, -75],
               [-177, 65, 75],
               [-177, 65, 180],
               [-177, 180, 75],
               [-177, 180, 180],
               [-177, 180, -75],
               [-67, 180, 75],
               [-67, 180, 180],
               [-67, 180, -75],
               [-65, -65, 103],
               [-65, -65, 180],
               [-65, -65, -70]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)



class Phenylalanine:
    Charge      = 0
    Index       = 11
    Filename    = fragbuild_dir + "residues/phe.xyz"
    ResName     = "F"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,11]]
    Rotamer     = [[62, 90],
               [-177, 80],
               [-65, -85],
               [-65, -30]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Proline:
    Charge      = 0
    Index       = 10
    Filename    = fragbuild_dir + "residues/pro.xyz"
    ResName     = "P"
    BB      = [1, 3, 9]
    SC      = []
    Rotamer     = []
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Serine:
    Charge      = 0
    Index       = 6
    Filename    = fragbuild_dir + "residues/ser.xyz"
    ResName     = "S"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9]]
    Rotamer     = [[62],
               [-177],
               [-65]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Threonine:
    Charge      = 0
    Index       = 7
    Filename    = fragbuild_dir + "residues/thr.xyz"
    ResName     = "T"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9]]
    Rotamer     = [[62],
               [-175],
               [-65]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Tryptophan:
    Charge      = 0
    Index       =  13
    Filename    = fragbuild_dir + "residues/trp.xyz"
    ResName     = "W"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,11]]
    Rotamer     = [[  62,  -90],
               [  62,   90],
               [-177, -105],
               [-177,   90],
               [ -65,  -90],
               [ -65,   -5],
               [ -65,   95]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Tyrosine:
    Charge      = 0
    Index       = 12
    Filename    = fragbuild_dir + "residues/tyr.xyz"
    ResName     = "Y"
    BB      = [1, 3, 6]
    SC      = [[1,3,5,9],[3,5,9,11]]
    Rotamer     = [[62,     90],
               [-177, 80],
               [-65, -85],
               [-65, -30]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)

class Valine:
    Charge      = 0
    Index       = 2
    Filename    = fragbuild_dir + "residues/val.xyz"
    ResName     = "V"
    BB      = [1, 3, 6]
    SC      = [[1, 3, 5, 8]]
    Rotamer     = [[63],
                  [175],
              [-60]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


class Valine:
    Charge      = 0
    Index       = 2
    Filename    = fragbuild_dir + "residues/val.xyz"
    ResName     = "V"
    BB      = [1, 3, 6]
    SC      = [[1, 3, 5, 8]]
    Rotamer     = [[63],
                  [175],
              [-60]]
    def __init__(self):
        self.Mol = pybel.readfile("xyz", self.Filename).next()
        self.AwesomeMol = ReadCoordzAwesome(self.Filename)


