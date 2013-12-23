FragBuilder
===========

FragBuilder is a tool to create, setup and analyze QM calculations on peptides.

![alt text](https://dl.dropboxusercontent.com/u/17435887/fragbuilder/model_angles_crop.png "Di-alanine peptide")
#### Requires:
 - Python 2.x
 - Openbabel 2.x.x (with Python bindings installed)
 - Numpy (pretty much any version should be fine)


#### Note: 
FragBuilder will run with most versions of Open Babel, however there was bug in Open Babel which prevented some dihedral angles to be set accurately. If you experience this while using fragbuilder you will have to update your Open Babel.

I have compiled a short guide to install the latest Open Babel here:
http://combichem.blogspot.dk/2013/12/compiling-open-babel-with-python.html

How to use
==========

##### 1) Clone fragbuilder from this repository

    git clone https://github.com/andersx/fragbuilder

##### 2) Export the frabuilder directory to your PYTHONPATH eg:

    export PYTHONPATH=/home/andersx/dev/fragbuilder:$PYTHONPATH

##### 3) In your Python:

    import fragbuilder

#### Documentation:

 - Look in examples/\*.py for inspiration
 - Documentation (i.e. the manual) is written in doc-strings. See the modules below.
 - Also see the paper -- a link is coming soon(TM)!

### The turbo-manual:

Below is list of all classes currently in the **fragbuilder** Python library. The manual for each class is accessed though the built-in help function. A short description of each class and the code to access the corresponding help page is as follows:

**fragbuilder.Peptide** class:

This class makes peptides, defines the geometry, optimizes geometries, writes files, etc.

    from fragbuilder import Peptide
    help(Peptide)

**fragbuilder.G09_opt** class:  
**fragbuilder.G09_NMR** class:  
**fragbuilder.G09_energy** class:

These classes will write an Gaussian 09 input file for either geometry optimization, NMR shielding calculations or single point energies.

    from fragbuilder import G09_opt, G09_NMR, G09_energy
    help(G09_opt)
    help(G09_NMR)
    help(G09_energy)

**fragbuilder.Basilisk_DBN** class:
 
This class wraps the dbn class from the original BASILISK library. The only real difference is that it returns angles in degrees in order to be compatible with the rest of FragBuilder

    from fragbuilder import Basilisk_DBN
    help(Basilisk_DBN)


**fragbuilder.PDB** class:
 
This class is a PDB-filereader which can read backbone and side chain torsion angles from a PDB file. 

    from fragbuilder import PDB
    help(PDB)

NOTE: The **fragbuilder.PDB** class is under development and currently requires Biopython to function (The rest of FragBuilder will work even if Biopython is not installed). I am in the process of removing either this class from FragBuilder or remove the dependency on Biopython.

Get help, report bugs, request features, etc:
===============================

You can use all the tools available here at GitHub or you can also contact me directly at this address: andersx [å] nano.ku.dk


How to cite:
===============================

#### FragBuilder:
The FragBuilder paper is currently pending peer-review at PeerJ. Meanwhile, you can cite the preprint, also available from PeerJ.
Anders S. Christensen, Thomas Hamelryck, Jan H. Jensen (2013) FragBuilder: An efficient Python library to setup quantum chemistry calculations on peptides models. PeerJ PrePrints 1:e169v1 http://dx.doi.org/10.7287/peerj.preprints.169v1

#### BASILISK:
The use of the functions from the BASILISK library (e.g. Peptide.sample\_bb\_angles(), Peptide.sample\_chi\_angles(), etc.) should cite:

- Tim Harder, Wouter Boomsma,  Martin Paluszewski, Jes Frellesen, Kristoffer E. Johansson, and Thomas Hamelryck  (2010). Beyond rotamers: a generative, probabilistic model of
side chains in proteins. BMC Bioinformatics, 11:306–318.


#### Open Babel:
Additional use of Open Babel should cite:

- Noel M. O'Boyle, Michael Banck, Craig A James, Chris Morley, Tim Vandermeersch and Geoffrey R Hutchison (2011) Open Babel: An open chemical toolbox. Journal of Cheminformatics, 3:33-46.



Licensing
=========
FragBuilder is licensed under the two-clause BSD license which mean you can pretty much do whatever you want as long as you cite properly and/or mentions the name of the authors. 

The BASILISK library which is included is licensed under the terms of the Gnu General Public License. Is also available from SourceForge: https://sourceforge.net/projects/basilisk-dbn/



