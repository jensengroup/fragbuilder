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

If you're interested in FragBuilder, you might also be interested in PeptideBuilder (very similar to FragBuilder):
- https://github.com/mtien/PeptideBuilder/
- https://peerj.com/articles/80/

How to use
==========

##### 1) Clone fragbuilder from this repository

    git clone https://github.com/andersx/fragbuilder

##### 2) Export the frabuilder directory to your PYTHONPATH eg:

    export PYTHONPATH=/home/andersx/dev/fragbuilder:$PYTHONPATH

##### 3) In your Python:

    import fragbuilder

Documentation
=============

First off, the paper on PeerJ is probably the first place to go in order to get familiar with Fragbuilder. Second, I recommend looking in the examples folder. This should cover most use-cases and functionality. I have also writte a couple of blog-posts (more to come!).

 - See the paper here: http://dx.doi.org/10.7287/peerj.preprints.169v2
 - Look in examples/\*.py for inspiration.
 - Installation help: http://combichem.blogspot.com/2013/12/fragbuilder-installation.html
 - Setup a scan of peptide conformations: http://combichem.blogspot.com/2014/01/fragbuilder-setting-up-scan-of-peptide.html

Below is list of all classes currently in the **fragbuilder** Python library. The manual for each class is accessed though the built-in help function. A short description of each class and the code to access the corresponding help page is as follows:

**fragbuilder.Peptide** class:

This class makes peptides, defines the geometry, optimizes geometries, writes files, etc.

[fragbuilder.Peptide documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.Peptide.html)

**fragbuilder.G09_opt** class:  
**fragbuilder.G09_NMR** class:  
**fragbuilder.G09_energy** class:

These classes will write an Gaussian 09 input file for either geometry optimization, NMR shielding calculations or single point energies.

[fragbuilder.G09_opt documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.G09_opt.html)  
[fragbuilder.G09_NMR documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.G09_NMR.html)  
[fragbuilder.G09_energy documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.G09_energy.html)

**fragbuilder.Basilisk_DBN** class:
 
This class wraps the dbn class from the original BASILISK library. The only real difference is that it returns angles in degrees in order to be compatible with the rest of FragBuilder.
For most practical purposes, this wrapper is not needed. Most likely you want to use the Peptide.sample_bb_angles(args) or Peptide.sample_chi_angles(args) functionality in order to invoke BASILISK library. 

[fragbuilder.Basilisk_DBN documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.Basilisk_DBN.html)


**fragbuilder.PDB** class:
 
This class is a PDB-filereader which can read backbone and side chain torsion angles from a PDB file. 

[fragbuilder.PDB documentation](https://rawgithub.com/jensengroup/fragbuilder/master/doc/fragbuilder.PDB.html)


Get help, report bugs, feature requests, etc:
===============================

You can use all the tools available here at GitHub or you can also contact me directly at this address: andersx [å] nano.ku.dk

I have the following functionality in the making (may be in the development branch):
 
- Handling of protonation states.
- Module that uses TorusDBN to sample backbone angles.
- Module to write Gaussian spin-spin coupling input files.



How to cite:
===============================

#### FragBuilder:
The FragBuilder paper is currently pending peer-review at PeerJ. Meanwhile, you can cite the preprint, also available from PeerJ.
- Anders S. Christensen, Thomas Hamelryck, Jan H. Jensen (2013) FragBuilder: An efficient Python library to setup quantum chemistry calculations on peptides models. PeerJ PrePrints 1:e169v2 http://dx.doi.org/10.7287/peerj.preprints.169v2

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



