FragBuilder
===========

FragBuilder is a tool to create, setup and analyze QM calculations on peptides.

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
 - Type help(fragbuilder.Peptide) in a Python shell to access most of the documentation. Similar for the other modules (G09_*, etc)
 

How to cite:
===============================

#### FragBuilder:
For now do something like this:
 - FragBuilder - An efficient Python library to setup quantum chemistry calculations on peptides models (2013) Anders S. Christensen, Jan H. Jensen, https://github.com/andersx/fragbuilder 
There will be a proper paper on arXiv shortly.

#### BASILISK:
The use of the functions from the BASILISK library (e.g. Peptide.sample\_bb\_angles(), Peptide.sample\_chi\_angles(), etc.) should cite:

- Tim Harder, Wouter Boomsma,  Martin Paluszewski, Jes Frellesen, Kristoffer E. Johansson, and Thomas Hamelryck,  (2010). Beyond rotamers: a generative, probabilistic model of
side chains in proteins. BMC Bioinformatics, 11:306–318.


#### Open BAbel:
Additional use of Open Babel should cite:

- Noel M. O'Boyle, Michael Banck, Craig A James, Chris Morley, Tim Vandermeersch and Geoffrey R Hutchison (2011) Open Babel: An open chemical toolbox. Journal of Cheminformatics, 3:33-46.



Licensing
=========
FragBuilder is licensed under the two-clause BSD license which mean you can pretty much do whatever you want as long as you cite properly and/or mentions the name of the authors. 

The BASILISK library which is included is licensed under the terms of the Gnu General Public License. Is also available from SourceForge: https://sourceforge.net/projects/basilisk-dbn/



