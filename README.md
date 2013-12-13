fragbuilder
===========

fragbuilder is a tool to create, setup and analyze QM calculations on peptides.

#### Requires:
Python 2.x
Openbabel 2.x.x (with Python bindings installed)
Numpy (pretty much any version should be fine)


Note: fragbuilder will run with most versions of Open Babel, however there was bug in Open Babel which prevented some dihedral angles to be set accurately. If you experience this while using fragbuilder you will have to update your Open Babel.

I have compiled a short guide to install the latest Open Babel here:
http://combichem.blogspot.dk/2013/12/compiling-open-babel-with-python.html

How to use
==========

0) Checkout a version of fragbuilder from this repository

1) Setup relevant file paths in settings.py

2) If you haven't got BASILISK or BackboneDBN installed, then disable import of these modules in fragbuilder.py

3) Export the frabuilder directory to your PYTHONPATH eg:

export PYTHONPATH=/home/andersx/programs/fragbuilder:$PYTHONPATH

4) In your Python:

import fragbuilder


5) Look in examples/\*.py for inspiration

6) Type help(fragbuilder) in a Python shell to access the documentation.

help(fragbuilder.Peptide) will give you the user documentation for most of fragbuilder's functionality.

How to cite:
===============================

For now do something like this:
"fragbuilder" (2013) Anders S. Christensen, Jan H. Jensen, https://github.com/andersx/fragbuilder 
There will be a proper paper on arXiv shortly.

The use of the functions from the BASILISK library (e.g. Peptide.sample\_bb\_angles(), Peptide.sample\_chi\_angles(), etc.) should cite:

Tim Harder, Wouter Boomsma,  Martin Paluszewski, Jes Frellesen, Kristoffer E. Johansson, and Thomas Hamelryck,  (2010). Beyond rotamers: a generative, probabilistic model of
side chains in proteins. BMC Bioinformatics, 11:306–318.



Additional use of Open Babel should cite:

Noel M. O'Boyle, Michael Banck, Craig A James, Chris Morley, Tim Vandermeersch and Geoffrey R Hutchison (2011) Open Babel: An open chemical toolbox. Journal of Cheminformatics, 3:33-46.



Licensing
=========
fragbuilder is licensed under the two-clause BSD license which mean you can pretty much do whatever you want as long as you cite properly and/or mentions the name of the authors. 

The BASILISK library which is included is licensed under the terms of the Gnu General Public License. Is also available from SourceForge: https://sourceforge.net/projects/basilisk-dbn/



