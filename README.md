fragbuilder
===========

fragbuilder is a tool to create, setup and analyze QM calculations on peptides.

#### Requires:
Openbabel 2.x.x (with Python bindings installed)

Guide to install Open Babel:
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


5) Look in examples/*.py for inspiration

6) There is not much documentation at this point


How to cite use of fragbuilder:
===============================
Do something like this:
"fragbuilder" (2013) Anders S. Christensen, Jan H. Jensen, https://github.com/andersx/fragbuilder 
There will be a proper paper on arXiv shortly.

Licensing
=========
fragbuilder is licensed under the two-clause BSD license which mean you can pretty much do whatever you want as long as you cite properly and/or mentions the name of the authors.
