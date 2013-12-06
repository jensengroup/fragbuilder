#   Mocapy: A parallelized Dynamic Bayesian Network toolkit
#
#   Copyright (C) 2004 Thomas Hamelryck 
#
#   This library is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with BASILISK.  If not, see <http://www.gnu.org/licenses/>.

"""
A parallelizable toolkit for learning and inference in Dynamic Baysian 
Networks (DBNs).

The following classes and functions are 'public', ie. the user can use them to build
a DBN or to do learning and inference. The other functions and classes are (mostly)
for internal use.

DBN definition:

    1. L{DBN}: DBN class
    2. L{load_dbn}: function to load a saved L{DBN} object
    3. L{mocapy_seed}: initialize random number generators for reproducibility

The node classes:

    1. L{DiscreteNode}: Discrete node, for symbols or integers

"""

try:
    import mpi
except ImportError:
    # Replace mpi module with dummy module
    # to run on a single processor without MPI
    import DummyMPI
    mpi=DummyMPI

from DBN import DBN, load_dbn, mocapy_seed
from DiscreteNode import DiscreteNode, make_random_cpd, normalize_cpd

# These need scipy
from VMNode import VMNode
