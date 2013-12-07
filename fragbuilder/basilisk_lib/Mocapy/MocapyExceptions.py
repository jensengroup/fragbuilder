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
Mocapy exception classes.
"""


class MocapyException(Exception): pass

class MocapyVectorException(MocapyException): pass

class MocapyDBNException(MocapyException): pass

class MocapyVMFException(MocapyException): pass

class MocapyGaussianException(MocapyException): pass

class MocapyDiscreteException(MocapyException): pass

class MocapyDirichletException(MocapyException): pass

class MocapyKentException(MocapyException): pass

class MocapyEMException(MocapyException): pass

class MocapyVMException(MocapyException): pass

