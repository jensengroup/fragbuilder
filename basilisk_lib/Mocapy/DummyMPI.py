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
Dummy pyMPI stubs to run without pyMPI on a single CPU.
"""


rank=0
size=1

SUM=0


def bcast(t):
    "Returns t"
    return t

def barrier():
    "Returns immediately"
    return

def reduce(x, operation):
    "Returns x"
    return x

def gather(l):
    "Returns l"
    return l

def scatter(l):
    "Returns l"
    return l


SUM = None
