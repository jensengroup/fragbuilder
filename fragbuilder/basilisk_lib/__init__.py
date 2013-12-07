#  BASILISK: A generative, probabilistic model of side chains in proteins
#
#  Copyright (C) 2010 	Tim Harder and Jes Frellsen 
#
#  BASILISK is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BASILISK is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BASILISK.  If not, see <http://www.gnu.org/licenses/>.

"""
Prerequisites:
---------------------------
BASILISK depends on Biopython being installed on the system for the
handling of protein structures. 
For more information on biopython see
Biopython: http://biopython.org

BASILISK includes a strip-down version of Mocapy 0.726.
For more information on Mocapy see
Mocapy: http://sourceforge.net/projects/mocapy/

"""

from basilisk_utils import *
from basilisk_dbn import *

__all__ = ["basilisk_utils", "basilisk_dbn"]

