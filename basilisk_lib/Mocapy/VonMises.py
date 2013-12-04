#   An extension to: Mocapy
#
#   Copyright (C) 2007 Jes Frellsen, Ida Moltke and Martin Thiim
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
Von Mises sampler and density.
"""



from math import pi, acos, cos, sin
from numpy import log, exp, i0            # i0 is modified Bessel function of first kind
from random import vonmisesvariate

i0_float = lambda x: float(i0(float(x)))

def to_radian(coords):
   if coords[1]>0:
      return acos(coords[0])
   else:
      return 2*pi-acos(coords[0])

def to_coords(angle):
   return array([cos(angle), sin(angle)])


def mod2pi(x):
   if x >= 2*pi:
      return mod2pi(x - 2*pi) 
   if x < 0:
      return mod2pi(x + 2*pi)
   else:
      return x


class VMSampler:
    """
    Sample from the Von Mises distribution.
    """
    def __init__(self, k, mu):
        """
        @param k: Kappa
        @type k: float

        @param mu: a 2-dimensional unit vector, VM mean direction
        @type mu: numpy array
        """
        # Make sure mu has the right shape
        if mu<0 or mu>2*pi:
            raise VMException, "Mu has to lie in the interval between 0 and 2pi"
        else:
	    # If so the parameters are stored (mu as radian)
           self.mu = mu
	   self.k = k
	   

    def __call__(self):
        """
        Generate a random vector sampled from the current Von Mises
        distribution.

        @return: random angle
        @rtype: float
        """

        try:
           value = vonmisesvariate(float(self.mu),float(self.k))
        except:
           print "ERROR: vonmisesvariate(%s,%s)" % (float(self.mu), float(self.k))
           raise
        
        return mod2pi(value)


class VMDens:
    """
    VM density.
    """
    def __init__(self, k, mu):
        """
        @param k: Kappa
        @type k: float

        @param mu: an angle in the interval [0,2pi[, VM mean direction
        @type mu: float
        """
	
	# Make sure mu is in the correct interval
        if mu<0 or mu>2*pi:
	    raise Exception, "Mu has to lie between 0 and 2pi"

	#Store mu and kappa
        self.mu = mu
        self.k = k

	# Calculate and store the logarithm of the numerator in the
        # density function (= the reciprocal of the normalising constant)
        self.log_denom = log(2*pi*i0_float(k))



    def __call__(self, x, log_space=0):
        """
        Return VM density at x

        @rtype: float
        """
        
        # Calculate log of the numerator in the density function
        lognum = self.k*cos(x-self.mu)       

        # Calculate log of the density 
        logdens=lognum-self.log_denom
        
        if log_space:
            return logdens
        else:
            return exp(logdens)
