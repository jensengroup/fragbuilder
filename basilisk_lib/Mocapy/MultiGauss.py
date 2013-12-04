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
Density function for a multimensional Gaussian distribution.
"""

from numpy.linalg import inv, det
from numpy import dot, pi, sqrt, power, transpose
from math import exp, log
from types import FloatType


class MultiGauss:
    """
    Density function for multidimensional Gaussian distribution.
    """
    def __init__(self, mean, sigma):
        """
        @param mean: mean of Gaussian density 
        @type mean: float or array of type 'd'

        @param sigma: sigma of Gaussian density
        @type sigma: float or array of type 'd'
        """
        self.sigma=sigma
        self.mean=mean
        if type(mean)==FloatType and type(sigma)==FloatType:
            # 1D Gaussian
            self.c=-log(sigma*sqrt(2*pi))
            self.c2=2*sigma*sigma
            self.one_dim=1
        else:
            # Multidimensional Gaussian
            assert(det(sigma)>0)
            self.dim=mean.shape[0]
            assert(sigma.shape==(self.dim,self.dim))
            sq_det=sqrt(det(self.sigma))
            self.b=1.0/(power(2*pi, self.dim/2.0)*sq_det)
            self.logb=log(self.b)
            self.inv_sigma=inv(self.sigma)
            self.one_dim=0

    def __call__(self, a, log_space=0):
        """
        Return likelihood of a.

        @param a: argument to Gaussian density function
        @type a: float or array of type 'd'

        @param log_space: if 1 return log(density)
        @type log_space: 0 or 1

        @return: Gaussian density at a
        @rtype: float
        """
        if self.one_dim:
            # 1D Gaussian
            d=(a-self.mean)
            x=self.c-d*d/self.c2
            if log_space:
                return x
            else:
                return exp(x)
        else:
            # Multidimensional Gaussian
            p=a-self.mean
            q=dot(self.inv_sigma, p)
            r=-0.5*dot(transpose(p), q)
            if log_space:
                return self.logb+r
            else:
                return self.b*exp(r)

    def __repr__(self):
        return "<%d-dim Gaussian, mean=%s>" % (self.dim, str(self.mean))

if __name__=="__main__":

    from numpy import *
    from numpy.random import *

    seed(1,2)

    mean=random(4)
    mgauss=MultiGauss(mean, random((4,4)))

    print mgauss(mean)
    print mgauss

