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
Von Mises estimator.
"""

from math import sqrt, pow

MAX_KAPPA_VALUE = 500

def estimate_kappa(N, r, mu_coords):
    if N==0:
        print "WARNING: There is no data to estimate mu and kappa from"
        return MAX_KAPPA_VALUE+1

    else:    
        sum_cos = float(r[0])
        sum_sin = float(r[1])
        N = float(N)
        c = (1.0 / N) * sum_cos
        s = (1.0 / N) * sum_sin
        p = sqrt(pow(c,2) + pow(s,2))

        if(p < 2.0 / 3):
            return p*((2-pow(p,2)) / (1 - pow(p,2)))
        else:
            return (p + 1) / (4*p*(1-p))

