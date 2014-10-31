# Copyright (c) 2013, Anders S. Christensen
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE# Unit vector normal to two vectors

import numpy as np
import math
from numpy.linalg import norm

RAD_TO_DEG = 180.0 / np.pi
DEG_TO_RAD = np.pi / 180.0

def normalvector(v1, v2):
    J = np.cross(v1, v2)
    return J/(np.sqrt(np.dot(J, J)))

# Vector perpendicular to (v1-v2) and (v3-v2) centered in v2.
def perp_vector(v1, v2, v3):

    J = np.cross(v2 - v3, v2 - v1)
    J /= norm(J)
    J += v2

    return J


# Bond angle between three coordinates
def bondangle(a,b,c):
    # In case numpy.dot() returns larger than 1
    # and we cannot take acos() to that number
    acos_out_of_bound = 1.0
    v1 = a - b
    v2 = c - b
    v1 = v1 / np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = v2 / np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    dot_product = np.dot(v1,v2)

    if dot_product > acos_out_of_bound:
        dot_product = acos_out_of_bound
    if dot_product < -1.0 * acos_out_of_bound:
        dot_product = -1.0 * acos_out_of_bound

    return np.arccos(dot_product)


# Distance between two vectors
def distance(a,b):
    v1 = a - b
    return np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


# Length of a vector
def length(v1):
    return np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)



# Dihedral angle between three points in space
def dihedral(a,b,c,d):
    b1 = b - a
    b2 = c - b
    b3 = d - c
    rb2 = np.sqrt(np.dot(b2,b2))
    b2xb3 = np.cross(b2,b3)
    b1xb2 = np.cross(b1,b2)
    dihedral = math.atan2( np.dot(rb2*b1,b2xb3), np.dot(b1xb2,b2xb3) )
    return dihedral


# Rotate a set of vectors
def rotate(V, J, T):
    x = V[0]
    y = V[1]
    z = V[2]
    u = J[0]
    v = J[1]
    w = J[2]
    a = (u*(u*x + v*y + w*z) + (x * (v*v + w*w) - u *(v*y + w*z))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(-w*y + v*z)*math.sin(T))/(u*u + v*v + w*w)
    b = (v*(u*x + v*y + w*z) + (y * (u*u + w*w) - v *(u*x + w*z))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(w*x - u*z)*math.sin(T))/(u*u + v*v + w*w)
    c = (w*(u*x + v*y + w*z) + (z * (u*u + v*v) - w *(u*x + v*y))*math.cos(T) + math.sqrt(u*u + v*v + w*w)*(-v*x + u*y)*math.sin(T))/(u*u + v*v + w*w)
    return np.array([a, b, c])

def to_rad(degrees):
    degrees = float(degrees)
    return degrees/180.0*math.pi






















