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
Miscellaneous array operations.
"""

from numpy import sqrt, sum, dot, transpose, eye
from random import choice, random
from types import FloatType, IntType



def random_sign():
    """
    Return random sign (ie. 1 or -1)

    @return: random sign
    @rtype: int
    """
    if random()<0.5:
        return -1
    else:
        return 1

def close(x, y, delta=1e-10):
    """
    Return 1 if x is close to y (within delta).
    If x and y are list-like, the test is done for each 
    element pair of x and y.

    @param x, y: test if x is close to y
    @type x,y: int, float, numpy array, list,...

    @param delta: closeness cutoff
    @type delta: float

    @return: 1 if x is close to y, 0 otherwise
    @rtype: int (0 or 1)
    """
    type_x=type(x)
    if type_x==FloatType or type_x==IntType:
        # Float or Int
        if abs(x-y)<delta:
            return 1
        else:
            return 0
    else:
        # Array, list,...
        for i in range(0, len(x)):
            if abs(x[i]-y[i])>delta:
                return 0
        return 1

def gen_refmat(p, q):
    """
    Return a (left multiplying) matrix that reflects p onto q.
    p and q can be of arbitrary dimension.

    @param p,q: vectors, same length
    @type p,q: numpy array (shape=n, with n>1)

    @return: matrix that reflects p onto q
    @rtype: numpy array (shape=(n,n))
    """
    q=q/norm(q)
    p=p/norm(p)
    pq=p-q
    npq=norm(pq)
    if npq<1e-5:
        n=p.shape[0]
        return eye(n)
    pq=pq/npq
    n=p.shape[0]
    pq.shape=(n,1)
    i=eye(n)
    ref=i-2*dot(pq, transpose(pq))
    return ref


def gen_rotmat(p,q):
    """
    Return a (left multiplying) matrix that rotates p onto q.
    p and q can be of arbitrary dimension.

    @param p,q: vectors, same length
    @type p,q: numpy array (shape=n, with n>1)

    @return: matrix that reflects p onto q
    @rtype: numpy array (shape=(n,n))
    """
    q=q/norm(q)
    p=p/norm(p)
    npq=norm(p-q)
    if npq<1e-5:
        n=p.shape[0]
        return eye(n)
    rot=dot(gen_refmat(q, -p), gen_refmat(p, -p))
    return rot


def argmax(a):
    """
    Return the index of a maximum in sequence a.
    This is an alternative implementation of argmax: 
    if there are equivalent maxima, it picks one at 
    random. Numeric's argmax always returns the index
    of the first maximum.

    @param a: array of numbers
    @type a: list

    @return: index of maximum of a
    @rtype: int
    """
    m=a[0]
    indices=[0]
    for i in range(1, len(a)):
        x=a[i]
        if x>m:
            # New maximum
            m=x
            indices=[i]
        elif x==m:
            # Same maximum: append index
            indices.append(i)
    # Pick a random index of a maximum
    return choice(indices)


def norm(a):
    """
    Return norm of array a. 
    Basically a shortcut to avoid creating a Vector object.

    @param a: vector
    @type a: numpy array (shape=n, with n>1)

    @return: vector norm
    @rtype: double
    """
    return sqrt(sum(a*a))

