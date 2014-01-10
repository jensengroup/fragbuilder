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
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


alphabet            = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
lower_case_alphabet = alphabet.lower()
numbers             = "0123456789"

regular_chars = alphabet + lower_case_alphabet + numbers + "._-"

bb_smiles = dict([('HN', "[$([H]N(C=O)C)]"     ),
                  ('NH', "[$([N](C=O)C)]"      ),
                  ('XX', "[$([N])]"      ),
                  ('CO', "[$([C](NC)=O)]"      ),
                  ('OC', "[$([O]=C(NC)C)]"     ),
                  ('CA', "[$([C](NC)C(=O))]"   ),
                  ('HA', "[$([H]C(NC)C(=O))]"  ),
                  ('CB', "[$([C]C(NC)C(=O))]"  )])


aa1="ACDEFGHIKLMNPQRSTVWY"
aa3=["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
     "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

d1_to_index={}
dindex_to_1={}
d3_to_index={}
dindex_to_3={}
d1_to_d3={}

for i in range(0, 20):
    n1=aa1[i]
    n3=aa3[i]
    d1_to_index[n1]=i
    dindex_to_1[i]=n1
    d3_to_index[n3]=i
    dindex_to_3[i]=n3

    d1_to_d3[n1]=n3

def index_to_one(index):
    "Amino acid index to single letter (eg 0 to A)"
    return dindex_to_1[index]

def one_to_index(s):
    "Amino acid single letter to index (eg A to 0)"
    return d1_to_index[s]

def index_to_three(i):
    "Amino acid index to three letter (eg 0 to ALA)"
    return dindex_to_3[i]

def three_to_index(s):
    "Amino acid three letter to index (eg ALA to 0)"
    return d3_to_index[s]

def three_to_one(s):
    "Amino acid three letter to single letter (eg ALA to A)"
    i=d3_to_index[s]
    return dindex_to_1[i]

def one_to_three(s):
    "Amino acid single letter to three letter (eg A to ALA)"
    i=d1_to_d3[s]
    return i
