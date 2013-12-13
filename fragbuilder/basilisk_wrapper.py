import basilisk_lib

from math_utils import DEG_TO_RAD, RAD_TO_DEG

import math 

class Basilisk_DBN:
    """
    A wrapper class for access to the BASILISK DBN class
    via fragbuilder.
    """
    def __init__(self):
        """
        Inititalize a DBN object from basilisk_lib.

        """

        self.dbn = basilisk_lib.basilisk_dbn()

    def _to_deg(self, a):
        """
        Return an angle in degrees in the range [-180, 180]. Input is 
        an angle in radians in the range [0, 2*pi].

        Arguments:
        a -- Input angle in radians.

        """

        if a > math.pi:
            a = a - 2.0 * math.pi
        return a / math.pi * 180.0

    def get_sample(self, aa, bb_angles=None):
        """
        Draws a sample from BASILISK. Returns a tuple of:
        ([chi_angles], [phi, psi], log_likelihood)

        The side chain angles can be sampled using dependence on the
        backbone angles (otherwise this dependence is integrated out).
        This is controlled by supplying a set of backbone angles as the
        bb_angles keywords arguments.

        NOTE: Angles are returned in units of degrees in the range
        of [-180.0, 180.0].

        Arguments:
        aa -- Amino acid type in one letter notation (e.g. Alanine is "A").

        Keyword arguments:
        bb_angles -- Optional [phi, psi] backbone angles for the residues
        (default None).

        """
        chis = []
        bbs = []
        ll = 0.0

        aa = basilisk_lib.basilisk_utils.d1_to_index[aa]

        #chis, bbs, ll = self.dbn.get_sample(aa)
        if bb_angles is None:
            chis, bbs, ll = self.dbn.get_sample(aa)
        else:
            phi = bb_angles[0] * DEG_TO_RAD
            psi = bb_angles[1] * DEG_TO_RAD

            chis, bbs, ll = self.dbn.get_sample(aa, phi, psi)

        chis = [self._to_deg(chi) for chi in chis]
        bbs = [self._to_deg(bb) for bb in bbs]

        return chis, bbs, ll
