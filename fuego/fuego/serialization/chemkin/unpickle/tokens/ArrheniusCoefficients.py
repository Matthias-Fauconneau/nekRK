#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

#
# It appears that this token class does not get used
# so I placed a Firewall in the constructor
#

from Token import Token
from RegularExpressions import namedNumbers_3


class ArrheniusCoefficients(Token):


    pattern = namedNumbers_3 % ("arrhenius_A", "arrhenius_beta", "arrhenius_E")


    def identify(self, auth):
        return auth.someArrheniusCoefficients(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        import journal
        journal.firewall("fuego").log("arrhenius coefficients token")

        numbers = map(groups.get, arrhenius)

        try:
            self.parameters = map(float, numbers)
        except ValueError:
            # this can't happen because the regexp requires three floats
            import journal
            msg = "Could not convert /%s/ into a list of numbers" % text
            journal.firewall("fuego").log(msg)
            return
            
        return

    
    def __str__(self):
        return "{Arrhenius coefficients: %s}" % self.parameters



# version
__id__ = "$Id: ArrheniusCoefficients.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
