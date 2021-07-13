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

from Token import Token


class ReactionUnitsA(Token):

    
    pattern = r"([Mm][Oo][Ll][Ee][Ss])|([Mm][Oo][Ll][Ee][Cc][Uu][Ll][Ee][Ss])"


    def identify(self, auth): 
        return auth.aReactionUnitsA(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.units_A = self.lexeme.lower()
        return


    def __str__(self):
        return "{Reaction units for A: %s}" % self.units_A


# version
__id__ = "$Id: ReactionUnitsA.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
