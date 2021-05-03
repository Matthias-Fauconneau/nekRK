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


class ReactionRLT(Token):


    pattern = r"[Rr][Ll][Tt]"


    def identify(self, auth): 
        return auth.aReactionRLT(self)


    def __str__(self):
        return "{RLT reaction}"


# version
__id__ = "$Id: ReactionRLT.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
