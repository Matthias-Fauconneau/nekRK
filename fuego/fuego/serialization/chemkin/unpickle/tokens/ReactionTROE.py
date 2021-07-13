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


class ReactionTROE(Token):


    pattern = r"[Tt][Rr][Oo][Ee]"


    def identify(self, auth):
        return auth.aReactionTROE(self)


    def __str__(self):
         return "{TROE fall off}" 


# version
__id__ = "$Id: ReactionTROE.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
