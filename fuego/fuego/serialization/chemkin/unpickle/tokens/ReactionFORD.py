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


class ReactionFORD(Token):


    pattern = r"[Ff][Oo][Rr][Dd]"


    def identify(self, auth):
        return auth.aReactionFORD(self)


    def __str__(self): 
        return "{Reaction order}"


# version
__id__ = "$Id: ReactionFORD.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
