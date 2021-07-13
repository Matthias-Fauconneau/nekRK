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


class ReactionSection(Token):


    pattern = r"[Rr][Ee][Aa][Cc]([Tt][Ii][Oo][Nn][Ss])?"


    def identify(self, auth):
        return auth.aReactionSection(self)


    def __str__(self):
         return "{Reaction section}"


# version
__id__ = "$Id: ReactionSection.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
