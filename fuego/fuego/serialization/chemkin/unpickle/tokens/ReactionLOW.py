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


class ReactionLOW(Token):


    pattern = r"[Ll][Oo][Ww]"


    def identify(self, auth):
        return auth.aReactionLOW(self)


    def __str__(self):
        return "{LOW fall off}"


# version
__id__ = "$Id: ReactionLOW.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
