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


class ReactionDuplicate(Token):
    
    pattern = r"[Dd][Uu][Pp]([Ll][Ii][Cc][Aa][Tt][Ee])?"


    def identify(self, auth):
        return auth.aReactionDuplicate(self)


    def __str__(self):
        return "{duplicate reaction}"


# version
__id__ = "$Id: ReactionDuplicate.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
