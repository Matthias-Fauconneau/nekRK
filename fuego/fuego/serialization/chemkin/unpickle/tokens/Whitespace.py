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
from RegularExpressions import whitespace


class Whitespace(Token):


    pattern = whitespace


    def identify(self, auth):
        return auth.aWhitespace(self)


    def __str__(self):
         return "{whitespace}"


# version
__id__ = "$Id: Whitespace.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
