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


class SpeciesSection(Token):


    pattern = r"[Ss][Pp][Ee][Cc]([Ii][Ee][Ss])?"


    def identify(self, auth):
        return auth.aSpeciesSection(self)


    def __str__(self): 
        return "{Species section}"


# version
__id__ = "$Id: SpeciesSection.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
