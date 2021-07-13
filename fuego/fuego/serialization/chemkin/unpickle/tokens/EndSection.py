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


class EndSection(Token):

    pattern = r"[Ee][Nn][Dd]"


    def identify(self, auth):
        return auth.anEndSection(self)


    def __str__(self):
         return "{section end}"


# version
__id__ = "$Id: EndSection.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
