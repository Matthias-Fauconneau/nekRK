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

from Declaration import Declaration


class SpeciesDeclaration(Declaration):


    def __init__(self, symbol, locator=None):
        Declaration.__init__(self, locator)
        self.symbol = symbol
        return


    def __str__(self):
        str = "symbol='%s', " % self.symbol + Declaration.__str__(self)
        return str


# version
__id__ = "$Id: SpeciesDeclaration.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

# End of file
