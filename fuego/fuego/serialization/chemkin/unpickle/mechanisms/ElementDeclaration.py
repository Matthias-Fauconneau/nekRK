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


class ElementDeclaration(Declaration):


    def __init__(self, symbol, weight=None, locator=None):
        Declaration.__init__(self, locator)
        
        self.symbol = symbol
        self.weight = weight

        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol
        if self.weight:
            str += "weight=%g" % self.weight

        str += ", " + Declaration.__str__(self)

        return str


# version
__id__ = "$Id: ElementDeclaration.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

# End of file
