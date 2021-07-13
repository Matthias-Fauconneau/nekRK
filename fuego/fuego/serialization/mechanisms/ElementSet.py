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

from EntitySet import EntitySet


class ElementSet(EntitySet):


    def element(self, symbol, atomicWeight=None, locator=None):
        from Element import Element
        element = Element(self.size(), symbol, atomicWeight, locator)
        self.insert(symbol, element)
        return element



# version
__id__ = "$Id: ElementSet.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

# End of file
