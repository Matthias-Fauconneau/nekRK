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


class ElementDb:


    def element(self, element):
        self._elements.append(element)
        self._index[element.symbol] = element
        return


    def size(self):
        return len(self._elements)


    def find(self, symbol=None):
        if symbol:
            return self._index.get(symbol)

        return self._elements


    def __init__(self):
        self._index = {}
        self._elements = []
        return


# version
__id__ = "$Id: ElementDb.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

# End of file
