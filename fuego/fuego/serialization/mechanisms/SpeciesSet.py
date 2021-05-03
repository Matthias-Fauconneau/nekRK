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


class SpeciesSet(EntitySet):


    def species(self, symbol, locator=None):
        from Species import Species
        species = Species(self.size(), symbol, locator)
        self.insert(symbol, species)
        return species



# version
__id__ = "$Id: SpeciesSet.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

# End of file
