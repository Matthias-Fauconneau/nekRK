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

from BaseScanner import BaseScanner


class Species(BaseScanner):


    def _tokenClasses(self):
        from fuego.serialization.chemkin.unpickle.tokens import speciesTokenClasses
        return speciesTokenClasses()


# version
__id__ = "$Id: Species.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
