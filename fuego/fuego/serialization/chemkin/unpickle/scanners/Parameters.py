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


class Parameters(BaseScanner):


    def _tokenClasses(self):
        from fuego.serialization.chemkin.unpickle.tokens import parameterTokenClasses
        return parameterTokenClasses()


# version
__id__ = "$Id: Parameters.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
