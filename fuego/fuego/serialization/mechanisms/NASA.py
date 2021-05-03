#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from Thermodynamics import Thermodynamics


class NASA(Thermodynamics):


    def __init__(self, lowT, highT, locator=None):
        Thermodynamics.__init__(self, lowT, highT, locator)
        self.parameters = []
        return


# version
__id__ = "$Id: NASA.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

#  End of file 
