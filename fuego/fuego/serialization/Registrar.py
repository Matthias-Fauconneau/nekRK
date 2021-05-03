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

from pyre.util.ResourceManager import ResourceManager


class Registrar(ResourceManager):


    def __init__(self):
        ResourceManager.__init__(self, "fuego")
        return


# version
__id__ = "$Id: Registrar.py,v 1.1.1.1 2007-09-13 18:17:28 aivazis Exp $"

#  End of file 
