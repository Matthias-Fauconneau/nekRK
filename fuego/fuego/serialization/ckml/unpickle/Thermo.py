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

from AbstractNode import AbstractNode
 

class Thermo(AbstractNode):


    tag = "thermo"


    def notify(self, parent):
        parent.onThermo(self._thermo)
        return


    def onNASA(self, lowT, highT, parameters):
        self._thermo.append(("NASA", lowT, highT, parameters))
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._thermo = []
        return
            

# version
__id__ = "$Id: Thermo.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
