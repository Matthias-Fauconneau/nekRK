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
 

class Efficiencies(AbstractNode):


    tag = "efficiencies"


    def notify(self, parent):
        parent.onEfficiencies(self._efficiencies)
        return


    def onEnhancement(self, species, coefficient):
        self._efficiencies.append((species, coefficient))
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._efficiencies = []
        return
            

# version
__id__ = "$Id: Efficiencies.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
