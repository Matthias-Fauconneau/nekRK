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
 

class Enhancement(AbstractNode):


    tag = "enhancement"


    def notify(self, parent):
        parent.onEnhancement(self._species, self._coefficient)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._species = attributes["species"]
        self._coefficient = float(attributes["factor"])

        return
            

# version
__id__ = "$Id: Enhancement.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
