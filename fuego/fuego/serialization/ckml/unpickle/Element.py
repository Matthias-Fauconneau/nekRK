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
 

class Element(AbstractNode):


    tag = "element"


    def notify(self, parent):
        parent.onElement(self._element)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        symbol = attributes["id"]
        atomicWeight = attributes.get("atomicWeight")
        locator = self.document.locator

        self._element = self.document.mechanism().newElement(symbol, atomicWeight, locator)
        return
            

# version
__id__ = "$Id: Element.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
