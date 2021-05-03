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
 

class NASA_a5(AbstractNode):


    tag = "a5"


    def notify(self, parent):
        parent.onNASA_a5(self._value)
        return


    def content(self, text):
        self._value = float(text)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._value = 0.0
        return
            

# version
__id__ = "$Id: NASA_a5.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
