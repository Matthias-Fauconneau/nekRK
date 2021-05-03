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
 

class ReverseLandauTeller(AbstractNode):


    tag = "rlt"


    def notify(self, parent):
        parent.onReverseLandauTeller(self._parameters)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        B = attributes["B"]
        C = attributes["C"]

        self._parameters = (B, C)
        return
            

# version
__id__ = "$Id: ReverseLandauTeller.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
