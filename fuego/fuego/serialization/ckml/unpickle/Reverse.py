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
 

class Reverse(AbstractNode):


    tag = "reverse"


    def notify(self, parent):
        parent.onReverse(self._parameters)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        A = attributes["A"]
        beta = attributes["beta"]
        E = attributes["E"]

        self._parameters = (A, beta, E)
        return
            

# version
__id__ = "$Id: Reverse.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
