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
 

class Duplicate(AbstractNode):


    tag = "duplicate"


    def notify(self, parent):
        parent.onDuplicate()
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        return
            

# version
__id__ = "$Id: Duplicate.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
