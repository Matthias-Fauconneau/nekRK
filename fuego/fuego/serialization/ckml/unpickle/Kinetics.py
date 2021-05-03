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
 

class Kinetics(AbstractNode):


    tag = "kinetics"


    def onElement(self, element):
        return


    def onSpecies(self, declaration):
        return


    def onReaction(self, reaction):
        return


    def notify(self, parent):
        parent.onKinetics(self._mechanism)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._mechanism = root.mechanism()
        return
            

# version
__id__ = "$Id: Kinetics.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
