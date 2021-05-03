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
 

class ReactionUnits(AbstractNode):


    tag = "units"


    def notify(self, parent):
        parent.onReactionUnits(self._reactionUnits)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._reactionUnits = {
            "activation": attributes.get("activation", "cal/mole"),
            "prefactor": attributes.get("prefactor", "mole/cm**3")
            }

        return
            

# version
__id__ = "$Id: ReactionUnits.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
