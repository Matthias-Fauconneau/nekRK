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
 

class ReactionRate(AbstractNode):


    tag = "rate"


    def notify(self, parent):
        parent.onReactionRate(self)
        return


    def onEfficiencies(self, parameters):
        self.efficiencies = parameters
        return


    def onArrhenius(self, parameters):
        self.arrhenius = parameters
        return


    def onReverse(self, parameters):
        self.rev = parameters
        return


    def onLowPressure(self, parameters):
        self.low = parameters
        return


    def onSRI(self, parameters):
        self.sri = parameters
        return


    def onTROE(self, parameters):
        self.troe = parameters
        return


    def onLandauTeller(self, parameters):
        self.lt = parameters
        return


    def onReverseLandauTeller(self, parameters):
        self.rlt = parameters
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self.efficiencies =[]

        self.arrhenius = None
        self.rev = None
        self.lt = None
        self.rlt = None


        self.low = None
        self.sri = None
        self.troe = None

        return
            

# version
__id__ = "$Id: ReactionRate.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
