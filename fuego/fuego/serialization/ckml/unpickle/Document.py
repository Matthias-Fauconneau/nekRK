#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.xml.Document import Document as DocumentNode


class Document(DocumentNode):

    tags = [
        "Kinetics",

        "Element",

        "Species", "Composition", "Atom", "Thermo",
        "NASA", "NASA_a1", "NASA_a2", "NASA_a3", "NASA_a4", "NASA_a5", "NASA_a6", "NASA_a7",

        "Reaction", "Duplicate", "ReactionUnits", "Reagents", "Reactant", "Product",
        "ReactionRate", "Efficiencies", "Enhancement",
        "Arrhenius", "Reverse",
        "LowPressure", "SRI", "TROE", "LandauTeller", "ReverseLandauTeller",
        ]


    def mechanism(self):
        return self._mechanism


    def onKinetics(self, mechanism):
        self._mechanism = mechanism
        return


    def __init__(self, mechanism, source):
        DocumentNode.__init__(self, source)
        self._mechanism = mechanism

        return


# version
__id__ = "$Id: Document.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

# End of file
