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

from Entity import Entity


class Species(Entity):


    def thermalParametrization(self, type, lowT, highT, locator, parameters):
        if type == "NASA":
            from NASA import NASA
            model = NASA(lowT, highT, locator)
            model.parameters = parameters
        else:
            import journal
            journal.firewall("fuego").log("unknown thermal parametrization type '%s'" % type)
            return

        self.thermo.append(model)
        return


    def __init__(self, id, symbol, locator=None):
        Entity.__init__(self, id, locator)
        self.symbol = symbol
        self.phase = None
        self.composition = []
        self.thermo = []
        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol
        str += ", phase='%s'" % self.phase
        str += ", composition=%s" % self.composition


        for p in self.thermo:
            str += ", thermo=([%g, %g]: %s)" % (p.lowT, p.highT, p.parameters)

        str += ", source=" + Entity.__str__(self)
        return str


# version
__id__ = "$Id: Species.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

# End of file
