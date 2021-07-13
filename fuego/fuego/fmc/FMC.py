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

from pyre.applications.Script import Script


class FMC(Script):


    class Inventory(Script.Inventory):

        import pyre.inventory

        out = pyre.inventory.outputFile(name="name")
        out.meta['tip'] = "the name of the output file"
        
        thermo = pyre.inventory.str(name="thermo")
        thermo.meta['tip'] = "the name of the thermodynamic properties database"
        
        mechanism = pyre.inventory.str(name="mechanism", default="GRIMech-3.0.ck2")
        mechanism.meta['tip'] = "the mechanism file to be processed"

        read = pyre.inventory.str(name="read")
        read.meta['tip'] = "the format of the input file; default is chemkin"
        
        write = pyre.inventory.str(name="write", default="c")
        write.meta['tip'] = "the format of the output file; default is c"

        import pyre.weaver
        pickler = pyre.inventory.facility("pickler", factory=pyre.weaver.weaver)
        pickler.meta['tip'] = 'the output generator'


    def main(self, *args, **kwds):

        import pyre.monitors
        import fuego.serialization

        timer = pyre.monitors.timer("fuego")
        if not self.read:
            print "loading '%s'" % (self.mechanism)
        else:
            print "loading '%s' using '%s' parser" % (self.mechanism, self.read)

        timer.start()

        mechanism = fuego.serialization.mechanism()
        if self.thermo:
            mechanism.externalThermoDatabase(self.thermo)
        mechanism = fuego.serialization.load(
            filename=self.mechanism, format=self.read, mechanism=mechanism)
    
        print "    ... done (%g sec)" % timer.stop()

        timer.reset()
        timer.start()
        print "saving in '%s' as '%s' format" % (self.out.name, self.write)
        lines = fuego.serialization.save(
            weaver=self.pickler, mechanism=mechanism, stream=self.out, format=self.write)
        print "    ... done (%g sec)" % timer.stop()

        return


    def __init__(self):
        Script.__init__(self, "fmc")
        self.out = None
        self.thermo = ""
        self.mechanism = ""
        self.read = ""
        self.write = ""
        self.pickler = None
        
        return


    def _configure(self):
        Script._configure(self)

        self.out = self.inventory.out
        self.thermo = self.inventory.thermo
        self.mechanism = self.inventory.mechanism
        self.read = self.inventory.read
        self.write = self.inventory.write
        self.pickler = self.inventory.pickler
        
        return


# version
__id__ = "$Id: FMC.py,v 1.1.1.1 2007-09-13 18:17:28 aivazis Exp $"

#
# End of file
