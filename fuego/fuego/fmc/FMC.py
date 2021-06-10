#!/usr/bin/env python
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
from pyre.applications.Application import Application
class FMC(Application):
    def run(self):
        import fuego
        import pyre.monitors
        save          = self.inventory.name
        input         = self.inventory.input
        output        = self.inventory.output
        mechanismFile = self.inventory.mechanism
        thermo        = self.inventory.thermo
        trans         = self.inventory.trans
        #AF
        save_chop   = save.split(".")[0]+"_1.cpp"
        save_header = "chemistry_file.H"
        mech_header = "mechanism.h"
        assert(not input)
        mechanism = fuego.serialization.mechanism()
        if thermo:
            mechanism.externalThermoDatabase(thermo)
        if trans:
            mechanism.externalTransDatabase(trans)
        mechanism = fuego.serialization.load(
            filename=mechanismFile, format=input, mechanism=mechanism)
        lines        = fuego.serialization.save(mechanism, output)
        outputFileHeader  = self._openOutput(save_header)
        count_lines = 0
        for line in lines:
            outputFileHeader.write(line)
            outputFileHeader.write('\n')
            count_lines += 1
        outputFile = self._openOutput(save)
        for line in lines:
            if ('#ifndef MECHANISM_h') in line:
                line_start_mech_header = count_lines
                break;
            outputFile.write(line)
            outputFile.write('\n')
            count_lines += 1
        MechHeaderFile = self._openOutput(mech_header)
        for line in lines:
            MechHeaderFile.write(line)
            MechHeaderFile.write('\n')

    def __init__(self):
        Application.__init__(self, "fmc")

    def _openOutput(self, name):
        if name == "stdout":
            import sys
            return sys.stdout
        return file(name, "w")

    class Inventory(Application.Inventory):
        import pyre.properties
        inventory = [
            pyre.properties.str("name", default="stdout"),
            pyre.properties.str("mechanism", default="GRIMech-3.0.ck2"),
            pyre.properties.str("thermo", default=""),
            pyre.properties.str("trans", default=""),
            pyre.properties.str("input", default=""),
            pyre.properties.str("output", default="c"),
        ]
