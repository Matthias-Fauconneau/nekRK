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


class Declaration:


    def __init__(self, locator):
        self.locator = locator
        return


    def __str__(self):
        if self.locator:
            return "source='%s', line=%d, column=%d" % (
                self.locator.filename, self.locator.line, self.locator.column)

        return "source=<unknown>"

# version
__id__ = "$Id: Declaration.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

# End of file
