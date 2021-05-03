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


def sections():
    from Sections import Sections
    return Sections()


def elements():
    from Elements import Elements
    return Elements()


def species():
    from Species import Species
    return Species()


def thermo():
    from Thermo import Thermo
    return Thermo()


def reactions():
    from Reactions import Reactions
    return Reactions()


def parameters():
    from Parameters import Parameters
    return Parameters()


# version
__id__ = "$Id: __init__.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

# End of file
