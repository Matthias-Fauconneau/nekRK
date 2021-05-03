#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2007  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# access to the chemistry related files


def mechanismPath():
    import pyre

    pathlist = FUEGO_MECHANISM_PATH

    if not pathlist:
        pathlist = [ defaultMechanismDirectory() ]
        
    return pathlist
    


def defaultMechanismDirectory():
    
    import os

    dir = os.path.abspath(os.path.join(home(), FUEGO_ETC_DIR, FUEGO_MECHANISM_DIR))

    return dir


def chemkinMechanismFile(filename):
    import os
    return os.path.join(chemkinMechanismDirectory(), filename)


def home():
    return __path__[0]
    

def copyright():
    return "fuego: Copyright (c) 1998-2007 Michael A.G. Aivazis"


FUEGO_ETC_DIR = "../etc/fuego"
FUEGO_MECHANISM_DIR = "mechanisms"
FUEGO_MECHANISM_PATH = [".", "./mechanisms"]


# version
__id__ = "$Id: __init__.py,v 1.1.1.1 2007-09-13 18:17:28 aivazis Exp $"

# End of file
