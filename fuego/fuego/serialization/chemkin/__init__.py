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


def format():
    return "chemkin"


def extensions():
    return [".ck2"]


def parser():
    from unpickle.parsers.ChemkinParser import ChemkinParser
    return ChemkinParser()


def pickler():
    from pickle.ChemkinPickler import ChemkinPickler
    return ChemkinPickler()


# version
__id__ = "$Id: __init__.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
