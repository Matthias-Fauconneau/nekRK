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
    return "native"


def extensions():
    return [".py"]


def parser():
    from NativeParser import NativeParser
    return NativeParser()


def pickler():
    from NativePickler import NativePickler
    return NativePickler()


# version
__id__ = "$Id: __init__.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

#  End of file 
