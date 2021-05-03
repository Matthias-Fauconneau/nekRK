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
    return "ckml"


def extensions():
    return [".ckml"]


def parser():
    from unpickle.CkmlParser import CkmlParser
    return CkmlParser()


def pickler():
    from pickle.CkmlPickler import CkmlPickler
    return CkmlPickler()


# version
__id__ = "$Id: __init__.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#  End of file 
