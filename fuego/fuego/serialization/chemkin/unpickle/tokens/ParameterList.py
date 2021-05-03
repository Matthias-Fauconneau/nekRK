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

import re
from Token import Token
from RegularExpressions import namedInlineParameters


class ParameterList(Token):

    pattern = namedInlineParameters % "parameter_list"
    whitespace = re.compile("\s+")


    def identify(self, auth):
        return auth.aParameterList(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        self.text = groups["parameter_list"]
        self.parameters = self.whitespace.split(self.text)
        self.count = len(self.parameters)

        return

    
    def __str__(self):
        return "{Parameter list: %s}" % self.parameters


# version
__id__ = "$Id: ParameterList.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
