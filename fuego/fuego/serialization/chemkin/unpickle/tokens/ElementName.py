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

from Token import Token
from RegularExpressions import whitespaceOpt, namedElement, namedInlineNumber


class ElementName(Token):


    _patternNames = ("element_name", "element_weight")
    _patternTemplate = namedElement + whitespaceOpt + r"(" + namedInlineNumber + r")?" \

    pattern = _patternTemplate % _patternNames


    def identify(self, auth):
        return auth.anElementName(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.name = groups["element_name"]
        self.weight = groups["element_weight"]
        return


    def __str__(self):
        str = "{element: name=" + self.name
        if self.weight:
            str = str + ", weight=" + self.weight
        str = str + "}"
        return str


# version
__id__ = "$Id: ElementName.py,v 1.1.1.1 2007-09-13 18:17:31 aivazis Exp $"

#
# End of file
