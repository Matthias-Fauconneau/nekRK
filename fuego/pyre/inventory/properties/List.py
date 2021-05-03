#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2005  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from pyre.inventory.Property import Property


class List(Property):


    def __init__(self, name, default=None, meta=None, validator=None):
        if default is None:
            default = []
        Property.__init__(self, name, "list", default, meta, validator)
        return


    def _cast(self, text):
        if isinstance(text, basestring):
            if text and text[0] in '[({':
                text = text[1:]
            if text and text[-1] in '])}':
                text = text[:-1]
                
            value = text.split(",")
        else:
            value = text

        if isinstance(value, list):
            return value
            
        raise TypeError("property '%s': could not convert '%s' to a list" % (self.name, text))
    

# version
__id__ = "$Id: List.py,v 1.2 2007-01-11 21:33:29 aivazis Exp $"

# End of file 
