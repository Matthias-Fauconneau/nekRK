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


class NativeParser:


    def parse(self, mechanism, file):
        context = {
            "mechanism": mechanism
            }

        exec file in context
        
        return mechanism



# version
__id__ = "$Id: NativeParser.py,v 1.1.1.1 2007-09-13 18:17:32 aivazis Exp $"

#  End of file 
