#!/usr/bin/env python
#
#--------------------------------------------------------------------------------
#
#                              Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007  All Rights Reserved
#
#--------------------------------------------------------------------------------

from pyre.xml.Node import Node


class AbstractNode(Node):


    def __init__(self, root, attributes):
        Node.__init__(self, root)
        return


# version
__id__ = "$Id: AbstractNode.py,v 1.1.1.1 2007-09-13 18:17:30 aivazis Exp $"

#
# End of file
