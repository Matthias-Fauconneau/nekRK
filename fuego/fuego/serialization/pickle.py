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

def save(weaver, mechanism, stream, format="chemkin"):

    # mga 20070913: FIXME: impossible to use with an driver:
    # who can configure a weaver manually?

    import journal
    journal.debug("fuego").log("pickling mechanism, format='%s')" % format)

    mill = pickler(format)
    if not pickler:
        journal.error("fuego").log("unknown mechanism file format '%s'" % format)
        return []
        
    weaver.renderer = mill
    return weaver.weave(mechanism, stream)


# factory methods for the serializers

def pickler(format="chemkin"):
    factory = registrar().find(format)
    if factory:
        # mga: 20070913
        return factory.pickler()
    
    return None


# the file format registrar
def registrar():
    global _registrar
    if not _registrar:
        from Registrar import Registrar
        _registrar = Registrar()

        import native
        _registrar.manage(native, native.format())

        import chemkin
        _registrar.manage(chemkin, chemkin.format())

        import ckml
        _registrar.manage(ckml, ckml.format())

        import c
        _registrar.manage(c, c.format())

        import html
        _registrar.manage(html, html.format())

        import python
        _registrar.manage(python, python.format())

    return _registrar


# access  to the registrar
def picklers():
    return registrar()


# the registrar singleton
_registrar = None


# version
__id__ = "$Id: pickle.py,v 1.1.1.1 2007-09-13 18:17:28 aivazis Exp $"

#  End of file 
