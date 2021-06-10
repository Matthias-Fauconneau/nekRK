#!/usr/bin/env python
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved

def save(mechanism, format="chemkin"):
    factory = registrar().retrieve(format)
    if not factory:
        return []
    pickler = factory.pickler()
    pickler.initialize()
    return pickler.pickle(mechanism)

# factory methods for the serializers

def pickler(format="chemkin"):
    factory = registrar().retrieve(format)
    if factory:
        return factory()
    return None

# the file format registrar
def registrar():
    global _registrar
    if not _registrar:
        from Registrar import Registrar
        _registrar = Registrar()
        import chemkin
        _registrar.register(chemkin, chemkin.format())
        import c
        _registrar.register(c, c.format())
    return _registrar

# access  to the registrar
def picklers():
    return registrar()

# the registrar singleton
_registrar = None
