#!/usr/bin/env python
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
from BaseParser import BaseParser
class Trans(BaseParser):
    def aTransLine(self, token):
        # dispatch to appropriate line info parser
        self._lineParsers[0](token)
        return 0
    # transitions
    def aTransSection(self, token):
        self._info.log("trans parser: section start")
        self._parse(self._scanner, self._tokenizer)
        return 0
    # other methods
    def __init__(self, mechanism, tokenizer):
        import pyre
        BaseParser.__init__(self, mechanism)
        self._tokenizer = tokenizer
        import fuego
        self._scanner = fuego.serialization.chemkin.unpickle.scanners.trans()
        # Private data
        self._transAll = 0
        self._parameters = []
        import re
        from fuego.serialization.chemkin.unpickle.tokens.RegularExpressions import species
        self._speciesScanner = re.compile(species)
        self._nextId = 0 #line ids are zero-based
        self._currentSpecies = None
        self._lineParsers = [ self._parseLine1 ]

    def _parseLine1(self, token):
        spec_tmp = token.text
        lin = token.id
        params_tmp = token.text_2
        spec = spec_tmp.strip()
        params = params_tmp.strip()
        match = self._speciesScanner.match(spec)
        if not match:
            msg = "Could not match a valid species name in '%s'" % spec
            self.onError(msg, self.locator())
        speciesName = match.group()
        species = self._mechanism.species(speciesName)
        if not species:
            msg = "trans section: undeclared species '%s'" % speciesName
            self.onWarning(msg, self.locator())
            species = self._mechanism.newSpecies(speciesName)
        species.locator(self.locator())
        self._currentSpecies = species
        self._parameters = [lin]
        self.EPS = params.split()[0]
        self.SIG = params.split()[1]
        self.DIP = params.split()[2]
        self.POL = params.split()[3]
        self.ZROT = params.split()[4]
        self._currentSpecies.transParametrization(str(lin).strip(), self.EPS, self.SIG, self. DIP, self.POL, self.ZROT, self.locator(), self._parameters )
