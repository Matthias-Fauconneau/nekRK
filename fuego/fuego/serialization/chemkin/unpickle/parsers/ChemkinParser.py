#!/usr/bin/env python
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved

from BaseParser import BaseParser

class ChemkinParser(BaseParser):
    def parse(self, mechanism, file):
        # Reset the token counts
        import pyre
        import fuego
        import journal
        from fuego.serialization.chemkin.unpickle.tokens.Token import Token

        Token._constructed = 0
        Token._destructed = 0

        self._mechanism = mechanism

        # prepare the parsing machinery
        scanner = fuego.serialization.chemkin.unpickle.scanners.sections()
        tokenizer = pyre.parsing.tokenizer(file)
        tokenizer._info = journal.debug("fuego")

        # section parsers
        from Elements import Elements
        self._elementParser = Elements(mechanism, tokenizer)

        from Species import Species
        self._speciesParser = Species(mechanism, tokenizer)

        from Thermo import Thermo
        self._thermoParser = Thermo(mechanism, tokenizer)

        #if doTrans/='n':
        from Trans import Trans
        self._transParser = Trans(mechanism, tokenizer)

        from Reactions import Reactions
        self._reactionParser = Reactions(mechanism, tokenizer)

        # enter the parsing loop
        return BaseParser.parse(self, scanner, tokenizer)


    # handlers for the section headers

    def anElementSection(self, token):
        return self._elementParser.anElementSection(token)

    def aSpeciesSection(self, token):
        return self._speciesParser.aSpeciesSection(token)

    def aThermoSection(self, token):
        return self._thermoParser.aThermoSection(token)

    def aTransSection(self, token):
        return self._transParser.aTransSection(token)

    def aReactionSection(self, token):
        return self._reactionParser.aReactionSection(token)


    # end-of-file handler

    def onEndOfFile(self):
        self._elementParser.onEndOfFile()
        self._speciesParser.onEndOfFile()
        self._thermoParser.onEndOfFile()
        #if doTrans/='n':
        self._transParser.onEndOfFile()
        self._reactionParser.onEndOfFile()
        return

    def __init__(self):

        BaseParser.__init__(self)

        # the table of declared species
        self._species = {}

        # section parsers
        self._elementParser = None
        self._speciesParser = None
        self._thermoParser = None
        self._transParser = None
        self._reactionParser = None
