#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007 All Rights Reserved

from pyre.weaver.mills.CMill import CMill

from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ
from pyre.handbook.constants.fundamental import avogadro
from pyre.handbook.constants.fundamental import gas_constant as R

smallnum = 1.0e-100
R = 8.314510 * J/mole/kelvin
Rc = 1.9872156 * cal/mole/kelvin
Patm = 1013250.0
sym  = "_"

class speciesDb:
        def __init__(self, id, name, mwt):
                self.id = id
                self.symbol = name
                self.weight = mwt
                return

class CPickler(CMill):
        def __init__(self):
                CMill.__init__(self)
                self.species = []
                self.n_species = 0
                return

        def set_species(self, mechanism):
                import pyre.handbook
                periodic = pyre.handbook.periodicTable()

                n_species = len(mechanism.species())
                self.species = [ 0.0 for x in range(n_species) ]

                for specie in mechanism.species():
                        weight = 0.0
                        for elem, coef in specie.composition:
                                aw = mechanism.element(elem).weight
                                if not aw:
                                        aw = periodic.symbol(elem.capitalize()).atomicWeight
                                weight += coef * aw

                        tempsp = speciesDb(specie.id, specie.symbol, weight)
                        self.species[specie.id] = tempsp

                self.n_species = n_species

        def _renderDocument(self, mechanism, options=None):
                self.set_species(mechanism)
                self._write('#ifdef __FG_ENABLE__')
                self._write('#define __FG_NSPECIES__ (%d)' % (len(mechanism.species())))
                inert_specie = next((i for i,s in enumerate(mechanism.species()) if s.symbol=='AR'), -1)
                if inert_specie == -1:
                    inert_specie = next((i for i,s in enumerate(mechanism.species()) if s.symbol=='N2'), -1)
                self._write('#define __FG_INERT_SPECIES__ (%d)' % inert_specie)
                self._write('#ifndef __FG_DIV__')
                self._write('#define __FG_DIV__(a,b) (a)/(b)')
                self._write('#endif')

                weight = 0.0 #species.molecularWeight()
                for elem, coef in mechanism.species()[inert_specie].composition:
                        aw = mechanism.element(elem).weight
                        if not aw:
                                import pyre.handbook
                                periodic = pyre.handbook.periodicTable()
                                aw = periodic.symbol(elem.capitalize()).atomicWeight
                        weight += coef * aw
                #self._write('const dfloat fg_rcp_molar_mass_inert_specie = %.16e;' % (1./weight))
                #self.names(mechanism)
                self.molecularWeight(mechanism)
                speciesInfo = self.analyzeThermodynamics(mechanism)
                self.gibbs(speciesInfo)
                self.cp(speciesInfo)
                self.speciesEnthalpy_RT(speciesInfo)
                self.rates(mechanism)
                #self._write('void rates(dfloat * sc, dfloat T, dfloat* gibbs0_RT, dfloat * wdot) {}') # Dummy placeholder because OCCA takes forever to parse the real thing
                self.mean_specific_heat_at_CP_R(mechanism)
                #self.concentrations(mechanism)
                self._write('#endif')

        def names(self, mechanism):
                self._write('const char* species[] = {')
                self._indent()
                for species in mechanism.species():
                        self._write('"%s",' % species.symbol)
                self._outdent()
                self._write('};')

        def concentrations(self, mechanism):
                self._write('__FG_DEVICE__ void fg_concentrations(dfloat rho, dfloat *const y, dfloat * c)')
                self._write('{')
                self._indent()
                for species in self.species:
                        self._write('c[%d] = rho * y[%d]/%.16e; ' % (species.id, species.id, species.weight) )
                self._outdent()
                self._write('}')

        def molecularWeight(self, mechanism):
                import pyre.handbook
                periodic = pyre.handbook.periodicTable()

                self._write('__FG_CONST__ dfloat fg_molar_mass[%s] = {'%(len(mechanism.species())))
                self._indent()
                for i, species in enumerate(mechanism.species()):
                    weight = 0.0 #species.molecularWeight()
                    for elem, coef in species.composition:
                            aw = mechanism.element(elem).weight
                            if not aw:
                                    aw = periodic.symbol(elem.capitalize()).atomicWeight
                            weight += coef * aw
                    self._write('%.16e ' % (weight))
                    if i != len(mechanism.species())-1: self._write(',') # OCCA forbids trailing comma
                self._outdent()
                self._write('};')
                self._write('__FG_CONST__ dfloat fg_rcp_molar_mass[%s] = {'%(len(mechanism.species())))
                self._indent()
                for i, species in enumerate(mechanism.species()):
                        weight = 0.0 #species.molecularWeight()
                        for elem, coef in species.composition:
                                        aw = mechanism.element(elem).weight
                                        if not aw:
                                                        aw = periodic.symbol(elem.capitalize()).atomicWeight
                                        weight += coef * aw
                        self._write('%.16e ' % (1./weight))
                        if i != len(mechanism.species())-1: self._write(',') # OCCA forbids trailing comma
                self._outdent()
                self._write('};')

        def rates(self, mechanism):
                n_species = len(mechanism.species())
                nReactions = len(mechanism.reaction())
                self._write('__FG_DEVICE__ void fg_rates(const dfloat *sc, const dfloat *tc, dfloat *wdot)')
                self._write('{')
                self._indent()
                self._write('dfloat qdot;')
                self.initialize_rates_calculation(mechanism)
                for reaction in mechanism.reaction():
                        self._write()
                        self._write(self.line('reaction %d: %s' % (reaction.id, reaction.equation())))
                        # compute the rates
                        self.forward_rates(mechanism, reaction)
                        self.reverse_rates(mechanism, reaction)

                        # store the progress rate
                        self._write("qdot = q_f - q_r;")

                        for symbol, coefficient in reaction.reactants:
                                #print reaction.id, mechanism.species(symbol).id
                                if coefficient==1:
                                        self._write("wdot[%d] -= qdot;" % (mechanism.species(symbol).id))
                                else:
                                        self._write("wdot[%d] -= %d * qdot;" % (mechanism.species(symbol).id, coefficient))

                        for symbol, coefficient in reaction.products:
                                #print reaction.id, mechanism.species(symbol).id
                                if coefficient==1:
                                        self._write("wdot[%d] += qdot;" % (mechanism.species(symbol).id))
                                else:
                                        self._write("wdot[%d] += %d * qdot;" % (mechanism.species(symbol).id, coefficient))

                self._write()
                self._write('return;')
                self._outdent()
                self._write('}')

        def initialize_rates_calculation(self, mechanism):
                n_species = len(mechanism.species())
                nReactions = len(mechanism.reaction())
                self._write('int id; ' + self.line('loop counter'))
                self._write('dfloat mixture;' + self.line('mixture concentration'))
                self._write('dfloat gibbs0_RT[%d];' % n_species)
                self._write('dfloat rcp_Kc;                      ' + self.line('equilibrium constant'))
                self._write('dfloat k_f;                     ' + self.line('forward reaction rate'))
                self._write('dfloat k_r;                     ' + self.line('reverse reaction rate'))
                self._write('dfloat q_f;                     ' + self.line('forward progress rate'))
                self._write('dfloat q_r;                     ' + self.line('reverse progress rate'))
                self._write('dfloat phi_f;                   ' + self.line('forward phase space factor'))
                self._write('dfloat phi_r;                   '+ self.line('reverse phase space factor'))
                self._write('dfloat alpha;                   ' + self.line('enhancement'))
                self._write('dfloat redP;                    ' + self.line('reduced pressure'))
                self._write('dfloat logPred;                 ' + self.line('log of above'))
                self._write('dfloat F;                       ' + self.line('fallof rate enhancement'))
                self._write('dfloat F_troe;                  ' + self.line('TROE intermediate'))
                self._write('dfloat logFcent;                ' + self.line('TROE intermediate'))
                self._write('dfloat troe;                    ' + self.line('TROE intermediate'))
                self._write('dfloat troe_c;                  ' + self.line('TROE intermediate'))
                self._write('dfloat troe_n;                  ' + self.line('TROE intermediate'))
                self._write('const dfloat refC = %.16e * tc[5];' % (atm.value/R.value))
                self._write('const dfloat rcp_refC = %.16e * tc[1];' % (R.value/atm.value))
                self._write('const dfloat T = tc[1];')
                self._write('mixture = 0.0;')
                self._write('for (id = 0; id < %d; ++id) {' % n_species)
                self._indent()
                self._write('mixture += sc[id];')
                self._outdent()
                self._write('}')
                self._write('gibbs_RT(gibbs0_RT, tc);')

        def forward_rates(self, mechanism, reaction):
                lt = reaction.lt
                if lt:
                        import journal
                        journal.firewall("fuego").log("Landau-Teller reactions are not supported yet")
                        return self._landau(reaction)

                dim = self.phaseSpaceUnits(reaction.reactants)

                phi_f = self.phaseSpace(mechanism, reaction.reactants, reaction)
                self._write("phi_f = %s;" % phi_f)

                arrhenius = self.arrhenius(reaction, reaction.arrhenius)

                thirdBody = reaction.thirdBody
                if not thirdBody:
                        uc = self.prefactorUnits(reaction.units["prefactor"], 1-dim)
                        self._write("k_f = %e * %s;" % (uc.value, arrhenius))
                        self._write("q_f = phi_f * k_f;")
                        return

                alpha = self.enhancement(mechanism, reaction)
                self._write("alpha = %s;" % alpha)

                sri = reaction.sri
                low = reaction.low
                troe = reaction.troe

                if not low:
                        uc = self.prefactorUnits(reaction.units["prefactor"], -dim)
                        self._write("k_f = %e * alpha * %s;" % (uc.value, arrhenius))
                        self._write("q_f = phi_f * k_f;")
                        return

                uc = self.prefactorUnits(reaction.units["prefactor"], 1-dim)
                self._write("k_f = %e * %s;" % (uc.value, arrhenius))
                k_0 = self.arrhenius(reaction, reaction.low, 0)

                redP = "__FG_DIV__(alpha, k_f) * " + k_0
                self._write("redP = 1.0e-12 * %s;" % redP)
                self._write("F = __FG_DIV__(redP, (1 + redP));")

                if sri:
                        self._write("logPred = __FG_LOG10_APPROX__(redP);")

                        self._write("X = 1.0 / (1.0 + logPred*logPred);")

                        SRI = "__FG_EXP__(X * log(%e*__FG_EXP_APPROX__(%e*tc[5]) + __FG_EXP_APPROX__(T*%e))" % (sri[0], -sri[1], 1./-sri[2])
                        if len(sri) > 3:
                                SRI += " * %e * __FG_EXP_APPROX__(%e*tc[0])" % (sri[3], sri[4])

                        self._write("F_sri = %s;" % SRI)
                        self._write("F *= Ftroe;")

                elif troe:
                        self._write("logPred = __FG_LOG10_APPROX__(redP);")

                        logF_cent = "logFcent = __FG_LOG10_APPROX__("
                        logF_cent += "(%e*__FG_EXP_APPROX__(T*(%e)))" % (1-troe[0], 1./-troe[1])
                        logF_cent += "+ (%e*__FG_EXP_APPROX__(T*(%e)))" % (troe[0], 1./-troe[2])
                        if len(troe) == 4:
                                logF_cent += "+ (__FG_EXP_APPROX__(%e*tc[5]))" % (-troe[3])
                        logF_cent += ');'
                        self._write(logF_cent)

                        d = .14
                        self._write("troe_c = -.4 - .67 * logFcent;")
                        self._write("troe_n = .75 - 1.27 * logFcent;")
                        self._write("troe = __FG_DIV__((troe_c + logPred), (troe_n - .14*(troe_c + logPred)));")
                        self._write("F_troe = __FG_POW_APPROX__(10, __FG_DIV__(logFcent, (1.0 + troe*troe)));")
                        self._write("F *= F_troe;")


                self._write("k_f *= F;")
                self._write("q_f = phi_f * k_f;")
                return


        def reverse_rates(self, mechanism, reaction):
                if not reaction.reversible:
                        self._write("q_r = 0.0;")
                        return

                phi_r = self.phaseSpace(mechanism, reaction.products, reaction)
                self._write("phi_r = %s;" % phi_r)

                if reaction.rlt:
                        import journal
                        journal.firewall("fuego").log("Landau-Teller reactions are not supported yet")
                        return

                if reaction.rev:
                        arrhenius = self.arrhenius(reaction, reaction.rev)
                        thirdBody = reaction.thirdBody
                        if thirdBody:
                                uc = self.prefactorUnits(reaction.units["prefactor"], -dim)
                                self._write("k_r = %e * alpha * %s;" % (uc.value, arrhenius))
                        else:
                                uc = self.prefactorUnits(reaction.units["prefactor"], 1-dim)
                                self._write("k_r = %e * %s;" % (uc.value, arrhenius))

                        self._write("q_f = phi_r * k_r;")
                        return

                rcp_Kc = self.rcp_Kc(mechanism, reaction)
                self._write("rcp_Kc = %s;" % rcp_Kc)

                self._write("k_r = k_f * rcp_Kc;")
                self._write("q_r = phi_r * k_r;")

                return


        def arrhenius(self, reaction, parameters, customExp = 1):
                A, beta, E = parameters
                if A == 0:
                        return "0.0"
                expr = "%e" % A
                if beta == 0 and E == 0:
                        return expr

                if customExp == 1:
                  expr +="*__FG_EXP__("
                else:
                  expr +="*exp("

                if beta != 0:
                        expr += "%e*tc[0]" % beta
                if E != 0:
                        uc = self.activationEnergyUnits(reaction.units["activation"])
                        expr += "%+e*tc[5]" % (- uc * E / Rc / kelvin) # catch unit conversion errors!
                expr += ')'

                return expr


        def prefactorUnits(self, code, exponent):

                if code == "mole/cm**3":
                        units = mole / cm**3
                elif code == "moles":
                        units = mole / cm**3
                elif code == "molecules":
                        from pyre.hadbook.constants.fundamental import avogadro
                        units = 1.0 / avogadro / cm**3
                else:
                        import journal
                        journal.firewall("fuego").log("unknown prefactor units '%s'" % code)
                        return 1

                return units ** exponent / second


        def activationEnergyUnits(self, code):
                if code == "cal/mole":
                        units = cal / mole
                elif code == "kcal/mole":
                        units = kcal /mole
                elif code == "joules/mole":
                        units = J / mole
                elif code == "kjoules/mole":
                        units = kJ / mole
                elif code == "kelvins":
                        units = Rc * kelvin
                else:
                        import journal
                        journal.firewall("fuego").log("unknown activation energy units '%s'" % code)
                        return 1

                return units

        def phaseSpace(self, mechanism, reagents, reaction):

                phi = []

                for symbol, coefficient in reagents:
                        #print reaction.id, mechanism.species(symbol).id
                        conc = "sc[%d]" % mechanism.species(symbol).id
                        phi += [conc] * coefficient

                return "*".join(phi)


        def phaseSpaceUnits(self, reagents):
                dim = 0
                for symbol, coefficient in reagents:
                        dim += coefficient

                return dim

        def enhancement(self, mechanism, reaction):
                thirdBody = reaction.thirdBody
                if not thirdBody:
                        import journal
                        journal.firewall("fuego").log("enhancement called for a reaction without a third body")
                        return

                species, coefficient = thirdBody
                efficiencies = reaction.efficiencies

                if not efficiencies:
                        if species == "<mixture>":
                                return "mixture"
                        #print reaction.id, mechanism.species(species).id
                        return "sc[%d]" % mechanism.species(species).id

                alpha = ["mixture"]
                for symbol, efficiency in efficiencies:
                        factor = efficiency - 1
                        #print reaction.id, mechanism.species(symbol).id
                        conc = "sc[%d]" % mechanism.species(symbol).id
                        if factor == 1:
                                alpha.append(conc)
                        else:
                                alpha.append("%e*%s" % (factor, conc))

                return " + ".join(alpha)

        def cp(self, speciesInfo):
                self._write(self.line('compute Cp/R at the given temperature'))
                self._write(self.line('tc contains precomputed powers of T, tc[0] = log(T)'))
                self.generateThermoRoutine("cp_R", self.cpNASA, speciesInfo)

        def gibbs(self, speciesInfo):
                self.generateThermoRoutine("gibbs_RT", self.gibbsNASA, speciesInfo)

        def speciesEnthalpy_RT(self, speciesInfo):
                self.generateThermoRoutine("fg_speciesEnthalpy_RT", self.enthalpyNASA, speciesInfo)

        def generateThermoRoutine(self, name, expressionGenerator, speciesInfo):
                lowT, highT, midpoints = speciesInfo
                self._write('__FG_DEVICE__ void %s(dfloat * species, const dfloat tc[])' % name)
                self._write('{')
                self._indent()
                self._write('dfloat T = tc[1];')
                for midT, speciesList in midpoints.items():
                        self._write(self.line('species with midpoint at T=%e kelvin' % midT))
                        self._write('if (T < %e) {' % midT)
                        self._indent()
                        for species, lowRange, highRange in speciesList:
                                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                                self._write('species[%d] =' % species.id)
                                self._indent()
                                expressionGenerator(lowRange.parameters)
                                self._outdent()
                        self._outdent()
                        self._write('} else {')
                        self._indent()
                        for species, lowRange, highRange in speciesList:
                                self._write(self.line('species %d: %s' % (species.id, species.symbol)))
                                self._write('species[%d] =' % species.id)
                                self._indent()
                                expressionGenerator(highRange.parameters)
                                self._outdent()
                        self._outdent()
                        self._write('}')
                self._outdent()
                self._write('}')

        def analyzeThermodynamics(self, mechanism):
                lowT = 0.0
                highT = 1000000.0

                midpoints = {}

                for species in mechanism.species():
                        models = species.thermo
                        if len(models) > 2:
                                import journal
                                print models
                                journal.firewall("fuego").log(
                                        "species '%s' has more than two thermo regions" % species.symbol)
                                return

                        m1 = models[0]
                        m2 = models[1]

                        if m1.lowT < m2.lowT:
                                lowRange = m1
                                highRange = m2
                        else:
                                lowRange = m2
                                highRange = m1

                        low = lowRange.lowT
                        mid = lowRange.highT
                        high = highRange.highT

                        if low > lowT:
                                lowT = low
                        if high < highT:
                                highT = high

                        midpoints.setdefault(mid, []).append((species, lowRange, highRange))

                return lowT, highT, midpoints

        def rcp_Kc(self, mechanism, reaction):
                dim = 0
                dG = ""
                terms = []
                for symbol, coefficient in reaction.reactants:
                        if coefficient == 1:
                                factor = ""
                        else:
                                factor = "%d * " % coefficient

                        #print reaction.id, mechanism.species(symbol).id
                        terms.append("%sgibbs0_RT[%d]" % (factor, mechanism.species(symbol).id))
                        dim -= coefficient
                dG += '(' + ' + '.join(terms) + ')'

                # flip the signs
                terms = []
                for symbol, coefficient in reaction.products:
                        if coefficient == 1:
                                factor = ""
                        else:
                                factor = "%d * " % coefficient
                        #print reaction.id, mechanism.species(symbol).id
                        terms.append("%sgibbs0_RT[%d]" % (factor, mechanism.species(symbol).id))
                        dim += coefficient
                dG += ' - (' + ' + '.join(terms) + ')'
                rcp_Kp = '__FG_EXP__(-(' + dG + '))'
                if dim == 0:
                        conversion = ""
                elif dim > 0:
                        conversion = "*".join(["rcp_refC"] * abs(dim)) + ' * '
                else:
                        conversion = "*".join(["refC"] * abs(dim)) + ' * '
                return conversion + rcp_Kp

        def cpNASA(self, parameters):
                self._write('%+15.8e' % parameters[0])
                self._write('%+15.8e * tc[1]' % parameters[1])
                self._write('%+15.8e * tc[2]' % parameters[2])
                self._write('%+15.8e * tc[3]' % parameters[3])
                self._write('%+15.8e * tc[4];' % parameters[4])
                return

        def gibbsNASA(self, parameters):
                self._write('%+15.8e * tc[5]' % parameters[5])
                self._write('%+15.8e' % (parameters[0] - parameters[6]))
                self._write('%+15.8e * tc[0]' % (-parameters[0]))
                self._write('%+15.8e * tc[1]' % (-parameters[1]/2))
                self._write('%+15.8e * tc[2]' % (-parameters[2]/6))
                self._write('%+15.8e * tc[3]' % (-parameters[3]/12))
                self._write('%+15.8e * tc[4];' % (-parameters[4]/20))
                return

        def enthalpyNASA(self, parameters):
                self._write('%+15.8e' % parameters[0])
                self._write('%+15.8e * tc[1]' % (parameters[1]/2))
                self._write('%+15.8e * tc[2]' % (parameters[2]/3))
                self._write('%+15.8e * tc[3]' % (parameters[3]/4))
                self._write('%+15.8e * tc[4]' % (parameters[4]/5))
                self._write('%+15.8e * tc[5];' % (parameters[5]))
                return

        def mean_specific_heat_at_CP_R(self, mechanism):
                self._write('__FG_DEVICE__ double fg_mean_specific_heat_at_CP_R(double T, const double* mole_fractions)')
                self._write('{')
                self._indent()
                #self._write('dfloat tc[] = { log(T), T, T*T, T*T*T, T*T*T*T, 1./T }; ')
                self._write('dfloat tc[6]; ')
                self._write('tc[0] = log(T); ')
                self._write('tc[1] = T; ')
                self._write('tc[2] = T*T; ')
                self._write('tc[3] = T*T*T; ')
                self._write('tc[4] = T*T*T*T; ')
                self._write('tc[5] = 1./T; ')
                self._write('dfloat result = 0; ')
                self._write('dfloat cpor[%d]; ' % self.n_species)
                self._write('cp_R(cpor, tc);')
                for species in self.species:
                        self._write('result += cpor[%d]*mole_fractions[%d]; ' % (species.id, species.id))
                self._write('return result;')
                self._outdent()
                self._write('}')
