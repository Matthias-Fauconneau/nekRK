# coding=utf8
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007 All Rights Reserved

from collections import OrderedDict
from weaver.mills.CMill import CMill
from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ
from pyre.handbook.constants.fundamental import avogadro
from pyre.handbook.constants.fundamental import gas_constant as R
import numpy as np

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
                self.lowT = 100.0
                self.highT = 10000.0

        def set_species(self, mechanism):
                import pyre.handbook
                periodic = pyre.handbook.periodicTable()

                n_species = len(mechanism.species())
                self.species = [ 0.0 for x in range(n_species) ]

                for specie in mechanism.species():
                        weight = 0.0
                        for element, number_of_atoms in specie.composition:
                                atomic_weight = mechanism.element(element).weight #?
                                if not atomic_weight:
                                        atomic_weight = periodic.symbol(element.capitalize()).atomicWeight # handbook
                                weight += number_of_atoms * atomic_weight

                        tempsp = speciesDb(specie.id, specie.symbol, weight)
                        self.species[specie.id] = tempsp

                self.n_species = n_species

        def _renderDocument(self, mechanism, options=None):
                self.set_species(mechanism)
                self._write('#define n_species %d' % (len(mechanism.species())))
                inert_specie = next((i for i,s in enumerate(mechanism.species()) if s.symbol=='AR'), -1)
                if inert_specie == -1:
                    inert_specie = next((i for i,s in enumerate(mechanism.species()) if s.symbol=='N2'), -1)
                self._write('const int inert_specie= %d;' % inert_specie)

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
                speciesTransport = self._analyzeTransport(mechanism)
                self.viscosity_and_thermal_conductivity(mechanism, speciesTransport, False, NTFit=50)
                self.binary_diffusion_coefficients(speciesTransport, False, NTFit=50)

        def names(self, mechanism):
                self._write('const char* species[] = {')
                self._indent()
                for species in mechanism.species():
                        self._write('"%s",' % species.symbol)
                self._outdent()
                self._write('};')

        def concentrations(self, mechanism):
                self._write('void fg_concentrations(dfloat rho, dfloat *const y, dfloat * c)')
                self._write('{')
                self._indent()
                for species in self.species:
                        self._write('c[%d] = rho * y[%d]/%.16e; ' % (species.id, species.id, species.weight) )
                self._outdent()
                self._write('}')

        def molecularWeight(self, mechanism):
                import pyre.handbook
                periodic = pyre.handbook.periodicTable()

                self._write('const dfloat fg_molar_mass[%s] = {'%(len(mechanism.species())))
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
                self._write('const dfloat fg_rcp_molar_mass[%s] = {'%(len(mechanism.species())))
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
                self._write('void fg_rates(dfloat * sc, dfloat tc[], dfloat * wdot)')
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
                                if coefficient==1:
                                        self._write("wdot[%d] -= qdot;" % (mechanism.species(symbol).id))
                                else:
                                        self._write("wdot[%d] -= %d * qdot;" % (mechanism.species(symbol).id, coefficient))

                        for symbol, coefficient in reaction.products:
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
                self._write('const dfloat refC = (%.16e / %.16e) * tc[5];' % (atm.value, R.value))
                self._write('const dfloat rcp_refC = 1/refC;')
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

                phi_f = self.phaseSpace(mechanism, reaction.reactants)
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
                k_0 = self.arrhenius(reaction, reaction.low)
                redP = "alpha / k_f * " + k_0
                self._write("redP = 1.0e-12 * %s;" % redP)
                self._write("F = redP / (1 + redP);")

                if sri:
                        self._write("logPred = log10(redP);")

                        self._write("X = 1.0 / (1.0 + logPred*logPred);")

                        SRI = "fgexp(X * log(%e*fgexp(%e*tc[5]) + fgexp(T*%e))" % (sri[0], -sri[1], 1./-sri[2])
                        if len(sri) > 3:
                                SRI += " * %e * fgexp(%e*tc[0])" % (sri[3], sri[4])

                        self._write("F_sri = %s;" % SRI)
                        self._write("F *= Ftroe;")

                elif troe:
                        self._write("logPred = log10(redP);")

                        logF_cent = "logFcent = log10("
                        logF_cent += "(%e*fgexp(T*(%e)))" % (1-troe[0], 1./-troe[1])
                        logF_cent += "+ (%e*fgexp(T*(%e)))" % (troe[0], 1./-troe[2])
                        if len(troe) == 4:
                                logF_cent += "+ (fgexp(%e*tc[5]))" % (-troe[3])
                        logF_cent += ');'
                        self._write(logF_cent)

                        d = .14
                        self._write("troe_c = -.4 - .67 * logFcent;")
                        self._write("troe_n = .75 - 1.27 * logFcent;")
                        self._write("troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));")
                        self._write("F_troe = pow(10, logFcent / (1.0 + troe*troe));")
                        self._write("F *= F_troe;")

                self._write("k_f *= F;")
                self._write("q_f = phi_f * k_f;")
                return


        def reverse_rates(self, mechanism, reaction):
                if not reaction.reversible:
                        self._write("q_r = 0.0;")
                        return

                phi_r = self.phaseSpace(mechanism, reaction.products)
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


        def arrhenius(self, reaction, parameters):
                A, beta, E = parameters
                if A == 0:
                        return "0.0"
                expr = "%e" % A
                if beta == 0 and E == 0:
                        return expr
                expr +="*fgexp("
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

        def phaseSpace(self, mechanism, reagents):
                phi = []
                for symbol, coefficient in reagents:
                        conc = "sc[%d]" % mechanism.species(symbol).id
                        phi += [conc] * int(coefficient)
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
                        return "sc[%d]" % mechanism.species(species).id

                alpha = ["mixture"]
                for symbol, efficiency in efficiencies:
                        factor = efficiency - 1
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
                self._write('void %s(dfloat * species, dfloat tc[])' % name)
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
                        terms.append("%sgibbs0_RT[%d]" % (factor, mechanism.species(symbol).id))
                        dim += coefficient
                dG += ' - (' + ' + '.join(terms) + ')'
                rcp_Kp = 'fgexp(-(' + dG + '))'
                if dim == 0:
                        conversion = ""
                elif dim > 0:
                        conversion = "*".join(["rcp_refC"] * abs(int(dim))) + ' * '
                else:
                        conversion = "*".join(["refC"] * abs(int(dim))) + ' * '
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
                self._write('dfloat fg_mean_specific_heat_at_CP_R(dfloat T, const dfloat* mole_fractions)')
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

        def _analyzeTransport(self, mechanism):
                transdata = OrderedDict()
                for specie in mechanism.species():
                        models = specie.trans
                        if len(models) > 2:
                                print 'species: ', specie
                                import pyre
                                pyre.debug.Firewall.hit("unsupported configuration in species.trans")
                                return
                        print(specie)
                        m1 = models[0]
                        lin = m1.parameters[0]
                        eps = m1.eps
                        sig = m1.sig
                        dip = m1.dip
                        pol = m1.pol
                        zrot = m1.zrot
                        transdata[specie] = [lin, eps, sig, dip, pol, zrot]
                return transdata

        def viscosity_and_thermal_conductivity(self, mechanism, speciesTransport, do_declarations, NTFit):
                #compute single constants in g/cm/s
                kb = 1.3806503e-16
                Na = 6.02214199e23
                RU = 8.31447e7
                #temperature increment
                dt = (self.highT-self.lowT) / (NTFit-1)
                #factor dependent upon the molecule
                m_crot = np.zeros(self.n_species)
                m_cvib = np.zeros(self.n_species)
                isatm = np.zeros(self.n_species)
                def for_():
                    for spec in speciesTransport:
                            if int(speciesTransport[spec][0]) == 0:
                                    m_crot[spec.id] = 0.0
                                    m_cvib[spec.id] = 0.0
                                    isatm[spec.id] = 0.0
                            elif int(speciesTransport[spec][0]) == 1:
                                    m_crot[spec.id] = 1.0
                                    m_cvib[spec.id] = 5.0 / 2.0
                                    isatm[spec.id] = 1.0
                            else:
                                    m_crot[spec.id] = 1.5
                                    m_cvib[spec.id] = 3.0
                                    isatm[spec.id] = 1.0
                for_()
                #viscosities coefs (4 per spec)
                ln_viscosity = OrderedDict()
                #conductivities coefs (4 per spec)
                ln_thermal_conductivity = OrderedDict()
                for transport_specie in speciesTransport:
                        spvisc = []
                        spcond = []
                        tlog = []
                        for n in range(NTFit):
                                t = self.lowT + dt*n
                                #variables
                                #eq. (2)
                                interaction_well_depth = float(speciesTransport[transport_specie][1])
                                diameter_Aengstroem = float(speciesTransport[transport_specie][2]) #Ångström
                                tr = t/interaction_well_depth
                                diameter = diameter_Aengstroem*1e-10
                                dst = (1./2 / kb#Why?
                                                * float(speciesTransport[transport_specie][3]#What is this?
                                        )**2 / (interaction_well_depth * diameter**3)) #What is this?
                                #viscosity of specie at t
                                #eq. (1)
                                molar_mass = self.species[transport_specie.id].weight
                                visc = (5.0 / 16.0) * np.sqrt(np.pi * molar_mass * kb * t / Na) / (self.om22_CHEMKIN(tr,dst) * np.pi * diameter**2)
                                #conductivity of spec at t
                                #eq. (30)
                                m_red = molar_mass / (2.0 * Na)
                                diffcoef = (3.0 / 16.0) * np.sqrt(2.0 * np.pi * kb**3 * t**3 / m_red) / (10.0 * np.pi * self.om11_CHEMKIN(tr,dst) * diameter**2)
                                #eq. (19)
                                cv_vib_R = (self._getCVdRspecies(mechanism, t, transport_specie) - m_cvib[transport_specie.id]) * isatm[transport_specie.id]
                                rho_atm = 10.0 * molar_mass /(RU * t)
                                f_vib = rho_atm * diffcoef / visc
                                #eq. (20)
                                A = 2.5 - f_vib
                                #eqs. (21) + (32-33)
                                cv_rot_R = m_crot[transport_specie.id]
                                #note: the T corr is not applied in CANTERA
                                B = (float(speciesTransport[transport_specie][5]) \
                                                * self.Fcorr(298.0, float(speciesTransport[transport_specie][1])) / self.Fcorr(t, float(speciesTransport[transport_specie][1])) \
                                                + (2.0 / np.pi) * ((5.0 / 3.0 ) * cv_rot_R  + f_vib))
                                #eq. (18)
                                f_rot = f_vib * (1.0 + 2.0 / np.pi * A / B )
                                #eq. (17)
                                cv_trans_R = 3.0 / 2.0
                                f_trans = 5.0 / 2.0 * (1.0 - 2.0 / np.pi * A / B * cv_rot_R / cv_trans_R )
                                if (int(speciesTransport[transport_specie][0]) == 0):
                                        cond = (visc * RU / molar_mass) * \
                                                        (5.0 / 2.0) * cv_trans_R
                                else:
                                        cond = (visc * RU / molar_mass) * \
                                                (f_trans * cv_trans_R + f_rot * cv_rot_R + f_vib * cv_vib_R)

                                #log transformation for polyfit
                                tlog.append(np.log(t))
                                spvisc.append(np.log(visc))
                                spcond.append(np.log(cond))

                        ln_viscosity[transport_specie.id] = np.polyfit(tlog, spvisc, 3) # log viscosity = P(log T)
                        ln_thermal_conductivity[transport_specie.id] = np.polyfit(tlog, spcond, 3)

                self._write('void fg_viscosity_and_thermal_conductivity(dfloat _p, dfloat T, const dfloat mass_fractions[], /*->*/ dfloat& viscosity, dfloat& thermal_conductivity) {')
                self._indent()
                self._write('dfloat mean_rcp_molar_mass = 0.;')
                self._write('for(int i=0; i<n_species; i++) {')
                self._write('mean_rcp_molar_mass += mass_fractions[i]*fg_rcp_molar_mass[i];')
                self._write('}')
                self._write('dfloat mean_molar_mass = 1./mean_rcp_molar_mass;')
                self._write('dfloat mole_fractions[n_species];')
                self._write('for(int i=0; i<n_species; i++) {')
                self._write('mole_fractions[i] = mass_fractions[i]*fg_rcp_molar_mass[i]*mean_molar_mass;')
                self._write('}')

                def evaluate_polynomial(P, x):
                    self._write('dfloat y = 0.;')
                    for i in range(4):
                            self._write('y += %.8E*pow(%s,%d);' % (P[3-i], x, i)) #?
                self._write('dfloat ln_T = log(T);')

                # Viscosity
                self._write('{')
                self._indent()
                self._write('dfloat sum = 0.;')
                for spec in self.species: #i|
                        i = spec.id
                        self._write('{')
                        self._indent()
                        #{
                        y = evaluate_polynomial(ln_thermal_conductivity[i], "ln_T")
                        self._write('dfloat viscosity_i = exp(y);')
                        self._write('sum += mole_fractions[%d] * pow(viscosity_i, 6.);'%(i))
                        #}
                        self._outdent()
                        self._write('}')
                #}
                self._write('viscosity = pow(sum, 1./6.);')
                self._outdent()
                self._write('}')

                # Thermal conductivity
                self._write('{')
                self._indent()
                self._write('dfloat sum = 0.;')
                for spec in self.species: #i|
                        i = spec.id
                        self._write('{')
                        self._indent()
                        #{
                        y = evaluate_polynomial(ln_thermal_conductivity[i], "ln_T")
                        self._write('dfloat thermal_conductivity_i = exp(y);')
                        self._write('sum += mole_fractions[%d] * pow(thermal_conductivity_i, 4.);'%(i)) #?
                        #}
                        self._outdent()
                        self._write('}')
                #}
                self._write('thermal_conductivity = pow(sum, 1./4.);') #?
                self._outdent()
                self._write('}')

                self._outdent()
                self._write('}')

        def binary_diffusion_coefficients(self, speciesTransport, do_declarations, NTFit) :
                #REORDERING OF SPECS
                specOrdered = []
                for i in range(self.n_species):
                        for spec in speciesTransport:
                                if spec.id == i:
                                        specOrdered.append(spec)
                                        break

                #compute single constants in g/cm/s
                kb = 1.3806503e-16
                Na = 6.02214199e23
                atmospheric_pressure = 101325.
                #temperature increment
                dt = (self.highT-self.lowT) / (NTFit-1)
                #diff coefs (4 per spec pair)
                binary_diffusion_coefficients = []
                for i, transport_specie_i in enumerate(specOrdered):
                        binary_diffusion_coefficients.append([])
                        if (i != transport_specie_i.id):
                                print "Problem in _diffcoefs computation"
                                stop
                        for j, transport_specie_j in enumerate(specOrdered[0:i+1]):
                                if (j != transport_specie_j.id):
                                        print "Problem in _diffcoefs computation"
                                        stop
                                diameter_Aengstroem_i = float(speciesTransport[transport_specie_i][2]) #Ångström
                                diameter_i = diameter_Aengstroem_i*1e-10
                                diameter_Aengstroem_j = float(speciesTransport[transport_specie_j][2]) #Ångström
                                diameter_j = diameter_Aengstroem_j*1e-10
                                sigma = (diameter_i + diameter_j)/2. * self.Xi(transport_specie_i, transport_specie_j, speciesTransport)**(1./6.) #What is this ?
                                molar_mass_i = self.species[transport_specie_i.id].weight
                                molar_mass_j = self.species[transport_specie_j.id].weight
                                m_red = molar_mass_i * molar_mass_j / (molar_mass_i + molar_mass_j) / Na #What is this ?
                                epsm_k = np.sqrt(float(speciesTransport[transport_specie_i][1]) * float(speciesTransport[transport_specie_j][1])) * self.Xi(transport_specie_i, transport_specie_j, speciesTransport)**2.0 #What is this ?
                                dst = float(speciesTransport[transport_specie_i][3]) * float(speciesTransport[transport_specie_j][3]) / (epsm_k * sigma**3) #What is this ?
                                if self.Xi_bool(transport_specie_i, transport_specie_j, speciesTransport)==False:
                                        dst = 0.0
                                #enter the loop on temperature
                                spdiffcoef = []
                                tlog = []
                                for n in range(NTFit):
                                    t = self.lowT + dt*n
                                    tr = t/ epsm_k
                                    #eq. (3)
                                    #note: these are "corrected" in CHEMKIN not in CANTERA... we chose not to
                                    difcoeff = 3.0 / 16.0 * 1 / atmospheric_pressure * (np.sqrt(2.0 * np.pi * t**3 * kb**3 / m_red) / ( np.pi * sigma * sigma * self.om11_CHEMKIN(tr,dst)))

                                    #log transformation for polyfit
                                    tlog.append(np.log(t))
                                    spdiffcoef.append(np.log(difcoeff))

                                binary_diffusion_coefficients[i].append(np.polyfit(tlog, spdiffcoef, 3))

                #use the symmetry for upper triangular terms
                #note: starting with this would be preferable (only one bigger loop)
                #note2: or write stuff differently !
                #for i,spec1 in enumerate(specOrdered):
                #    for j,spec2 in enumerate(specOrdered[i+1:]):
                #        binary_diffusion_coefficients[i].append(binary_diffusion_coefficients[spec2.id][spec1.id])

                # No idea what this code does
                self._write('void fg_Pele_Ddiag(const dfloat wbar, const dfloat Xloc[n_species], const dfloat Yloc[n_species], dfloat logT[3], dfloat* Ddiag) {')
                self._indent()
                for i,spec1 in enumerate(specOrdered):
                    self._write('{')
                    self._indent()
                    self._write('dfloat term1 = 0.0, term2 = 0.0;')
                    for j,spec2 in enumerate(specOrdered[0:i]):
                        self._write('{')
                        self._indent()
                        self._write('dfloat dbintemp = %.8E + %.8E * logT[0] + %.8E * logT[1] + %.8E * logT[2];' % (# Why is the coefficient order being reversed here ?
                            binary_diffusion_coefficients[i][j][3], binary_diffusion_coefficients[i][j][2], binary_diffusion_coefficients[i][j][1], binary_diffusion_coefficients[i][j][0]))
                        self._write('term1 += Yloc[j];');
                        self._write('term2 += Xloc[j] / exp(dbintemp);')
                        self._write('}')
                        self._outdent()
                    self._write('Ddiag[%d] = fg_molar_mass[%d] * term1 / term2 / wbar;' % (i,i))
                    self._write('}')
                    self._outdent()
                self._outdent()
                self._write('}')
                return

        def astar(self, tslog):
                aTab = [.1106910525E+01, -.7065517161E-02,-.1671975393E-01,
                                .1188708609E-01,  .7569367323E-03,-.1313998345E-02,
                                .1720853282E-03]

                B = aTab[6]
                for i in range(6):
                        B = aTab[5-i] + B*tslog

                return B

        def bstar(self, tslog):
                bTab = [.1199673577E+01, -.1140928763E+00,-.2147636665E-02,
                                .2512965407E-01, -.3030372973E-02,-.1445009039E-02,
                                .2492954809E-03]

                B = bTab[6]
                for i in range(6):
                        B = bTab[5-i] + B*tslog

                return B

        def cstar(self, tslog):
                cTab = [ .8386993788E+00,  .4748325276E-01, .3250097527E-01,
                                -.1625859588E-01, -.2260153363E-02, .1844922811E-02,
                                -.2115417788E-03]

                B = cTab[6]
                for i in range(6):
                        B = cTab[5-i] + B*tslog

                return B

        def Xi(self, spec1, spec2, speciesTransport):
                dipmin = 1e-20
                #1 is polar, 2 is nonpolar
                #err in eq. (11) ?
                if (float(speciesTransport[spec2][3]) < dipmin) and (float(speciesTransport[spec1][3]) > dipmin):
                        xi = 1.0 + 1.0/4.0 * self.redPol(spec2, speciesTransport)*self.redDip(spec1, speciesTransport) *\
                                        self.redDip(spec1, speciesTransport) *\
                                        np.sqrt(float(speciesTransport[spec1][1]) / float(speciesTransport[spec2][1]))
                #1 is nonpolar, 2 is polar
                elif (float(speciesTransport[spec2][3]) > dipmin) and (float(speciesTransport[spec1][3]) < dipmin):
                        xi = 1.0 + 1.0/4.0 * self.redPol(spec1, speciesTransport)*self.redDip(spec2, speciesTransport) *\
                                        self.redDip(spec2, speciesTransport) *\
                                        np.sqrt(float(speciesTransport[spec2][1]) / float(speciesTransport[spec1][1]))
                #normal case, either both polar or both nonpolar
                else:
                        xi = 1.0

                return xi

        def Xi_bool(self, spec1, spec2, speciesTransport):
                dipmin = 1e-20
                #1 is polar, 2 is nonpolar
                #err in eq. (11) ?
                if (float(speciesTransport[spec2][3]) < dipmin) and (float(speciesTransport[spec1][3]) > dipmin):
                        xi_b = False
                #1 is nonpolar, 2 is polar
                elif (float(speciesTransport[spec2][3]) > dipmin) and (float(speciesTransport[spec1][3]) < dipmin):
                        xi_b = False
                #normal case, either both polar or both nonpolar
                else:
                        xi_b = True

                return xi_b

        def redPol(self, spec, speciesTransport):
                return (float(speciesTransport[spec][4]) / float(speciesTransport[spec][2])**3.0)


        def redDip(self, spec, speciesTransport):
                #compute single constants in g/cm/s
                kb = 1.3806503e-16
                #conversion coefs
                AtoCM = 1.0e-8
                DEBYEtoCGS = 1.0e-18
                convert = DEBYEtoCGS / np.sqrt( kb * AtoCM**3.0 )
                return convert * float(speciesTransport[spec][3]) / \
                                np.sqrt(float(speciesTransport[spec][1]) * float(speciesTransport[spec][2])**3.0)

        def Fcorr(self, t, eps_k):
                thtwo = 3.0 / 2.0
                return 1 + np.pi**(thtwo) / 2.0 * np.sqrt(eps_k / t) + \
                                (np.pi**2 / 4.0 + 2.0) * (eps_k / t) + \
                                (np.pi * eps_k / t)**(thtwo)

        def om11(self, tr, dst):
                # This is an overhaul of CANTERA version 2.3
                #range of dst
                dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

                #range of tr
                trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

                #tab of astar corresp. to (tr, dst)
                #CANTERA
                astarTab = [1.0065, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840,
                                1.0231, 1.0660, 1.0380, 1.0400, 1.0430, 1.0500, 1.0520, 1.0510,
                                1.0424, 1.0450, 1.0480, 1.0520, 1.0560, 1.0650, 1.0660, 1.0640,
                                1.0719, 1.0670, 1.0600, 1.0550, 1.0580, 1.0680, 1.0710, 1.0710,
                                1.0936, 1.0870, 1.0770, 1.0690, 1.0680, 1.0750, 1.0780, 1.0780,
                                1.1053, 1.0980, 1.0880, 1.0800, 1.0780, 1.0820, 1.0840, 1.0840,
                                1.1104, 1.1040, 1.0960, 1.0890, 1.0860, 1.0890, 1.0900, 1.0900,
                                1.1114, 1.1070, 1.1000, 1.0950, 1.0930, 1.0950, 1.0960, 1.0950,
                                1.1104, 1.1070, 1.1020, 1.0990, 1.0980, 1.1000, 1.1000, 1.0990,
                                1.1086, 1.1060, 1.1020, 1.1010, 1.1010, 1.1050, 1.1050, 1.1040,
                                1.1063, 1.1040, 1.1030, 1.1030, 1.1040, 1.1080, 1.1090, 1.1080,
                                1.1020, 1.1020, 1.1030, 1.1050, 1.1070, 1.1120, 1.1150, 1.1150,
                                1.0985, 1.0990, 1.1010, 1.1040, 1.1080, 1.1150, 1.1190, 1.1200,
                                1.0960, 1.0960, 1.0990, 1.1030, 1.1080, 1.1160, 1.1210, 1.1240,
                                1.0943, 1.0950, 1.0990, 1.1020, 1.1080, 1.1170, 1.1230, 1.1260,
                                1.0934, 1.0940, 1.0970, 1.1020, 1.1070, 1.1160, 1.1230, 1.1280,
                                1.0926, 1.0940, 1.0970, 1.0990, 1.1050, 1.1150, 1.1230, 1.1300,
                                1.0934, 1.0950, 1.0970, 1.0990, 1.1040, 1.1130, 1.1220, 1.1290,
                                1.0948, 1.0960, 1.0980, 1.1000, 1.1030, 1.1120, 1.1190, 1.1270,
                                1.0965, 1.0970, 1.0990, 1.1010, 1.1040, 1.1100, 1.1180, 1.1260,
                                1.0997, 1.1000, 1.1010, 1.1020, 1.1050, 1.1100, 1.1160, 1.1230,
                                1.1025, 1.1030, 1.1040, 1.1050, 1.1060, 1.1100, 1.1150, 1.1210,
                                1.1050, 1.1050, 1.1060, 1.1070, 1.1080, 1.1110, 1.1150, 1.1200,
                                1.1072, 1.1070, 1.1080, 1.1080, 1.1090, 1.1120, 1.1150, 1.1190,
                                1.1091, 1.1090, 1.1090, 1.1100, 1.1110, 1.1130, 1.1150, 1.1190,
                                1.1107, 1.1110, 1.1110, 1.1110, 1.1120, 1.1140, 1.1160, 1.1190,
                                1.1133, 1.1140, 1.1130, 1.1140, 1.1140, 1.1150, 1.1170, 1.1190,
                                1.1154, 1.1150, 1.1160, 1.1160, 1.1160, 1.1170, 1.1180, 1.1200,
                                1.1172, 1.1170, 1.1170, 1.1180, 1.1180, 1.1180, 1.1190, 1.1200,
                                1.1186, 1.1190, 1.1190, 1.1190, 1.1190, 1.1190, 1.1200, 1.1210,
                                1.1199, 1.1200, 1.1200, 1.1200, 1.1200, 1.1210, 1.1210, 1.1220,
                                1.1223, 1.1220, 1.1220, 1.1220, 1.1220, 1.1230, 1.1230, 1.1240,
                                1.1243, 1.1240, 1.1240, 1.1240, 1.1240, 1.1240, 1.1250, 1.1250,
                                1.1259, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260,
                                1.1273, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1280,
                                1.1297, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1290,
                                1.1339, 1.1340, 1.1340, 1.1350, 1.1350, 1.1340, 1.1340, 1.1320,
                                1.1364, 1.1370, 1.1370, 1.1380, 1.1390, 1.1380, 1.1370, 1.1350,
                                1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187,
                                1.14187]


                #Find for each fixed tr the poly of deg 6 in dst approx astar values
                #store the poly coefs in m_apoly
                m_apoly = []
                for i in range(37):
                        dstDeg = 6
                        #Polynomial coefficients, highest power first
                        polycoefs = np.polyfit(dstTab,astarTab[8*(i+1):8*(i+2)],dstDeg)
                        m_apoly.append(polycoefs)

                #Find 3 referenced temp points around tr
                for i in range(37):
                        if tr<trTab[i]:
                                break
                i1 = max(i-1, 0)
                i2 = i1+3
                if (i2 > 36):
                        i2 = 36
                        i1 = i2 - 3
                #compute astar value for these 3 points
                values = []
                for j in range(i1,i2):
                        if (dst == 0.0):
                                values.append(astarTab[8*(j+1)])
                        else:
                                poly6 = np.poly1d(m_apoly[j])
                                values.append(poly6(dst))

                #interpolate to find real tr value
                trTab_log = []
                for j in range(len(trTab)):
                        trTab_log.append(np.log(trTab[j]))

                astar_interp = self.quadInterp(np.log(tr), trTab_log[i1:i2], values)
                return self.om22(tr,dst)/astar_interp

        def om11_CHEMKIN(self, tr, dst):
                # This is an overhaul of CANTERA version 2.3
                #range of dst
                dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

                #range of tr
                trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

                #tab of omega11 corresp. to (tr, dst)
                #CANTERA
                omegaTab = [4.008, 4.002, 4.655, 5.52, 6.454, 8.214, 9.824, 11.31,
                                3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160,
                                2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483 ,
                                2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539 ,
                                2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918 ,
                                1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747 ,
                                1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137 ,
                                1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872 ,
                                1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657 ,
                                1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780,
                                1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199,
                                1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99,
                                1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827,
                                1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7,
                                1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592,
                                1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397,
                                .95  ,  .9538,  .9656,  .9852, 1.012 , 1.082 , 1.168 , 1.265,
                                .9131,  .9162,  .9256,  .9413,  .9626, 1.019 , 1.09  , 1.17,
                                .8845,  .8871,  .8948,  .9076,  .9252,  .972 , 1.03  , 1.098,
                                .8428,  .8446,  .850 ,  .859 ,  .8716,  .9053,  .9483,  .9984,
                                .813 ,  .8142,  .8183,  .825 ,  .8344,  .8598,  .8927,  .9316,
                                .7898,  .791 ,  .794 ,  .7993,  .8066,  .8265,  .8526,  .8836,
                                .7711,  .772 ,  .7745,  .7788,  .7846,  .8007,  .822 ,  .8474,
                                .7555,  .7562,  .7584,  .7619,  .7667,  .78  ,  .7976,  .8189,
                                .7422,  .743 ,  .7446,  .7475,  .7515,  .7627,  .7776,  .796 ,
                                .72022, .7206,  .722 ,  .7241,  .7271,  .7354,  .7464,  .76  ,
                                .7025,  .703 ,  .704 ,  .7055,  .7078,  .7142,  .7228,  .7334,
                                .68776, .688,   .6888,  .6901,  .6919,  .697 ,  .704 ,  .7125,
                                .6751,  .6753,  .676 ,  .677 ,  .6785,  .6827,  .6884,  .6955,
                                .664 ,  .6642,  .6648,  .6657,  .6669,  .6704,  .6752,  .681,
                                .6414,  .6415,  .6418,  .6425,  .6433,  .6457,  .649 ,  .653,
                                .6235,  .6236,  .6239,  .6243,  .6249,  .6267,  .629 ,  .632,
                                .60882, .6089,  .6091,  .6094,  .61  ,  .6112,  .613 ,  .6154,
                                .5964,  .5964,  .5966,  .597 ,  .5972,  .5983,  .600 ,  .6017,
                                .5763,  .5763,  .5764,  .5766,  .5768,  .5775,  .5785,  .58,
                                .5415,  .5415,  .5416,  .5416,  .5418,  .542 ,  .5424,  .543,
                                .518 ,  .518 ,  .5182,  .5184,  .5184,  .5185,  .5186,  .5187]

                #First test on tr
                if (tr > 75.0):
                        omeg12 = 0.623 - 0.136e-2*tr + 0.346e-5*tr*tr - 0.343e-8*tr*tr*tr
                else:
                        #Find tr idx in trTab
                        if (tr <= 0.2):
                                ii = 1
                        else:
                                ii = 36
                        for i in range(1,37):
                                if (tr > trTab[i-1]) and (tr <= trTab[i]):
                                        ii = i
                                        break
                        #Find dst idx in dstTab
                        if (abs(dst) >= 1.0e-5):
                                if (dst <= 0.25):
                                        kk = 1
                                else:
                                        kk = 6
                                for i in range(1,7):
                                        if (dstTab[i-1] < dst) and (dstTab[i] >= dst):
                                                kk = i
                                                break
                                #Find surrounding values and interpolate
                                #First on dst
                                vert = np.zeros(3)
                                for i in range(3):
                                        arg = np.zeros(3)
                                        val = np.zeros(3)
                                        for k in range(3):
                                            arg[k] = dstTab[kk-1+k]
                                            val[k] = omegaTab[8*(ii-1+i) + (kk-1+k)]
                                        vert[i] = self.qinterp(dst, arg, val)
                                #Second on tr
                                arg = np.zeros(3)
                                for i in range(3):
                                    arg[i] = trTab[ii-1+i]
                                omeg12 = self.qinterp(tr, arg, vert)
                        else:
                                arg = np.zeros(3)
                                val = np.zeros(3)
                                for i in range(3):
                                    arg[i] = trTab[ii-1+i]
                                    val[i] = omegaTab[8*(ii-1+i)]
                                omeg12 =self. qinterp(tr, arg, val)
                return omeg12

        def om22_CHEMKIN(self, tr, dst):
                # This is an overhaul of CANTERA version 2.3
                #range of dst
                dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

                #range of tr
                trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

                #tab of omega22 corresp. to (tr, dst)
                #CANTERA
                omegaTab = [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89,
                                3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618,
                                2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874,
                                2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895,
                                2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249,
                                2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786,
                                1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435,
                                1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156,
                                1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933,
                                1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746,
                                1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451,
                                1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228,
                                1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053,
                                1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912,
                                1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795,
                                1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578,
                                1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428,
                                0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319,
                                0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236,
                                0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121,
                                0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044,
                                0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
                                0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
                                0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
                                0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
                                0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
                                0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
                                0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
                                0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
                                0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642,
                                0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339,
                                0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
                                0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
                                0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
                                0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
                                0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
                                0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885]

                #First test on tr
                if (tr > 75.0):
                        omeg12 = 0.703 - 0.146e-2*tr + 0.357e-5*tr*tr - 0.343e-8*tr*tr*tr
                else:
                        #Find tr idx in trTab
                        if (tr <= 0.2):
                                ii = 1
                        else:
                                ii = 36
                        for i in range(1,37):
                                if (tr > trTab[i-1]) and (tr <= trTab[i]):
                                        ii = i
                                        break
                        #Find dst idx in dstTab
                        if (abs(dst) >= 1.0e-5):
                                if (dst <= 0.25):
                                        kk = 1
                                else:
                                        kk = 6
                                for i in range(1,7):
                                        if (dstTab[i-1] < dst) and (dstTab[i] >= dst):
                                                kk = i
                                                break
                                #Find surrounding values and interpolate
                                #First on dst
                                vert = np.zeros(3)
                                for i in range(3):
                                        arg = np.zeros(3)
                                        val = np.zeros(3)
                                        for k in range(3):
                                            arg[k] = dstTab[kk-1+k]
                                            val[k] = omegaTab[8*(ii-1+i) + (kk-1+k)]
                                        vert[i] = self.qinterp(dst, arg, val)
                                #Second on tr
                                arg = np.zeros(3)
                                for i in range(3):
                                    arg[i] = trTab[ii-1+i]
                                omeg12 = self.qinterp(tr, arg, vert)
                        else:
                                arg = np.zeros(3)
                                val = np.zeros(3)
                                for i in range(3):
                                    arg[i] = trTab[ii-1+i]
                                    val[i] = omegaTab[8*(ii-1+i)]
                                omeg12 =self. qinterp(tr, arg, val)
                return omeg12

        def om22(self, tr, dst):
                # This is an overhaul of CANTERA version 2.3
                #range of dst
                dstTab = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5]

                #range of tr
                trTab = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                                5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
                                18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0]

                #tab of omega22 corresp. to (tr, dst)
                #CANTERA
                omegaTab = [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89,
                                3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618,
                                2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874,
                                2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895,
                                2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249,
                                2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786,
                                1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435,
                                1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156,
                                1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933,
                                1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746,
                                1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451,
                                1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228,
                                1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053,
                                1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912,
                                1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795,
                                1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578,
                                1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428,
                                0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319,
                                0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236,
                                0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121,
                                0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044,
                                0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
                                0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
                                0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
                                0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
                                0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
                                0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
                                0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
                                0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
                                0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642,
                                0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339,
                                0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
                                0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
                                0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
                                0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
                                0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
                                0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885]


                #Find for each fixed tr the poly of deg 6 in dst approx omega22 values
                #store the poly coefs in m_o22poly
                m_o22poly = []
                for i in range(37):
                        dstDeg = 6
                        #Polynomial coefficients, highest power first
                        polycoefs = np.polyfit(dstTab,omegaTab[8*i:8*(i+1)],dstDeg)
                        m_o22poly.append(polycoefs)

                #Find 3 referenced temp points around tr
                for i in range(37):
                        if tr<trTab[i]:
                                break
                i1 = max(i-1, 0)
                i2 = i1+3
                if (i2 > 36):
                        i2 = 36
                        i1 = i2 - 3
                #compute omega22 value for these 3 points
                values = []
                for j in range(i1,i2):
                        if (dst == 0.0):
                                values.append(omegaTab[8*j])
                        else:
                                poly6 = np.poly1d(m_o22poly[j])
                                values.append(poly6(dst))

                #interpolate to find real tr value
                trTab_log = []
                for j in range(len(trTab)):
                        trTab_log.append(np.log(trTab[j]))
                #print trTab_log[i1:i2], values
                om22_interp = self.quadInterp(np.log(tr), trTab_log[i1:i2], values)
                return om22_interp

        def qinterp(self, x0, x, y):
                val1 = y[0] + (x0-x[0])*(y[1]-y[0]) / (x[1]-x[0])
                val2 = y[1] + (x0-x[1])*(y[2]-y[1]) / (x[2]-x[1])
                fac1 = (x0-x[0]) / (x[1]-x[0]) / 2.0
                fac2 = (x[2]-x0) / (x[2]-x[1]) / 2.0
                if (x0 >= x[1]):
                    val = (val1*fac2+val2) / (1.0+fac2)
                else:
                    val = (val1+val2*fac1) / (1.0+fac1)
                return val

        def quadInterp(self, x0, x, y):
                dx21 = x[1] - x[0]
                dx32 = x[2] - x[1]
                dx31 = dx21 + dx32
                dy32 = y[2] - y[1]
                dy21 = y[1] - y[0]
                a = (dx21*dy32 - dy21*dx32)/(dx21*dx31*dx32)
                return a*(x0 - x[0])*(x0 - x[1]) + (dy21/dx21)*(x0 - x[1]) + y[1]


        def _getCVdRspecies(self, mechanism, t, species):
                models = mechanism.species(species.symbol).thermo
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

                if t < mid:
                        parameters = lowRange.parameters
                else:
                        parameters = highRange.parameters

                return ((parameters[0] - 1.0) + parameters[1] * t + parameters[2] * t * t \
                                + parameters[3] * t * t * t + parameters[4] * t * t * t * t)

