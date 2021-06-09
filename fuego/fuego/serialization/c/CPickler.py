# coding=utf8
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2007 All Rights Reserved
from __future__ import print_function

def product_of_exponentiations(c, v):
    #from itertools import partition
    def partition(pred, iterable):
        evaluations = ((pred(x), x) for x in iterable)
        from itertools import tee
        t1, t2 = tee(evaluations)
        return (
                (x for (cond, x) in t1 if cond),
                (x for (cond, x) in t2 if not cond),
        )
    (num, div) = partition(lambda (_,c): c>0, filter(lambda (_,c): c!=0, enumerate(c)))
    from itertools import repeat, chain
    num = '*'.join(chain(*map(lambda (i,c): repeat("%s[%d]"%(v,i),max(0,c)), num)))
    div = '*'.join(chain(*map(lambda (i,c): repeat("%s[%d]"%(v,i),max(0,-c)), div)))
    if (num=='') and (div==''): return '1.'
    elif div=='': return num
    elif num=='': return '1./(%s)'%(div)
    else: return '%s/(%s)'%(num, div)
#}

import journal
from weaver.mills.CMill import CMill
from pyre.units.pressure import atm
from pyre.units.SI import meter, second, mole, kelvin
from pyre.units.length import cm
from pyre.units.energy import cal, kcal, J, kJ
K = 1.380649e-23 #* J/kelvin
light_speed = 299792458.
mu_0 = 1.2566370621e-6 #  H/m (Henry=kg⋅m²/(s²A²))
epsilon_0 = 1./(light_speed*light_speed*mu_0) # F/m (Farad=s⁴A²/(m²kg)

NA = 6.02214076e23 #/mole
#from pyre.handbook.constants.fundamental import avogadro as NA
R = K*NA
#from pyre.handbook.constants.fundamental import gas_constant as R
from numpy import square as sq
from numpy import dot
cb = lambda x: x*x*x
from numpy import pi #π = pi
from numpy import sqrt
from numpy import log as ln
from numpy import log2
from numpy import polyfit as polynomial_regression

class NASA7:
    def piece(self, T):
        return self.pieces[0 if T < self.temperature_split else 1]
    def molar_heat_capacity_at_constant_pressure_R(self, T):
        a = self.piece(T)
        return a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T

import sys
header_T_star = [sys.float_info.epsilon, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 50., 75., 100., 500.]
header_delta_star = [0., 1./4., 1./2., 3./4., 1., 3./2., 2., 5./2.]
omega_star_22 = [
[4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89], #FIXME
[4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89],
[3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618],
[2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874],
[2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895],
[2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249],
[2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786],
[1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435],
[1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156],
[1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933],
[1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746],
[1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451],
[1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228],
[1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053],
[1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912],
[1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795],
[1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578],
[1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428],
[0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319],
[0.96988, 0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236],
[0.92676, 0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121],
[0.89616, 0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044],
[0.87272, 0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893],
[0.85379, 0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482],
[0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916],
[0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901],
[0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504],
[0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212],
[0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983],
[0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797],
[0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642],
[0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339],
[0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112],
[0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932],
[0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784],
[0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546],
[0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147],
[0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885],
[0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885], #FIXME
];
A_star = [
[1.0065, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840, 1.0840],
[1.0231, 1.0660, 1.0380, 1.0400, 1.0430, 1.0500, 1.0520, 1.0510],
[1.0424, 1.0450, 1.0480, 1.0520, 1.0560, 1.0650, 1.0660, 1.0640],
[1.0719, 1.0670, 1.0600, 1.0550, 1.0580, 1.0680, 1.0710, 1.0710],
[1.0936, 1.0870, 1.0770, 1.0690, 1.0680, 1.0750, 1.0780, 1.0780],
[1.1053, 1.0980, 1.0880, 1.0800, 1.0780, 1.0820, 1.0840, 1.0840],
[1.1104, 1.1040, 1.0960, 1.0890, 1.0860, 1.0890, 1.0900, 1.0900],
[1.1114, 1.1070, 1.1000, 1.0950, 1.0930, 1.0950, 1.0960, 1.0950],
[1.1104, 1.1070, 1.1020, 1.0990, 1.0980, 1.1000, 1.1000, 1.0990],
[1.1086, 1.1060, 1.1020, 1.1010, 1.1010, 1.1050, 1.1050, 1.1040],
[1.1063, 1.1040, 1.1030, 1.1030, 1.1040, 1.1080, 1.1090, 1.1080],
[1.1020, 1.1020, 1.1030, 1.1050, 1.1070, 1.1120, 1.1150, 1.1150],
[1.0985, 1.0990, 1.1010, 1.1040, 1.1080, 1.1150, 1.1190, 1.1200],
[1.0960, 1.0960, 1.0990, 1.1030, 1.1080, 1.1160, 1.1210, 1.1240],
[1.0943, 1.0950, 1.0990, 1.1020, 1.1080, 1.1170, 1.1230, 1.1260],
[1.0934, 1.0940, 1.0970, 1.1020, 1.1070, 1.1160, 1.1230, 1.1280],
[1.0926, 1.0940, 1.0970, 1.0990, 1.1050, 1.1150, 1.1230, 1.1300],
[1.0934, 1.0950, 1.0970, 1.0990, 1.1040, 1.1130, 1.1220, 1.1290],
[1.0948, 1.0960, 1.0980, 1.1000, 1.1030, 1.1120, 1.1190, 1.1270],
[1.0965, 1.0970, 1.0990, 1.1010, 1.1040, 1.1100, 1.1180, 1.1260],
[1.0997, 1.1000, 1.1010, 1.1020, 1.1050, 1.1100, 1.1160, 1.1230],
[1.1025, 1.1030, 1.1040, 1.1050, 1.1060, 1.1100, 1.1150, 1.1210],
[1.1050, 1.1050, 1.1060, 1.1070, 1.1080, 1.1110, 1.1150, 1.1200],
[1.1072, 1.1070, 1.1080, 1.1080, 1.1090, 1.1120, 1.1150, 1.1190],
[1.1091, 1.1090, 1.1090, 1.1100, 1.1110, 1.1130, 1.1150, 1.1190],
[1.1107, 1.1110, 1.1110, 1.1110, 1.1120, 1.1140, 1.1160, 1.1190],
[1.1133, 1.1140, 1.1130, 1.1140, 1.1140, 1.1150, 1.1170, 1.1190],
[1.1154, 1.1150, 1.1160, 1.1160, 1.1160, 1.1170, 1.1180, 1.1200],
[1.1172, 1.1170, 1.1170, 1.1180, 1.1180, 1.1180, 1.1190, 1.1200],
[1.1186, 1.1190, 1.1190, 1.1190, 1.1190, 1.1190, 1.1200, 1.1210],
[1.1199, 1.1200, 1.1200, 1.1200, 1.1200, 1.1210, 1.1210, 1.1220],
[1.1223, 1.1220, 1.1220, 1.1220, 1.1220, 1.1230, 1.1230, 1.1240],
[1.1243, 1.1240, 1.1240, 1.1240, 1.1240, 1.1240, 1.1250, 1.1250],
[1.1259, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260],
[1.1273, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1280],
[1.1297, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1290],
[1.1339, 1.1340, 1.1340, 1.1350, 1.1350, 1.1340, 1.1340, 1.1320],
[1.1364, 1.1370, 1.1370, 1.1380, 1.1390, 1.1380, 1.1370, 1.1350],
[1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187, 1.14187],
]

def arrhenius(rate_constant):
    A, temperature_exponent, activation_temperature = rate_constant.preexponential_factor, rate_constant.temperature_exponent, rate_constant.activation_temperature
    return "exp2(%e * rcp_T + %e * log_T + %e)"%(-activation_temperature/ln(2), temperature_exponent, log2(A))

def rcp_arrhenius(rate_constant):
    A, temperature_exponent, activation_temperature = rate_constant.preexponential_factor, rate_constant.temperature_exponent, rate_constant.activation_temperature
    return "exp2(%e * rcp_T - %e * log_T - %e)"%(activation_temperature/ln(2), temperature_exponent, log2(A))

def rates(reactions):
    species_len = len(reactions[0].reactants)
    def reaction((reaction_index, reaction)):
        reactants, products, net, sum_net, rate_constant, efficiencies = \
            reaction.reactants, reaction.products, reaction.net, reaction.sum_net, reaction.rate_constant, reaction.efficiencies
        def dot((specie, efficiency)):
            if efficiency == 1.: return "concentrations[%d]" % (specie)
            else: return "%e*concentrations[%d]" % (efficiency, specie)
        if reaction.type == "Elementary":
            c = "const dfloat c = %s"%(arrhenius(rate_constant))
        elif reaction.type == "ThreeBody":
            c = "const dfloat c = %s * (%s)"%(arrhenius(rate_constant), " + ".join(map(dot, enumerate(efficiencies))))
        elif reaction.type == "PressureModification":
            c = "const dfloat Pr = %s * (%s);\n    "%(arrhenius(reaction.k0), " + ".join(map(dot, enumerate(efficiencies))))
            c += "const dfloat c = Pr / (%s * Pr + 1.)"%(rcp_arrhenius(rate_constant))
        elif reaction.type == "Falloff":
            c = "const dfloat Pr = %s * (%s);\n    "%(arrhenius(reaction.k0), " + ".join(map(dot, enumerate(efficiencies))))
            A, T3, T1, T2 = reaction.troe.A, reaction.troe.T3, reaction.troe.T1, reaction.troe.T2
            c += "const dfloat logFcent = log2(%e * exp2(%e*T) + %e * exp2(%e*T) + exp2(%e*rcp_T));\n    "%(1.-A, 1./(-ln(2)*T3), A, 1./(-ln(2)*T1), (-T2/ln(2)))
            c += "const dfloat logPr_c = log2(Pr) - 0.67*logFcent - %e;\n    "%(0.4*log2(10))
            c += "const dfloat f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-%e);\n    "%(0.75*log2(10.));
            c += "const dfloat c = Pr / (%s * Pr + 1.) * exp2(logFcent/(f1*f1+1.))"%(rcp_arrhenius(rate_constant))
        else: exit(reaction.type)

        Rf = product_of_exponentiations(reactants, 'concentrations')
        if reaction.reversible:
            rcp_equilibrium_constant = product_of_exponentiations(net, "exp_Gibbs0_RT");
            from sys import exit
            if -sum_net == 0: pass
            elif -sum_net == 1: rcp_equilibrium_constant += "* P0_RT"
            elif -sum_net == -1: rcp_equilibrium_constant += "* rcp_P0_RT"
            else: exit("Σnet %d"%sum_net)
            Rr = "%s * %s"%(rcp_equilibrium_constant, product_of_exponentiations(products, "concentrations"))
            R = "%s - %s"%(Rf, Rr)
        else:
            R = "%s /*irreversible*/"%(Rf)
        return "{%s;\n    cR[%d] = c * (%s);}"%(c, reaction_index, R)
    def specie(specie):
        def expr((specie, net)):
            if net == 1: return "cR[%d]"%(specie)
            elif net == -1: return "-cR[%d]"%(specie)
            else: return "%d*cR[%d]"%(net, specie)
        rate = '+'.join(map(expr, filter(lambda (_, net): net != 0, enumerate(map(lambda reaction: reaction.net[specie], reactions)))))
        return "molar_rates[%d] = %s;"%(specie, rate)
    return "void fg_rates(const dfloat log_T, const dfloat T, const dfloat T2, const dfloat T4, const dfloat rcp_T, const dfloat rcp_T2, const dfloat P0_RT, const dfloat rcp_P0_RT, const dfloat exp_Gibbs0_RT[], const dfloat concentrations[], dfloat* molar_rates) {\n    %s\n}\n"%(
    "\n    ".join(
        ["dfloat cR[%d];"%(len(reactions))] +
        #map(lambda reaction_index, reaction: expr0(reaction_index, reaction), enumerate(reactions)) +
        map(reaction, enumerate(reactions)) +
        map(specie, range(species_len-1))
    ))

class CPickler(CMill):
    def interaction_well_depth(self, a, b):
        well_depth_J = self.well_depth_J
        return sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.xi(a, b))

    def T_star(self, a, b, T):
        return T * K / self.interaction_well_depth(a, b)

    def reduced_dipole_moment(self, a, b):
        (well_depth_J, permanent_dipole_moment, diameter) = (self.well_depth_J, self.permanent_dipole_moment, self.diameter)
        return permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * pi * epsilon_0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))

    def collision_integral(self, table, a, b, T):
        print(header_T_star, a, b, T, self.T_star(a, b, T))
        ln_T_star = ln(self.T_star(a, b, T))
        delta_star = self.reduced_dipole_moment(a, b)
        header_ln_T_star = map(ln, header_T_star)
        interpolation_start_index = min((1+next(i for i,header_ln_T_star in enumerate(header_ln_T_star[slice(1, -1)]) if ln_T_star < header_ln_T_star))-1, len(header_ln_T_star)-3);
        header_ln_T_star_slice = header_ln_T_star[slice(interpolation_start_index, -1)][slice(3)]
        polynomials = table[slice(interpolation_start_index, -1)][slice(3)]
        evaluate_polynomial = lambda P, x: dot(P, map(lambda k: pow(x, k), range(len(P))))
        quadratic_interpolation = lambda x, y, x0: ((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
        return quadratic_interpolation(header_ln_T_star_slice, map(lambda P: evaluate_polynomial(P, delta_star), polynomials), ln_T_star)
    #}

    def omega_star_22(self, a, b, T):
            return self.collision_integral(omega_star_22, a, b, T)

    def omega_star_11(self, a, b, T):
            return self.omega_star_22(a, b, T)/self.collision_integral(A_star, a, b, T)

    def transport_polynomials(self):
        N = 50
        (temperature_min, temperature_max) = (300., 3000.)
        T = map(lambda n: temperature_min + float(n) / float(N-1) * (temperature_max-temperature_min), range(N))
        class TransportPolynomials:
            pass
        transport_polynomials = TransportPolynomials()
        transport_polynomials.sqrt_viscosity_T14 = map(lambda a: polynomial_regression(map(ln, T), map(lambda T: self.viscosity(a, T) / sqrt(sqrt(T)), T), 3), range(len(self.species)))
        transport_polynomials.thermal_conductivity_T12 = map(lambda a: polynomial_regression(map(ln, T), map(lambda T: self.thermal_conductivity(a, T) / sqrt(T), T), 3), range(len(self.species)))
        transport_polynomials.binary_thermal_diffusion_coefficients_T32 = map(
            lambda a: map(lambda b: polynomial_regression(map(ln, T), map(lambda T: self.binary_thermal_diffusion_coefficient(a, b, T) / pow(T, 3./2.), T), 3), range(len(self.species))), range(len(self.species)))
        return transport_polynomials

    def thermodynamic_function(self, name, expression):
        temperature_splits = {}
        for index, specie in enumerate(self.thermodynamics):
            temperature_splits.setdefault(specie.temperature_split, []).append(index)

        self._write('void %s(const dfloat log_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {' % name)
        self._indent()
        for temperature_split, species in temperature_splits.items():
            self._write('if (T < %s) {' % temperature_split)
            self._indent()
            for specie in species:
                    self._write('species[%d] = %s;' % (specie, expression(self.thermodynamics[specie].pieces[0])))
            self._outdent()
            self._write('} else {')
            self._indent()
            for specie in species:
                    self._write('species[%d] = %s;' % (specie, expression(self.thermodynamics[specie].pieces[1])))
            self._outdent()
            self._write('}')
        self._outdent()
        self._write('}\n')

    def viscosity_function(self, sqrt_viscosity_T14):
        self._write('dfloat fg_viscosity(dfloat T, const dfloat mole_fractions[]) {')
        self._indent()
        self._write('dfloat T_14 =	sqrt(sqrt(T));')
        self._write('dfloat ln_T = log(T);')
        self._write('dfloat ln_T_2 = ln_T*ln_T; ')
        self._write('dfloat ln_T_3 = ln_T_2*ln_T; ')
        self._write('return pow(0.')
        for i in range(len(self.species)):
            P = sqrt_viscosity_T14[i]
            self._write('+ mole_fractions[%d] * pow((%.8E + %.8E * ln_T + %.8E * ln_T_2 + %.8E * ln_T_3)*T_14, 12.)'%(i, P[0], P[1], P[2], P[3])) #2.*6.
        self._write(', 1./6.);')
        self._outdent()
        self._write('}')

    def thermal_conductivity_function(self, thermal_conductivity_T12):
        self._write('dfloat fg_thermal_conductivity(dfloat T, const dfloat mole_fractions[]) {')
        self._indent()
        self._write('dfloat T_12 = sqrt(T);')
        self._write('dfloat ln_T = log(T);')
        self._write('dfloat ln_T_2 = ln_T*ln_T; ')
        self._write('dfloat ln_T_3 = ln_T_2*ln_T; ')
        self._write('return pow(0.')
        for i in range(len(self.species)):
            P = thermal_conductivity_T12[i]
            self._write('+ mole_fractions[%d] * pow((%.8E + %.8E * ln_T + %.8E * ln_T_2 + %.8E * ln_T_3)*T_12, 4.)'%(i, P[0], P[1], P[2], P[3]))
        self._write(', 1./4.);')
        self._outdent()
        self._write('}')

    def mixture_diffusion_coefficients_function(self, binary_thermal_diffusion_coefficients_T32):
        self._write('void fg_mixture_diffusion_coefficients(const dfloat mole_fractions[n_species], const dfloat mass_fractions[n_species], dfloat T, dfloat* Ddiag) {')
        self._indent()
        self._write('dfloat T_12 = sqrt(T);')
        self._write('dfloat ln_T = log(T);')
        self._write('dfloat ln_T_2 = ln_T*ln_T; ')
        self._write('dfloat ln_T_3 = ln_T_2*ln_T; ')
        self._write('dfloat T_32 = T*T_12;')
        for k in range(len(self.species)):
            self._write('Ddiag[%d] = (1. - mass_fractions[%d]) * mole_fractions[%d] / ( 0.' % (k,k,k))
            for j in range(len(self.species)):
                if k != j:
                    a = max(k, j)
                    b = min(k, j)
                    P = binary_thermal_diffusion_coefficients_T32[a][b]
                    self._write('+ mole_fractions[%d] / ((%.8E + %.8E * ln_T + %.8E * ln_T_2 + %.8E * ln_T_3)*T_32)' % (j, P[0], P[1], P[2], P[3]))
            self._write(');')
        self._outdent()
        self._write('}')

    def _renderDocument(self, mechanism, options=None):
        species = mechanism.species()
        self.species = species
        import pyre.handbook
        molar_mass = map(lambda s: sum(map(lambda (element, count): count * pyre.handbook.periodicTable().symbol(element.capitalize()).atomicWeight/1e3, s.composition)), species)
        f = open("species", "w")
        to_string = lambda x: '%s'%x
        f.write("molar_mass: "+ to_string(molar_mass) +"\n")
        self.internal_degrees_of_freedom = map(lambda s: [0., 1., 3./2.][s.trans[0].parameters[0]], species) #transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. }
        f.write("internal_degrees_of_freedom: "+ to_string(self.internal_degrees_of_freedom) +"\n")
        self.heat_capacity_ratio = map(lambda s: 1. + 2./[3., 5., 6.][s.trans[0].parameters[0]], species) #1. + 2. / match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. };
        f.write("heat_capacity_ratio: "+ to_string(self.heat_capacity_ratio) +"\n")
        self.well_depth_J = map(lambda s: float(s.trans[0].eps)*K, species) #well_depth_K * K
        f.write("well_depth_J: "+ to_string(self.well_depth_J) +"\n")
        self.diameter = map(lambda s: float(s.trans[0].sig)*1e-10, species) #diameter_Å*1e-10
        f.write("diameter: "+ to_string(self.diameter) +"\n")
        Cm_per_Debye = 3.33564e-30 #C·m (Coulomb=A⋅s)
        self.permanent_dipole_moment = map(lambda s: float(s.trans[0].dip)*Cm_per_Debye, species)
        f.write("permanent_dipole_moment: "+ to_string(self.permanent_dipole_moment) +"\n")
        self.polarizability = map(lambda s: float(s.trans[0].pol)*1e-30, species) # polarizability_Å3*1e-30
        f.write("polarizability: "+ to_string(self.polarizability) +"\n")
        self.rotational_relaxation = map(lambda s: float(s.trans[0].zrot), species)
        f.write("rotational_relaxation: "+ to_string(self.rotational_relaxation) +"\n")
        def block_expr(s):
            self = NASA7()
            [high, low] = s.thermo # /!\ Reverse order
            assert(low.highT == high.lowT)
            self.temperature_split = low.highT
            self.pieces = [low.parameters, high.parameters]
            return self
        self.thermodynamics = map(block_expr, species)
        f.write("thermodynamics: "+ to_string(self.thermodynamics) +"\n")
        #f.write(to_string(self))
        f.close()
        species_names = map(lambda s: s.symbol, species)
        self.species_names = species_names

        def from_fuego(species_names, reaction):
            assert(reaction.units["prefactor"]=='mole/cm**3')
            assert(reaction.units["activation"]=='cal/mole')
            preexponential_factor, temperature_exponent, activation_energy_cal = reaction.arrhenius

            reaction.type = "Elementary"
            if reaction.thirdBody:
                reaction.type = "ThreeBody"
            if reaction.low:
                reaction.type = "PressureModification"
            if reaction.troe:
                reaction.type = "Falloff"

            reaction.reactants = map(lambda specie: sum(map(lambda (_, coefficient): coefficient, filter(lambda (s, _): s == specie, reaction.reactants))), species_names)
            reaction.products = map(lambda specie: sum(map(lambda (_, coefficient): coefficient, filter(lambda (s, _): s == specie, reaction.products))), species_names)
            reaction.net = map(lambda (reactant, product): -reactant + product, zip(reaction.reactants, reaction.products))
            reaction.sum_net = sum(reaction.net)

            class RateConstant(): pass
            reaction.rate_constant = RateConstant()
            concentration_cm3_unit_conversion_factor_exponent = sum(reaction.reactants)
            if reaction.type == "ThreeBody":
                concentration_cm3_unit_conversion_factor_exponent += 1
            reaction.rate_constant.preexponential_factor = preexponential_factor * pow(1e-6, concentration_cm3_unit_conversion_factor_exponent-1)
            reaction.rate_constant.temperature_exponent = temperature_exponent
            J_per_cal = 4.184
            reaction.rate_constant.activation_temperature = activation_energy_cal * J_per_cal / (K*NA)

            if reaction.thirdBody:
                specie, coefficient = reaction.thirdBody
                if not reaction.efficiencies:
                        assert(specie != "<mixture>")
                        reaction.efficiencies = {specie: coefficient}
            else:
                    assert(not reaction.efficiencies)
            reaction.efficiencies = map(lambda specie: dict(reaction.efficiencies).get(specie, 1.), species_names)
            if reaction.low:
                preexponential_factor, temperature_exponent, activation_energy_cal = reaction.low
                reaction.k0 = RateConstant()
                reaction.k0.preexponential_factor = preexponential_factor * pow(1e-6, sum(reaction.reactants)-1+1),
                reaction.k0.temperature_exponent = temperature_exponent
                J_per_cal = 4.184
                reaction.k0.activation_temperature = activation_energy_cal * J_per_cal / (K*NA)
            if reaction.troe:
                A, T3, T1 = reaction.troe[0], reaction.troe[1], reaction.troe[2]
                if len(reaction.troe) == 4: T2 = reaction.troe[3]
                else: T2 = 0
                class Troe(): pass
                reaction.troe = Troe()
                reaction.troe.A = A
                reaction.troe.T3 = T3
                reaction.troe.T1 = T1
                reaction.troe.T2 = T2
            #}
            return reaction

        self.reactions = map(lambda reaction: from_fuego(species_names, reaction), mechanism.reaction())
        # Species
        self._write('#define n_species %d' % (len(species)))
        self._write('const dfloat fg_molar_mass[%s] = {%s};'%(len(molar_mass), ', '.join(map(to_string, molar_mass))))
        self._write('const dfloat fg_rcp_molar_mass[%s] = {%s};'%(len(molar_mass), ', '.join(map(to_string, map(lambda w: 1./w, molar_mass)))))

        def molar_heat_capacity_at_constant_pressure_R_expression(a):
            return '%+15.8e %+15.8e * T %+15.8e * T_2 %+15.8e * T_3 %+15.8e * T_4' % (a[0], a[1], a[2], a[3], a[4])
        self.thermodynamic_function("fg_molar_heat_capacity_at_constant_pressure_R", molar_heat_capacity_at_constant_pressure_R_expression)
        def enthalpy_RT_expression(a):
            return '%+15.8e %+15.8e * T %+15.8e * T_2 %+15.8e * T_3 %+15.8e * T_4 %+15.8e * rcp_T' % (a[0], a[1]/2, a[2]/3, a[3]/4, a[4]/5, a[5])
        self.thermodynamic_function("fg_enthalpy_RT", enthalpy_RT_expression)
        def exp_Gibbs_RT_expression(a):
            return 'exp2(%+15.8e * rcp_T %+15.8e %+15.8e * log_T %+15.8e * T %+15.8e * T_2 %+15.8e * T_3 %+15.8e * T_4)' % \
                              (a[5]/ln(2),             (a[0] - a[6])/ln(2), -a[0], -a[1]/2/ln(2), (1./3.-1./2.)*a[2]/ln(2), (1./4.-1./3.)*a[3]/ln(2), (1./5.-1./4.)*a[4]/ln(2))
        self.thermodynamic_function("fg_exp_Gibbs_RT", exp_Gibbs_RT_expression)

        self.molar_mass = molar_mass
        #transport_polynomials = self.transport_polynomials() # {sqrt_viscosity_T14, thermal_conductivity_T12, binary_thermal_diffusion_coefficients_T32}
        #f.write(to_string(transport_polynomials.binary_thermal_diffusion_coefficients_T32))
        #self.viscosity_function(transport_polynomials.sqrt_viscosity_T14)
        #self.thermal_conductivity_function(transport_polynomials.thermal_conductivity_T12)
        #self.mixture_diffusion_coefficients_function(transport_polynomials.binary_thermal_diffusion_coefficients_T32)

        # Reactions
        self._write(rates(self.reactions))

    def viscosity(self, a, T):
        return 5./16. * sqrt(pi * self.molar_mass[a]/NA * K*T) / (self.omega_star_22(a, a, T) * pi * sq(self.diameter[a]))

    def thermal_conductivity(self, a, T):
        (molar_mass, thermodynamics, diameter, well_depth_J, rotational_relaxation, internal_degrees_of_freedom) = (self.molar_mass, self.thermodynamics, self.diameter, self.well_depth_J, self.rotational_relaxation, self.internal_degrees_of_freedom)

        f_internal = molar_mass[a]/NA/(K * T) * self.binary_thermal_diffusion_coefficient(a,a,T) / self.viscosity(a, T);
        T_star= self.T_star(a, a, T)
        fz = lambda T_star: 1. + pow(pi, 3./2.) / sqrt(T_star) * (1./2. + 1./T_star) + (1./4. * sq(pi) + 2.) / T_star
        # Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
        c1 = 2./pi * (5./2. - f_internal)/(rotational_relaxation[a] * fz(298.*K / well_depth_J[a]) / fz(T_star) + 2./pi * (5./3. * internal_degrees_of_freedom[a] + f_internal))
        f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.))
        f_rotation = f_internal * (1. + c1)
        Cv_internal = molar_heat_capacity_at_constant_pressure_R(thermodynamics[a], T) - 5./2. - internal_degrees_of_freedom[a]
        return (self.viscosity(a, T)/(molar_mass[a]/NA))*K*(f_translation * 3./2. + f_rotation * internal_degrees_of_freedom[a] + f_internal * Cv_internal)

    def reduced_mass(self, a, b):
        molar_mass = self.molar_mass
        return molar_mass[a]/NA * molar_mass[b]/NA / (molar_mass[a]/NA + molar_mass[b]/NA)

    def xi(self, a, b):
        (permanent_dipole_moment, polarizability, diameter, well_depth_J) = (self.permanent_dipole_moment, self.polarizability, self.diameter, self.well_depth_J)

        def expr():
            (polar, non_polar) = (a,b) if permanent_dipole_moment[a] != 0. else (b,a)
            return 1. + 1./4. * polarizability[non_polar]/cb(diameter[non_polar]) * permanent_dipole_moment[polar]/sqrt(well_depth_J[polar]*cb(diameter[polar])) * sqrt(well_depth_J[polar]/well_depth_J[non_polar])
        return 1. if (permanent_dipole_moment[a]>0.) == (permanent_dipole_moment[b]>0.) else expr()
    #}

    def reduced_diameter(self, a, b):
        diameter = self.diameter

        return (diameter[a] + diameter[b])/2. * pow(self.xi(a, b), -1./6.)

    def binary_thermal_diffusion_coefficient(self, a, b, T):
        return 3./16. * sqrt(2.*pi/self.reduced_mass(a,b)) * pow(K*T, 3./2.) / (pi*sq(self.reduced_diameter(a,b))*self.omega_star_11(a, b, T))
