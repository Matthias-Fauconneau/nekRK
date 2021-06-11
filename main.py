#!/bin/env python
import sys
#from more_itertools import partition
def partition(pred, iterable):
    t1, t2 = tee(iterable)
    return filterfalse(pred, t1), filter(pred, t2)

from numpy import square as sq
from numpy import dot
cb = lambda x: x*x*x
from numpy import pi #π = pi
from numpy import sqrt
from numpy import log as ln
from numpy import log2
from numpy import polyfit as polynomial_regression

K = 1.380649e-23 #* J/kelvin
light_speed = 299792458.
mu_0 = 1.2566370621e-6 #  H/m (Henry=kg⋅m²/(s²A²))
epsilon_0 = 1./(light_speed*light_speed*mu_0) # F/m (Farad=s⁴A²/(m²kg)
NA = 6.02214076e23 #/mole
R = K*NA
Cm_per_Debye = 3.33564e-30 #C·m (Coulomb=A⋅s)
standard_atomic_weights = {'H': 1.008, 'C': 12.011, 'N': 14.0067, 'O': 15.999, 'Ar': 39.95}

class NASA7:
    def piece(self, T):
        return self.pieces[0 if T < self.temperature_split else 1]
    def molar_heat_capacity_at_constant_pressure_R(self, T):
        a = self.piece(T)
        return a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T

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

class Species:
    def interaction_well_depth(self, a, b):
        well_depth_J = self.well_depth_J
        return sqrt(well_depth_J[a]*well_depth_J[b]) * sq(self.xi(a, b))

    def T_star(self, a, b, T):
        return T * K / self.interaction_well_depth(a, b)

    def reduced_dipole_moment(self, a, b):
        (well_depth_J, permanent_dipole_moment, diameter) = (self.well_depth_J, self.permanent_dipole_moment, self.diameter)
        return permanent_dipole_moment[a]*permanent_dipole_moment[b] / (8. * pi * epsilon_0 * sqrt(well_depth_J[a]*well_depth_J[b]) * cb((diameter[a] + diameter[b])/2.))

    def collision_integral(self, table, a, b, T):
        ln_T_star = ln(self.T_star(a, b, T))
        header_ln_T_star = list(map(ln, header_T_star))
        interpolation_start_index = min((1+next(i for i,header_ln_T_star in enumerate(header_ln_T_star[1:]) if ln_T_star < header_ln_T_star))-1, len(header_ln_T_star)-3);
        header_ln_T_star_slice = header_ln_T_star[interpolation_start_index:][:3]
        assert(len(header_ln_T_star_slice) == 3)
        polynomials = table[interpolation_start_index:][:3]
        assert(len(polynomials) == 3)
        evaluate_polynomial = lambda P, x: dot(P, list(map(lambda k: pow(x, k), range(len(P)))))
        quadratic_interpolation = lambda x, y, x0: ((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
        delta_star = self.reduced_dipole_moment(a, b)
        return quadratic_interpolation(header_ln_T_star_slice, list(map(lambda P: evaluate_polynomial(P, delta_star), polynomials)), ln_T_star)
    #}

    def omega_star_22(self, a, b, T):
            return self.collision_integral(omega_star_22, a, b, T)

    def omega_star_11(self, a, b, T):
            return self.omega_star_22(a, b, T)/self.collision_integral(A_star, a, b, T)

    viscosity = lambda self, a, T: 5./16. * sqrt(pi * self.molar_mass[a]/NA * K*T) / (self.omega_star_22(a, a, T) * pi * sq(self.diameter[a]))

    def thermal_conductivity(self, a, T):
        (molar_mass, thermodynamics, diameter, well_depth_J, rotational_relaxation, internal_degrees_of_freedom) = (self.molar_mass, self.thermodynamics, self.diameter, self.well_depth_J, self.rotational_relaxation, self.internal_degrees_of_freedom)

        f_internal = molar_mass[a]/NA/(K * T) * self.binary_thermal_diffusion_coefficient(a,a,T) / self.viscosity(a, T);
        T_star= self.T_star(a, a, T)
        fz = lambda T_star: 1. + pow(pi, 3./2.) / sqrt(T_star) * (1./2. + 1./T_star) + (1./4. * sq(pi) + 2.) / T_star
        # Scaling factor for temperature dependence of rotational relaxation: Kee, Coltrin [2003:12.112, 2017:11.115]
        c1 = 2./pi * (5./2. - f_internal)/(rotational_relaxation[a] * fz(298.*K / well_depth_J[a]) / fz(T_star) + 2./pi * (5./3. * internal_degrees_of_freedom[a] + f_internal))
        f_translation = 5./2. * (1. - c1 * internal_degrees_of_freedom[a]/(3./2.))
        f_rotation = f_internal * (1. + c1)
        Cv_internal = thermodynamics[a].molar_heat_capacity_at_constant_pressure_R(T) - 5./2. - internal_degrees_of_freedom[a]
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
        return (self.diameter[a] + self.diameter[b])/2. * pow(self.xi(a, b), -1./6.)

    def binary_thermal_diffusion_coefficient(self, a, b, T):
        return 3./16. * sqrt(2.*pi/self.reduced_mass(a,b)) * pow(K*T, 3./2.) / (pi*sq(self.reduced_diameter(a,b))*self.omega_star_11(a, b, T))

    def transport_polynomials(self):
        N = 50
        (temperature_min, temperature_max) = (300., 3000.)
        T = list(map(lambda n: temperature_min + float(n) / float(N-1) * (temperature_max-temperature_min), range(N)))
        ln_T = list(map(ln, T))
        class TransportPolynomials:
            pass
        transport_polynomials = TransportPolynomials()
        transport_polynomials.sqrt_viscosity_T14 = map(lambda a: polynomial_regression(ln_T, list(map(lambda T: self.viscosity(a, T) / sqrt(sqrt(T)), T)), 3), range(self.len))
        sys.exit("\n".join(map(lambda p: f'{p}', transport_polynomials.sqrt_viscosity_T14)))
        transport_polynomials.thermal_conductivity_T12 = map(lambda a: polynomial_regression(ln_T, map(lambda T: self.thermal_conductivity(a, T) / sqrt(T), T), 3), range(self.len))
        print("transport_polynomials.binary_thermal_diffusion_coefficients_T32 using Python (slow version)")
        import time
        start_time = time.time()
        transport_polynomials.binary_thermal_diffusion_coefficients_T32 = map(
            lambda a: map(lambda b: polynomial_regression(map(ln, T), map(lambda T: self.binary_thermal_diffusion_coefficient(a, b, T) / pow(T, 3./2.), T), 3), range(len(self.species))), range(self.len))
        print("%ds"%(time.time() - start_time))
        return transport_polynomials

import ruamel.yaml
model = ruamel.yaml.YAML(typ='safe').load(open('/usr/share/cantera/data/LiDryer.yaml'))
species = model['species']
for inert in ['AR','N2']:
    def position(iterator, predicate):
        for index, item in enumerate(iterator):
                if predicate(item):
                    return index
        return None
    inert = position(species, lambda e: e['name']==inert)
    if inert:
        species.append(species.pop(inert)) # Move inert specie to last index
        break
else:
    print("Missing inert specie", file=sys.stderr)

self = Species()
self.len = len(species)
p = lambda f: list(map(f, species))
self.molar_mass = p(lambda s: sum(map(lambda element: element[1] * standard_atomic_weights[element[0]]/1e3, s['composition'].items())))
self.internal_degrees_of_freedom = p(lambda s: {'atom': 0, 'linear': 1, 'nonlinear': 3/2}[s['transport']['geometry']])
self.heat_capacity_ratio = p(lambda s: 1. + 2./{'atom': 3, 'linear': 5, 'nonlinear': 6}[s['transport']['geometry']])
self.well_depth_J = p(lambda s: s['transport']['well-depth']*K)
self.diameter = p(lambda s: s['transport']['diameter']*1e-10) #Å
self.permanent_dipole_moment = p(lambda s: s['transport'].setdefault('dipole',0)*Cm_per_Debye)
self.polarizability = p(lambda s: s['transport'].setdefault('polarizability',0)*1e-30) # Å³
self.rotational_relaxation = p(lambda s: float(s['transport'].setdefault('rotational-relaxation',0)))
def from_model(s):
    self = NASA7()
    self.temperature_split = s['temperature-ranges'][1]
    self.pieces = s['data']
    return self
self.thermodynamics = list(map(lambda s: from_model(s['thermo']), species))
self.species_names = list(map(lambda s: s['name'], species))

transport_polynomials = self.transport_polynomials() # {sqrt_viscosity_T14, thermal_conductivity_T12, binary_thermal_diffusion_coefficients_T32}

temperature_splits = {}
for index, specie in enumerate(thermodynamics):
    temperature_splits.setdefault(specie.temperature_split, []).append(index)

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

    reaction.reactants = map(lambda specie: sum(map(lambda _, coefficient: coefficient, filter(lambda s, _: s == specie, reaction.reactants))), species_names)
    reaction.products = map(lambda specie: sum(map(lambda _, coefficient: coefficient, filter(lambda s, _: s == specie, reaction.products))), species_names)
    reaction.net = map(lambda reactant, product: -reactant + product, zip(reaction.reactants, reaction.products))
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

select = lambda predicate, consequent, alternative: consequent if predicate else alternative
mul = lambda c, v: select(c==0, '', f"{select(c == 1, '', select(c == -1, '-', f'{c}'))}*{v}")
code = lambda f, l: '\n    '.join(map(f, l))
def product_of_exponentiations(c, v):
    for _c in c: assert(_c == int(_c))
    c = map(lambda c: int(c), c)
    (num, div) = partition(lambda _, c: c>0, filter(lambda _, c: c!=0, enumerate(c)))
    from itertools import repeat, chain
    num = '*'.join(chain(*map(lambda i, c: repeat(f'{v}[{i}',max(0,c)), num)))
    div = '*'.join(chain(*map(lambda i, c: repeat(f'{v}[{i}',max(0,-c)), div)))
    if (num=='') and (div==''): return '1.'
    elif div=='': return num
    elif num=='': return f'1./{div}'
    else: return f'{num}/{div}'
#}

piece = lambda expression, i: code(lambda i, nasa7: f'_[{i}] = {expression(nasa7.pieces[i])};', enumerate(thermodynamics))
thermodynamics = lambda e: code(lambda temperature_split, species: f'if (T < {temperature_split}) {{\n{piece(e, 0)}\n}}\n else\n {{\n{piece(e, 1)}}}', temperature_splits.items())

arrhenius = lambda rate_constant: f'exp2({-rate_constant.activation_temperature/ln(2)} * rcp_T + {rate_constant.temperature_exponent} * log_T + {log2(rate_constant.A)})'
rcp_arrhenius = lambda rate_constant: f'exp2({rate_constant.activation_temperature/ln(2)} * rcp_T - {rate_constant.temperature_exponent} * log_T - {log2(rate_constant.A)})'

def reaction(id, r):
    efficiency = f"({'+'.join(map(lambda specie, efficiency: f'concentrations[{specie}]' if efficiency == 1. else f'{efficiency}*concentrations[{specie}]', enumerate(r.efficiencies)))})"
    if r.type == "Elementary":
        c = f'c = {arrhenius(rate_constant)}'
    elif r.type == "ThreeBody":
        c = f'c = {arrhenius(rate_constant)} * {efficiency}'
    elif r.type == "PressureModification":
        c = f'''Pr = {arrhenius(r.k0)} * {efficiency};
        c = Pr / ({rcp_arrhenius(rate_constant)} * Pr + 1.)'''
    elif r.type == "Falloff":
        A, T3, T1, T2 = r.troe.A, r.troe.T3, r.troe.T1, r.troe.T2
        c = f'''Pr = {arrhenius(reaction.k0)} * {efficiency};
        logFcent = log2({1.-A} * exp2({1./(-ln(2)*T3)}*T) + {A} * exp2({1./(-ln(2)*T1)}*T) + exp2({-T2/ln(2)}*rcp_T));
        logPr_c = log2(Pr) - 0.67*logFcent - {0.4*log2(10)};
        f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-{0.75*log2(10.)});
        c = Pr / ({rcp_arrhenius(rate_constant)} * Pr + 1.) * exp2(logFcent/(f1*f1+1.))'''
    else: sys.exit(reaction.type)
    Rf = product_of_exponentiations(reactants, 'concentrations')
    if reaction.reversible:
        rcp_equilibrium_constant = product_of_exponentiations(net, 'exp_Gibbs0_RT');
        if -sum_net == 0: pass
        elif -sum_net == 1: rcp_equilibrium_constant += '* P0_RT'
        elif -sum_net == -1: rcp_equilibrium_constant += '* rcp_P0_RT'
        else: sys.exit(f'Σnet {sum_net}')
        Rr = f'{rcp_equilibrium_constant} * {product_of_exponentiations(products, "concentrations")}'
        R = f'{Rf} - {Rr}'
    else:
        R = Rf
    return f'''{c};
    const float cR{id} = c * {R};'''
#}

print(
f'''#define n_species {len(self)}
const float fg_molar_mass[{len(self.molar_mass)}] = {', '.join(map(lambda x: '%s'%x, self.molar_mass))};
const float fg_rcp_molar_mass[{len(self.molar_mass)}] = {', '.join(map(lambda x: '%s'%x, map(lambda w: 1./w, self.molar_mass)))};
void fg_molar_heat_capacity_at_constant_pressure_R(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {{
 {thermodynamics(lambda a: f'{a[0]} + {a[1]} * T + {a[2]} * T_2 + {a[3]} * T_3 + {a[4]} * T_4')}
}}
void fg_enthalpy_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {{
 {thermodynamics(lambda a: f'{a[0]} + {a[1]/2} * T + {a[2]/3} * T_2 + {a[3]/4} * T_3 + {a[4]/5} * T_4 + {a[5]} * rcp_T')}
}}
void fg_exp_Gibbs_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {{
 {thermodynamics(lambda a: f'exp2({a[5]/ln(2)} * rcp_T + {(a[0] - a[6])/ln(2)} + {-a[0]} * log_T + {-a[1]/2/ln(2)} * T + {(1./3.-1./2.)*a[2]/ln(2)} * T_2 + {(1./4.-1./3.)*a[3]/ln(2)} * T_3 + {(1./5.-1./4.)*a[4]/ln(2)} * T_4)')}
}}
float fg_viscosity(float T, const float mole_fractions[]) {{
    float T_14 =	sqrt(sqrt(T));
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T;
    float ln_T_3 = ln_T_2*ln_T;
    return pow({'+'.join(map(lambda i, P: f'+ mole_fractions[{i}] * pow(({P[0]} + {P[1]} * ln_T + {P[2]} * ln_T_2 + {P[3]} * ln_T_3)*T_14, {2*6})', enumerate(sqrt_vis_T14)))}), 1./6.);
}}

float fg_thermal_conductivity(float T, const float mole_fractions[]) {{
    float T_12 =	sqrt(T);
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T;
    float ln_T_3 = ln_T_2*ln_T;
    return pow({'+'.join(map(lambda i, P: f'T_12 * pow(({P[i][0]} + {P[i][1]}*ln_T + {P[i][2]}*ln_T_2 + {P[i][3]}*ln_T_3), 4) * mole_fractions[{i}]', enumerate(conductivity_T12)))}, 1./4.);
}}

float fg_mixture_diffusion_coefficients(float T, const float mole_fractions[]) {{
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T;
    float ln_T_3 = ln_T_2*ln_T;
    float T_32 = T*sqrt(T);
    {lines(lambda k: f"Ddiag[{k}] = (1. - mass_fractions[{k}]) * mole_fractions[{k}] / ({map(lambda j, P: f'mole_fractions[{i}] / (({P[i][0]} + {P[i][1]}*ln_T + {P[i][2]}*ln_T_2 + {P[i][3]}*ln_T_3)*T_32)', M[k])});", enumerate(M))}
}}

void fg_rates(const float log_T, const float T, const float T2, const float T4, const float rcp_T, const float rcp_T2, const float P0_RT, const float rcp_P0_RT, const float exp_Gibbs0_RT[], const float concentrations[], float* molar_rates) {{
    float c, Pr, logFcent, logPr_c, f1;
    {lines(reaction, enumerate(reactions))}
    {lines(lambda specie: f"molar_rates[{specie}] = {'+'.join(filter(None, map(lambda i, r: mul(r.net[specie],f'cR{i}'), enumerate(reactions))))};", range(species.len()-1))}
}}
''')
