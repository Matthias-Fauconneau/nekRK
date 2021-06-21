#!/bin/env python
from sys import float_info, argv, stderr
#from more_itertools import partition
def partition(pred, iterable):
    from itertools import tee, filterfalse
    t1, t2 = tee(iterable)
    return list(filterfalse(pred, t1)), list(filter(pred, t2))

from numpy import square as sq
cb = lambda x: x*x*x
from numpy import dot
from numpy import pi #π = pi
from numpy import sqrt
from numpy import log2
from numpy import log as ln
from numpy import polynomial
polynomial_regression = lambda X, Y, degree=3: polynomial.polynomial.polyfit(X, Y, deg=degree, w=[1/sq(y) for y in Y])
from numpy import linspace

K = 1.380649e-23 #* J/kelvin
light_speed = 299792458.
mu_0 = 1.2566370621e-6 #  H/m (Henry=kg⋅m²/(s²A²))
epsilon_0 = 1./(light_speed*light_speed*mu_0) # F/m (Farad=s⁴A²/(m²kg)
NA = 6.02214076e23 #/mole
R = K*NA
Cm_per_Debye = 3.33564e-30 #C·m (Coulomb=A⋅s)
J_per_cal = 4.184
standard_atomic_weights = {'H': 1.008, 'C': 12.011, 'N': 14.0067, 'O': 15.999, 'Ar': 39.95}

class NASA7:
    def piece(self, T):
        return self.pieces[0 if T < self.temperature_split else 1]
    def molar_heat_capacity_at_constant_pressure_R(self, T):
        a = self.piece(T)
        return a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T

header_T_star = [float_info.epsilon, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 50., 75., 100., 500.]
header_delta_star = [0., 1./4., 1./2., 3./4., 1., 3./2., 2., 5./2.]
collision_integrals_Omega_star_22 = [
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
collision_integrals_A_star = [
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

# Least square fits polynomials in δ⃰, for each T⃰  row of the collision integrals tables
polynomial_regression_delta_star = lambda table: [polynomial_regression(header_delta_star, T_star_row, degree=6) for T_star_row in table]
Omega_star_22 = polynomial_regression_delta_star(collision_integrals_Omega_star_22);
A_star = polynomial_regression_delta_star(collision_integrals_A_star);

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
        evaluate_polynomial = lambda P, x: dot(P, [pow(x, k) for k in range(len(P))])
        quadratic_interpolation = lambda x, y, x0: ((x[1]-x[0])*(y[2]-y[1])-(y[1]-y[0])*(x[2]-x[1]))/((x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1]))*(x0 - x[0])*(x0 - x[1]) + ((y[1]-y[0])/(x[1]-x[0]))*(x0-x[1]) + y[1]
        delta_star = self.reduced_dipole_moment(a, b)
        for P in polynomials: assert len(P) == 7, len(P)
        return quadratic_interpolation(header_ln_T_star_slice, [evaluate_polynomial(P, delta_star) for P in polynomials], ln_T_star)
    #}

    def Omega_star_22(self, a, b, T):
            return self.collision_integral(Omega_star_22, a, b, T)

    def Omega_star_11(self, a, b, T):
            return self.Omega_star_22(a, b, T)/self.collision_integral(A_star, a, b, T)

    viscosity = lambda self, a, T: 5./16. * sqrt(pi * self.molar_mass[a]/NA * K*T) / (self.Omega_star_22(a, a, T) * pi * sq(self.diameter[a]))

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
        return 3./16. * sqrt(2.*pi/self.reduced_mass(a,b)) * pow(K*T, 3./2.) / (pi*sq(self.reduced_diameter(a,b))*self.Omega_star_11(a, b, T))

    def transport_polynomials(self):
        T = linspace(300., 3000., 50)
        class TransportPolynomials:
            pass
        _ = TransportPolynomials()
        _.sqrt_viscosity_T14 = [polynomial_regression(ln(T), [sqrt(self.viscosity(a, T)) / sqrt(sqrt(T)) for T in T]) for a in range(self.len)]
        _.thermal_conductivity_T12 = [polynomial_regression(ln(T), [self.thermal_conductivity(a, T) / sqrt(T) for T in T]) for a in range(self.len)]
        #print("_.binary_thermal_diffusion_coefficients_T32 using Python (slow version)", file=stderr)
        import time
        start_time = time.time()
        _.binary_thermal_diffusion_coefficients_T32 = [[polynomial_regression(ln(T), [self.binary_thermal_diffusion_coefficient(a, b, T) / (T*sqrt(T)) for T in T]) for b in range(self.len)] for a in range(self.len)]
        #print(f'{int(time.time() - start_time)}s', file=stderr)
        return _

import ruamel.yaml
model = ruamel.yaml.YAML(typ='safe').load(open(argv[1]))
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
    print("Missing inert specie", file=stderr)

def from_model(species):
    self = Species()
    self.len = len(species)
    p = lambda f: list(map(f, species))
    self.names = p(lambda s: s['name'])
    self.molar_mass = p(lambda s: sum([element[1] * standard_atomic_weights[element[0]]/1e3 for element in s['composition'].items()]))
    def from_model(s): self = NASA7(); self.temperature_split = s['temperature-ranges'][1]; self.pieces = s['data']; return self
    self.thermodynamics = p(lambda s: from_model(s['thermo']))
    self.internal_degrees_of_freedom = p(lambda s: {'atom': 0, 'linear': 1, 'nonlinear': 3/2}[s['transport']['geometry']])
    self.heat_capacity_ratio = p(lambda s: 1. + 2./{'atom': 3, 'linear': 5, 'nonlinear': 6}[s['transport']['geometry']])
    self.well_depth_J = p(lambda s: s['transport']['well-depth']*K)
    self.diameter = p(lambda s: s['transport']['diameter']*1e-10) #Å
    self.permanent_dipole_moment = p(lambda s: s['transport'].get('dipole',0)*Cm_per_Debye)
    self.polarizability = p(lambda s: s['transport'].get('polarizability',0)*1e-30) # Å³
    self.rotational_relaxation = p(lambda s: float(s['transport'].get('rotational-relaxation',0)))
    return self
species = from_model(species)

def from_model(species_names, r):
    class Reaction(): pass
    reaction = Reaction()
    reaction.type = r.get('type', 'elementary')
    reaction.reversible = r.get('reversible', True)
    import re
    [reaction.reactants, reaction.products] = [[sum([c for (_, c) in filter(lambda s: s[0] == specie, side)]) for specie in species.names] for side in
                            [[(s.split(' ')[1], int(s.split(' ')[0])) if ' ' in s else (s, 1) for s in [s.strip() for s in side.removesuffix('+ M').removesuffix('(+M)').split(' + ')]] for side in [s.strip() for s in re.split('<?=>', r['equation'])]]]
    reaction.net = [-reactant + product for reactant, product in zip(reaction.reactants, reaction.products)]
    reaction.sum_net = sum(reaction.net)

    class RateConstant(): pass
    reaction.rate_constant = RateConstant()
    concentration_cm3_unit_conversion_factor_exponent = sum(reaction.reactants)
    if reaction.type == "three-body":
        concentration_cm3_unit_conversion_factor_exponent += 1
    rate_constant = r.get('rate-constant', r.get('high-P-rate-constant'))
    reaction.rate_constant.preexponential_factor = rate_constant['A'] * pow(1e-6, concentration_cm3_unit_conversion_factor_exponent-1)
    reaction.rate_constant.temperature_exponent = rate_constant['b']
    reaction.rate_constant.activation_temperature = rate_constant['Ea'] * J_per_cal / (K*NA)
    if r.get('efficiencies'): reaction.efficiencies = [r['efficiencies'].get(specie, r.get('default-efficiency', 1)) for specie in species.names]
    k0 = r.get('low-P-rate-constant')
    if k0:
        reaction.k0 = RateConstant()
        reaction.k0.preexponential_factor = k0['A'] * pow(1e-6, sum(reaction.reactants)-1+1)
        reaction.k0.temperature_exponent = k0['b']
        reaction.k0.activation_temperature = k0['Ea'] * J_per_cal / (K*NA)
    if reaction.type == 'falloff':
        class Troe(): pass
        reaction.troe = Troe()
        reaction.troe.A = r['Troe']['A']
        reaction.troe.T3 = r['Troe']['T3']
        reaction.troe.T1 = r['Troe']['T1']
        reaction.troe.T2 = r['Troe'].get('T2', 0)
    #}
    return reaction
reactions = [from_model(species.names, reaction) for reaction in model['reactions']]

temperature_splits = {}
for index, specie in enumerate(species.thermodynamics):
    temperature_splits.setdefault(specie.temperature_split, []).append(index)

transport_polynomials = species.transport_polynomials() # {sqrt_viscosity_T14, thermal_conductivity_T12, binary_thermal_diffusion_coefficients_T32}

code = lambda lines: '\n\t'.join(lines)
mul = lambda c, v: None if c==0 else f"{'' if c == 1 else '-' if c == -1 else f'{c}*'}{v}"
def product_of_exponentiations(c, v):
    for _c in c: assert(_c == int(_c))
    c = [int(c) for c in c]
    (div, num) = partition(lambda c: c[1]>0, filter(lambda c: c[1]!=0, enumerate(c)))
    from itertools import repeat, chain
    flatten = lambda t: [item for sublist in t for item in sublist]
    num = '*'.join(flatten([repeat(f'{v}[{i}]',max(0,c)) for i, c in num]))
    div   = '*'.join(flatten([repeat(f'{v}[{i}]',max(0,-c)) for i ,c in div]))
    if (num=='') and (div==''): return '1.'
    elif div=='': return num
    elif num=='': return f'1./({div})'
    else: return f'{num}/({div})'
#}

piece = lambda specie_indices, expression, piece: code([f'_[{specie}] = {expression(species.thermodynamics[specie].pieces[piece])};' for specie in specie_indices])
thermodynamics = lambda e: code([f'if (T < {temperature_split}) {{\n\t{piece(species, e, 0)}\n }} else {{\n\t{piece(species, e, 1)}\n }}' for temperature_split, species in temperature_splits.items()])

arrhenius = lambda r: f'exp2({-r.activation_temperature/ln(2)} * rcp_T + {r.temperature_exponent} * log_T + {log2(r.preexponential_factor)})'
rcp_arrhenius = lambda r: f'exp2({r.activation_temperature/ln(2)} * rcp_T - {r.temperature_exponent} * log_T - {log2(r.preexponential_factor)})'

def reaction(id, r):
    if hasattr(r, 'efficiencies'):
        efficiency = f"({'+'.join([f'concentrations[{specie}]' if efficiency == 1. else f'{efficiency}*concentrations[{specie}]' for specie, efficiency in enumerate(r.efficiencies)])})"
    if r.type == "elementary":
        c = f'c = {arrhenius(r.rate_constant)}'
    elif r.type == "three-body":
        c = f'c = {arrhenius(r.rate_constant)} * {efficiency}'
    elif r.type == "pressure-modification":
        c = f'''Pr = {arrhenius(r.k0)} * {efficiency};
        c = Pr / ({rcp_arrhenius(r.rate_constant)} * Pr + 1.)'''
    elif r.type == "falloff":
        A, T3, T1, T2 = r.troe.A, r.troe.T3, r.troe.T1, r.troe.T2
        c = f'''
        Pr = {arrhenius(r.k0)} * {efficiency};
        logFcent = log2({1.-A} * exp2({1./(-ln(2)*T3)}*T) + {A} * exp2({1./(-ln(2)*T1)}*T) + exp2({-T2/ln(2)}*rcp_T));
        logPr_c = log2(Pr) - 0.67*logFcent - {0.4*log2(10)};
        f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-{0.75*log2(10.)});
        c = Pr / ({rcp_arrhenius(r.rate_constant)} * Pr + 1.) * exp2(logFcent/(f1*f1+1.))'''
    else: exit(r.type)
    Rf = product_of_exponentiations(r.reactants, 'concentrations')
    if r.reversible:
        rcp_equilibrium_constant = product_of_exponentiations(r.net, 'exp_Gibbs0_RT');
        if -r.sum_net == 0: pass
        elif -r.sum_net == 1: rcp_equilibrium_constant += '* P0_RT'
        elif -r.sum_net == -1: rcp_equilibrium_constant += '* rcp_P0_RT'
        else: raise(f'Σnet {sum_net}')
        Rr = f'{rcp_equilibrium_constant} * {product_of_exponentiations(r.products, "concentrations")}'
        R = f'({Rf} - {Rr})'
    else:
        R = Rf
    return f'''{c};
    const float cR{id} = c * {R};'''
#}

line= '\n\t'
print(
f"""float sq(float x) {{ return x*x; }}
#define n_species {species.len}
const float fg_molar_mass[{species.len}] = {{{', '.join([f'{w}' for w in species.molar_mass])}}};
const float fg_rcp_molar_mass[{species.len}] = {{{', '.join([f'{1./w}' for w in species.molar_mass])}}};
void fg_molar_heat_capacity_at_constant_pressure_R(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {{
 {thermodynamics(lambda a: f'{a[0]} + {a[1]} * T + {a[2]} * T_2 + {a[3]} * T_3 + {a[4]} * T_4')}
}}
void fg_enthalpy_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {{
 {thermodynamics(lambda a: f'{a[0]} + {a[1]/2} * T + {a[2]/3} * T_2 + {a[3]/4} * T_3 + {a[4]/5} * T_4 + {a[5]} * rcp_T')}
}}
void fg_exp_Gibbs_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {{
 {thermodynamics(lambda a: f'exp2({a[5]/ln(2)} * rcp_T + {(a[0] - a[6])/ln(2)} + {-a[0]} * log_T + {-a[1]/2/ln(2)} * T + {(1./3.-1./2.)*a[2]/ln(2)} * T_2 + {(1./4.-1./3.)*a[3]/ln(2)} * T_3 + {(1./5.-1./4.)*a[4]/ln(2)} * T_4)')}
}}
float fg_viscosity_T_12(float ln_T, float ln_T_2, float ln_T_3, const float mole_fractions[]) {{
    {code([f'float sqrt_viscosity_T14_{k} = {P[0]} + {P[1]}*ln_T + {P[2]}*ln_T_2 + {P[3]}*ln_T_3;' for k, P in enumerate(transport_polynomials.sqrt_viscosity_T14)])}
 return
    {('+'+line).join([f'''mole_fractions[{k}] * sq(sqrt_viscosity_T14_{k}) / (
        {('+'+line+'	').join([(lambda sqrt_a: f'mole_fractions[{j}] * sq({sqrt_a} + {sqrt_a*sqrt(sqrt(species.molar_mass[j]/species.molar_mass[k]))} * sqrt_viscosity_T14_{k}/sqrt_viscosity_T14_{j})')(sqrt(1/sqrt(8) * 1/sqrt(1. + species.molar_mass[k]/species.molar_mass[j]))) for j in range(species.len)])}
    )''' for k in range(species.len)])};
}}

float fg_thermal_conductivity_T_12_2(float ln_T, float ln_T_2, float ln_T_3, const float mole_fractions[]) {{
 {code([f'float conductivity_T12_{k} = {P[0]} + {P[1]}*ln_T + {P[2]}*ln_T_2 + {P[3]}*ln_T_3;' for k, P in enumerate(transport_polynomials.thermal_conductivity_T12)])}
 return (
    {'+'.join([f'mole_fractions[{k}]*conductivity_T12_{k}' for k in range(species.len)])})
 + 1./ (
    {'+'.join([f'mole_fractions[{k}]/conductivity_T12_{k}' for k in range(species.len)])});
}}

void fg_P_T_32_mixture_diffusion_coefficients(float ln_T, float ln_T_2, float ln_T_3, const float mole_fractions[], const float mass_fractions[], float* _) {{
 {code(f'''_[{k}] = (1. - mass_fractions[{k}]) / (
    {('+'+line).join([(lambda P: f"mole_fractions[{j}] / ({P[0]} + {P[1]}*ln_T + {P[2]}*ln_T_2 + {P[3]}*ln_T_3)")(transport_polynomials.binary_thermal_diffusion_coefficients_T32[k if k>j else j][j if k>j else k]) for j in list(range(k))+list(range(k+1,species.len))])});''' for k in range(species.len))}
}}

void fg_rates(const float log_T, const float T, const float T_2, const float T_4, const float rcp_T, const float rcp_T2, const float P0_RT, const float rcp_P0_RT, const float exp_Gibbs0_RT[], const float concentrations[], float* _) {{
 float c, Pr, logFcent, logPr_c, f1;
    {code([reaction(i, r) for i, r in enumerate(reactions)])}
    {code([f"_[{specie}] = {'+'.join(filter(None, [mul(r.net[specie],f'cR{i}') for i, r in enumerate(reactions)]))};" for specie in range(species.len-1)])}
}}
""".replace('+ -','- ').replace('+-','-').replace('concentrations','C').replace('exp_Gibbs0_RT','G'))
