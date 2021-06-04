#define n_species 9
const dfloat fg_molar_mass[9] = {2.01594, 31.9988, 18.01534, 1.00797, 15.9994, 17.00737, 33.00677, 34.01474, 28.0134};
const dfloat fg_rcp_molar_mass[9] = {0.496046509321, 0.0312511719189, 0.0555082501912, 0.992093018641, 0.0625023438379, 0.0587980387326, 0.0302968148656, 0.0293990193663, 0.0356972020533};
void fg_molar_heat_capacity_at_constant_pressure_R(const dfloat Ts[], dfloat* species) {
    if (Ts[1] < 5000.0) {
        species[0] = +2.99142337e+00 +7.00064411e-04 * Ts[1] -5.63382869e-08 * Ts[2] -9.23157818e-12 * Ts[3] +1.58275179e-15 * Ts[4];
        species[1] = +3.69757819e+00 +6.13519689e-04 * Ts[1] -1.25884199e-07 * Ts[2] +1.77528148e-11 * Ts[3] -1.13643531e-15 * Ts[4];
        species[2] = +2.67214561e+00 +3.05629289e-03 * Ts[1] -8.73026011e-07 * Ts[2] +1.20099639e-10 * Ts[3] -6.39161787e-15 * Ts[4];
        species[3] = +2.50000000e+00 +0.00000000e+00 * Ts[1] +0.00000000e+00 * Ts[2] +0.00000000e+00 * Ts[3] +0.00000000e+00 * Ts[4];
        species[4] = +2.54205966e+00 -2.75506191e-05 * Ts[1] -3.10280335e-09 * Ts[2] +4.55106742e-12 * Ts[3] -4.36805150e-16 * Ts[4];
        species[7] = +4.57316685e+00 +4.33613639e-03 * Ts[1] -1.47468882e-06 * Ts[2] +2.34890357e-10 * Ts[3] -1.43165356e-14 * Ts[4];
        species[8] = +2.92664000e+00 +1.48797700e-03 * Ts[1] -5.68476100e-07 * Ts[2] +1.00970400e-10 * Ts[3] -6.75335100e-15 * Ts[4];
    } else {
        species[0] = +3.29812431e+00 +8.24944174e-04 * Ts[1] -8.14301529e-07 * Ts[2] -9.47543433e-11 * Ts[3] +4.13487224e-13 * Ts[4];
        species[1] = +3.21293640e+00 +1.12748635e-03 * Ts[1] -5.75615047e-07 * Ts[2] +1.31387723e-09 * Ts[3] -8.76855392e-13 * Ts[4];
        species[2] = +3.38684249e+00 +3.47498246e-03 * Ts[1] -6.35469633e-06 * Ts[2] +6.96858127e-09 * Ts[3] -2.50658847e-12 * Ts[4];
        species[3] = +2.50000000e+00 +0.00000000e+00 * Ts[1] +0.00000000e+00 * Ts[2] +0.00000000e+00 * Ts[3] +0.00000000e+00 * Ts[4];
        species[4] = +2.94642878e+00 -1.63816649e-03 * Ts[1] +2.42103170e-06 * Ts[2] -1.60284319e-09 * Ts[3] +3.89069636e-13 * Ts[4];
        species[7] = +3.38875365e+00 +6.56922581e-03 * Ts[1] -1.48501258e-07 * Ts[2] -4.62580552e-09 * Ts[3] +2.47151475e-12 * Ts[4];
        species[8] = +3.29867700e+00 +1.40824000e-03 * Ts[1] -3.96322200e-06 * Ts[2] +5.64151500e-09 * Ts[3] -2.44485500e-12 * Ts[4];
    }
    if (Ts[1] < 6000.0) {
        species[5] = +2.86472886e+00 +1.05650448e-03 * Ts[1] -2.59082758e-07 * Ts[2] +3.05218674e-11 * Ts[3] -1.33195876e-15 * Ts[4];
    } else {
        species[5] = +4.12530561e+00 -3.22544939e-03 * Ts[1] +6.52764691e-06 * Ts[2] -5.79853643e-09 * Ts[3] +2.06237379e-12 * Ts[4];
    }
    if (Ts[1] < 3500.0) {
        species[6] = +4.01721090e+00 +2.23982013e-03 * Ts[1] -6.33658150e-07 * Ts[2] +1.14246370e-10 * Ts[3] -1.07908535e-14 * Ts[4];
    } else {
        species[6] = +4.30179801e+00 -4.74912051e-03 * Ts[1] +2.11582891e-05 * Ts[2] -2.42763894e-08 * Ts[3] +9.29225124e-12 * Ts[4];
    }
}
void fg_enthalpy_RT(const dfloat Ts[], dfloat* species) {
    if (Ts[1] < 5000.0) {
        species[0] = +2.99142337e+00 +3.50032206e-04 * Ts[1] -1.87794290e-08 * Ts[2] -2.30789455e-12 * Ts[3] +3.16550358e-16 * Ts[4] -8.35033997e+02 * Ts[5];
        species[1] = +3.69757819e+00 +3.06759845e-04 * Ts[1] -4.19613997e-08 * Ts[2] +4.43820370e-12 * Ts[3] -2.27287062e-16 * Ts[4] -1.23393018e+03 * Ts[5];
        species[2] = +2.67214561e+00 +1.52814644e-03 * Ts[1] -2.91008670e-07 * Ts[2] +3.00249098e-11 * Ts[3] -1.27832357e-15 * Ts[4] -2.98992090e+04 * Ts[5];
        species[3] = +2.50000000e+00 +0.00000000e+00 * Ts[1] +0.00000000e+00 * Ts[2] +0.00000000e+00 * Ts[3] +0.00000000e+00 * Ts[4] +2.54716270e+04 * Ts[5];
        species[4] = +2.54205966e+00 -1.37753096e-05 * Ts[1] -1.03426778e-09 * Ts[2] +1.13776685e-12 * Ts[3] -8.73610300e-17 * Ts[4] +2.92308027e+04 * Ts[5];
        species[7] = +4.57316685e+00 +2.16806820e-03 * Ts[1] -4.91562940e-07 * Ts[2] +5.87225893e-11 * Ts[3] -2.86330712e-15 * Ts[4] -1.80069609e+04 * Ts[5];
        species[8] = +2.92664000e+00 +7.43988500e-04 * Ts[1] -1.89492033e-07 * Ts[2] +2.52426000e-11 * Ts[3] -1.35067020e-15 * Ts[4] -9.22797700e+02 * Ts[5];
    } else {
        species[0] = +3.29812431e+00 +4.12472087e-04 * Ts[1] -2.71433843e-07 * Ts[2] -2.36885858e-11 * Ts[3] +8.26974448e-14 * Ts[4] -1.01252087e+03 * Ts[5];
        species[1] = +3.21293640e+00 +5.63743175e-04 * Ts[1] -1.91871682e-07 * Ts[2] +3.28469308e-10 * Ts[3] -1.75371078e-13 * Ts[4] -1.00524902e+03 * Ts[5];
        species[2] = +3.38684249e+00 +1.73749123e-03 * Ts[1] -2.11823211e-06 * Ts[2] +1.74214532e-09 * Ts[3] -5.01317694e-13 * Ts[4] -3.02081133e+04 * Ts[5];
        species[3] = +2.50000000e+00 +0.00000000e+00 * Ts[1] +0.00000000e+00 * Ts[2] +0.00000000e+00 * Ts[3] +0.00000000e+00 * Ts[4] +2.54716270e+04 * Ts[5];
        species[4] = +2.94642878e+00 -8.19083245e-04 * Ts[1] +8.07010567e-07 * Ts[2] -4.00710797e-10 * Ts[3] +7.78139272e-14 * Ts[4] +2.91476445e+04 * Ts[5];
        species[7] = +3.38875365e+00 +3.28461290e-03 * Ts[1] -4.95004193e-08 * Ts[2] -1.15645138e-09 * Ts[3] +4.94302950e-13 * Ts[4] -1.76631465e+04 * Ts[5];
        species[8] = +3.29867700e+00 +7.04120000e-04 * Ts[1] -1.32107400e-06 * Ts[2] +1.41037875e-09 * Ts[3] -4.88971000e-13 * Ts[4] -1.02090000e+03 * Ts[5];
    }
    if (Ts[1] < 6000.0) {
        species[5] = +2.86472886e+00 +5.28252240e-04 * Ts[1] -8.63609193e-08 * Ts[2] +7.63046685e-12 * Ts[3] -2.66391752e-16 * Ts[4] +3.68362875e+03 * Ts[5];
    } else {
        species[5] = +4.12530561e+00 -1.61272470e-03 * Ts[1] +2.17588230e-06 * Ts[2] -1.44963411e-09 * Ts[3] +4.12474758e-13 * Ts[4] +3.34630913e+03 * Ts[5];
    }
    if (Ts[1] < 3500.0) {
        species[6] = +4.01721090e+00 +1.11991006e-03 * Ts[1] -2.11219383e-07 * Ts[2] +2.85615925e-11 * Ts[3] -2.15817070e-15 * Ts[4] +1.11856713e+02 * Ts[5];
    } else {
        species[6] = +4.30179801e+00 -2.37456025e-03 * Ts[1] +7.05276303e-06 * Ts[2] -6.06909735e-09 * Ts[3] +1.85845025e-12 * Ts[4] +2.94808040e+02 * Ts[5];
    }
}
void fg_gibbs_RT(const dfloat Ts[], dfloat* species) {
    if (Ts[1] < 5000.0) {
        species[0] = -8.35033997e+02 * Ts[5] +4.34653354e+00 -2.99142337e+00 * Ts[0] +3.50032206e-04 * Ts[1] +9.38971448e-09 * Ts[2] +7.69298182e-13 * Ts[3] -7.91375895e-17 * Ts[4];
        species[1] = -1.23393018e+03 * Ts[5] +5.08412600e-01 -3.69757819e+00 * Ts[0] +3.06759845e-04 * Ts[1] +2.09806998e-08 * Ts[2] -1.47940123e-12 * Ts[3] +5.68217655e-17 * Ts[4];
        species[2] = -2.98992090e+04 * Ts[5] -4.19067120e+00 -2.67214561e+00 * Ts[0] +1.52814644e-03 * Ts[1] +1.45504335e-07 * Ts[2] -1.00083033e-11 * Ts[3] +3.19580894e-16 * Ts[4];
        species[3] = +2.54716270e+04 * Ts[5] +2.96011764e+00 -2.50000000e+00 * Ts[0] +0.00000000e+00 * Ts[1] -0.00000000e+00 * Ts[2] -0.00000000e+00 * Ts[3] -0.00000000e+00 * Ts[4];
        species[4] = +2.92308027e+04 * Ts[5] -2.37824845e+00 -2.54205966e+00 * Ts[0] -1.37753096e-05 * Ts[1] +5.17133892e-10 * Ts[2] -3.79255618e-13 * Ts[3] +2.18402575e-17 * Ts[4];
        species[7] = -1.80069609e+04 * Ts[5] +4.07202989e+00 -4.57316685e+00 * Ts[0] +2.16806820e-03 * Ts[1] +2.45781470e-07 * Ts[2] -1.95741964e-11 * Ts[3] +7.15826780e-16 * Ts[4];
        species[8] = -9.22797700e+02 * Ts[5] -3.05388800e+00 -2.92664000e+00 * Ts[0] +7.43988500e-04 * Ts[1] +9.47460167e-08 * Ts[2] -8.41420000e-12 * Ts[3] +3.37667550e-16 * Ts[4];
    } else {
        species[0] = -1.01252087e+03 * Ts[5] +6.59221840e+00 -3.29812431e+00 * Ts[0] +4.12472087e-04 * Ts[1] +1.35716922e-07 * Ts[2] +7.89619527e-12 * Ts[3] -2.06743612e-14 * Ts[4];
        species[1] = -1.00524902e+03 * Ts[5] -2.82180119e+00 -3.21293640e+00 * Ts[0] +5.63743175e-04 * Ts[1] +9.59358412e-08 * Ts[2] -1.09489769e-10 * Ts[3] +4.38427696e-14 * Ts[4];
        species[2] = -3.02081133e+04 * Ts[5] +7.96609640e-01 -3.38684249e+00 * Ts[0] +1.73749123e-03 * Ts[1] +1.05911606e-06 * Ts[2] -5.80715106e-10 * Ts[3] +1.25329424e-13 * Ts[4];
        species[3] = +2.54716270e+04 * Ts[5] +2.96011761e+00 -2.50000000e+00 * Ts[0] +0.00000000e+00 * Ts[1] -0.00000000e+00 * Ts[2] -0.00000000e+00 * Ts[3] -0.00000000e+00 * Ts[4];
        species[4] = +2.91476445e+04 * Ts[5] -1.75662000e-02 -2.94642878e+00 * Ts[0] -8.19083245e-04 * Ts[1] -4.03505283e-07 * Ts[2] +1.33570266e-10 * Ts[3] -1.94534818e-14 * Ts[4];
        species[7] = -1.76631465e+04 * Ts[5] -3.39660955e+00 -3.38875365e+00 * Ts[0] +3.28461290e-03 * Ts[1] +2.47502097e-08 * Ts[2] +3.85483793e-10 * Ts[3] -1.23575738e-13 * Ts[4];
        species[8] = -1.02090000e+03 * Ts[5] -6.51695000e-01 -3.29867700e+00 * Ts[0] +7.04120000e-04 * Ts[1] +6.60537000e-07 * Ts[2] -4.70126250e-10 * Ts[3] +1.22242750e-13 * Ts[4];
    }
    if (Ts[1] < 6000.0) {
        species[5] = +3.68362875e+03 * Ts[5] -2.83691187e+00 -2.86472886e+00 * Ts[0] +5.28252240e-04 * Ts[1] +4.31804597e-08 * Ts[2] -2.54348895e-12 * Ts[3] +6.65979380e-17 * Ts[4];
    } else {
        species[5] = +3.34630913e+03 * Ts[5] +4.81573857e+00 -4.12530561e+00 * Ts[0] -1.61272470e-03 * Ts[1] -1.08794115e-06 * Ts[2] +4.83211369e-10 * Ts[3] -1.03118689e-13 * Ts[4];
    }
    if (Ts[1] < 3500.0) {
        species[6] = +1.11856713e+02 * Ts[5] +2.32108750e-01 -4.01721090e+00 * Ts[0] +1.11991006e-03 * Ts[1] +1.05609692e-07 * Ts[2] -9.52053083e-12 * Ts[3] +5.39542675e-16 * Ts[4];
    } else {
        species[6] = +2.94808040e+02 * Ts[5] +5.85135560e-01 -4.30179801e+00 * Ts[0] -2.37456025e-03 * Ts[1] -3.52638152e-06 * Ts[2] +2.02303245e-09 * Ts[3] -4.64612562e-13 * Ts[4];
    }
}
dfloat fg_viscosity(dfloat T, const dfloat mole_fractions[]) {
    dfloat T_14 =	sqrt(sqrt(T));
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((-8.43124438E-29 + -2.81221390E-28 * ln_T + -6.77956489E-28 * ln_T_2 + 1.43298122E-28 * ln_T_3)*T_14, 12.)
    + mole_fractions[1] * pow((-2.04236623E-28 + -7.03273352E-28 * ln_T + -1.86569998E-27 * ln_T_2 + -1.33706385E-27 * ln_T_3)*T_14, 12.)
    + mole_fractions[2] * pow((-2.12764345E-72 + -7.84145313E-72 * ln_T + -2.46555433E-71 * ln_T_2 + -5.32725773E-71 * ln_T_3)*T_14, 12.)
    + mole_fractions[3] * pow((-9.87847960E-29 + -3.43034438E-28 * ln_T + -9.31547260E-28 * ln_T_2 + -8.66426326E-28 * ln_T_3)*T_14, 12.)
    + mole_fractions[4] * pow((-2.38512773E-28 + -8.14279582E-28 * ln_T + -2.10766357E-27 * ln_T_2 + -1.02515451E-27 * ln_T_3)*T_14, 12.)
    + mole_fractions[5] * pow((-2.45911236E-28 + -8.39537843E-28 * ln_T + -2.17304150E-27 * ln_T_2 + -1.05695393E-27 * ln_T_3)*T_14, 12.)
    + mole_fractions[6] * pow((-2.07428434E-28 + -7.14264111E-28 * ln_T + -1.89485714E-27 * ln_T_2 + -1.35795949E-27 * ln_T_3)*T_14, 12.)
    + mole_fractions[7] * pow((-2.10571869E-28 + -7.25088293E-28 * ln_T + -1.92357240E-27 * ln_T_2 + -1.37853843E-27 * ln_T_3)*T_14, 12.)
    + mole_fractions[8] * pow((-1.76755675E-28 + -6.06958851E-28 * ln_T + -1.59757718E-27 * ln_T_2 + -1.02837655E-27 * ln_T_3)*T_14, 12.)
    , 1./6.);
}
dfloat fg_thermal_conductivity(dfloat T, const dfloat mole_fractions[]) {
    dfloat T_12 = sqrt(T)
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((-2.23320216E-28 + -1.14772401E-27 * ln_T + -5.87799312E-27 * ln_T_2 + -3.03916863E-26 * ln_T_3)*T_12, 4.)
    + mole_fractions[1] * pow((-2.18878344E-29 + -2.01946778E-28 * ln_T + -1.48262867E-27 * ln_T_2 + -9.81186276E-27 * ln_T_3)*T_12, 4.)
    + mole_fractions[2] * pow((-3.59810621E-29 + -3.76015847E-29 * ln_T + 5.45804495E-28 * ln_T_2 + 6.35615723E-27 * ln_T_3)*T_12, 4.)
    + mole_fractions[3] * pow((4.06603914E-29 + -1.53712163E-27 * ln_T + -1.66239047E-26 * ln_T_2 + -1.27840726E-25 * ln_T_3)*T_12, 4.)
    + mole_fractions[4] * pow((5.60932767E-30 + -2.32064800E-28 * ln_T + -2.49613606E-27 * ln_T_2 + -1.91648175E-26 * ln_T_3)*T_12, 4.)
    + mole_fractions[5] * pow((-7.52789207E-29 + -3.97495239E-28 * ln_T + -2.08892157E-27 * ln_T_2 + -1.10551162E-26 * ln_T_3)*T_12, 4.)
    + mole_fractions[6] * pow((-7.07042966E-29 + -3.01378413E-28 * ln_T + -1.23275378E-27 * ln_T_2 + -4.88657015E-27 * ln_T_3)*T_12, 4.)
    + mole_fractions[7] * pow((-9.93405991E-29 + -3.93326612E-28 * ln_T + -1.42686938E-27 * ln_T_2 + -4.56541721E-27 * ln_T_3)*T_12, 4.)
    + mole_fractions[8] * pow((-2.98243791E-29 + -1.81178461E-28 * ln_T + -1.06773262E-27 * ln_T_2 + -6.18993987E-27 * ln_T_3)*T_12, 4.)
    , 1./4.);
}
void fg_mixture_diffusion_coefficients(const dfloat mole_fractions[n_species], const dfloat mass_fractions[n_species], dfloat T, dfloat* Ddiag) {
    dfloat T_12 = sqrt(T)
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    dfloat T_32 = T*T_12;
    Ddiag[0] = (1. - mass_fractions[0]) * mole_fractions[0] / ( 0.
    + mole_fractions[1] / ((1.12118490E-29 + -1.71594299E-28 * ln_T + -2.02767051E-27 * ln_T_2 + -1.59827303E-26 * ln_T_3)*T_32)
    + mole_fractions[2] / ((1.67469938E-29 + -2.17803488E-28 * ln_T + -2.63851525E-27 * ln_T_2 + -2.09320359E-26 * ln_T_3)*T_32)
    + mole_fractions[3] / ((3.16346812E-29 + -4.68461507E-28 * ln_T + -5.56206992E-27 * ln_T_2 + -4.38967779E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((1.42973456E-29 + -2.26454889E-28 * ln_T + -2.66308388E-27 * ln_T_2 + -2.09645889E-26 * ln_T_3)*T_32)
    + mole_fractions[5] / ((1.42498566E-29 + -2.25702714E-28 * ln_T + -2.65423839E-27 * ln_T_2 + -2.08949546E-26 * ln_T_3)*T_32)
    + mole_fractions[6] / ((1.12016983E-29 + -1.71438944E-28 * ln_T + -2.02583474E-27 * ln_T_2 + -1.59682601E-26 * ln_T_3)*T_32)
    + mole_fractions[7] / ((1.11921407E-29 + -1.71292668E-28 * ln_T + -2.02410624E-27 * ln_T_2 + -1.59546356E-26 * ln_T_3)*T_32)
    + mole_fractions[8] / ((1.06335867E-29 + -1.64545141E-28 * ln_T + -1.94134182E-27 * ln_T_2 + -1.52959726E-26 * ln_T_3)*T_32)
    );
    Ddiag[1] = (1. - mass_fractions[1]) * mole_fractions[1] / ( 0.
    + mole_fractions[0] / ((1.12118490E-29 + -1.71594299E-28 * ln_T + -2.02767051E-27 * ln_T_2 + -1.59827303E-26 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.73575478E-30 + -6.88536746E-29 * ln_T + -8.45481726E-28 * ln_T_2 + -6.73044504E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.26027830E-29 + -3.02592794E-28 * ln_T + -3.64857378E-27 * ln_T_2 + -2.89105021E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.22422501E-30 + -7.39246747E-29 * ln_T + -8.83692461E-28 * ln_T_2 + -6.98659487E-27 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.11997744E-30 + -7.24495339E-29 * ln_T + -8.66058690E-28 * ln_T_2 + -6.84717984E-27 * ln_T_3)*T_32)
    + mole_fractions[6] / ((3.46446413E-30 + -4.76616380E-29 * ln_T + -5.72224199E-28 * ln_T_2 + -4.52916710E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((3.43910338E-30 + -4.73127428E-29 * ln_T + -5.68035373E-28 * ln_T_2 + -4.49601245E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.43160456E-30 + -4.76391732E-29 * ln_T + -5.71150114E-28 * ln_T_2 + -4.51902339E-27 * ln_T_3)*T_32)
    );
    Ddiag[2] = (1. - mass_fractions[2]) * mole_fractions[2] / ( 0.
    + mole_fractions[0] / ((1.67469938E-29 + -2.17803488E-28 * ln_T + -2.63851525E-27 * ln_T_2 + -2.09320359E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.73575478E-30 + -6.88536746E-29 * ln_T + -8.45481726E-28 * ln_T_2 + -6.73044504E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((3.41281989E-29 + -4.01254458E-28 * ln_T + -4.94525371E-27 * ln_T_2 + -3.94027333E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((8.48774448E-30 + -1.04109554E-28 * ln_T + -1.27363857E-27 * ln_T_2 + -1.01292709E-26 * ln_T_3)*T_32)
    + mole_fractions[5] / ((8.35346892E-30 + -1.02462547E-28 * ln_T + -1.25348969E-27 * ln_T_2 + -9.96902648E-27 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.70303533E-30 + -6.84688808E-29 * ln_T + -8.40739566E-28 * ln_T_2 + -6.69266096E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.67312088E-30 + -6.81097372E-29 * ln_T + -8.36329587E-28 * ln_T_2 + -6.65755557E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((5.55911099E-30 + -6.71959768E-29 * ln_T + -8.24133133E-28 * ln_T_2 + -6.55851747E-27 * ln_T_3)*T_32)
    );
    Ddiag[3] = (1. - mass_fractions[3]) * mole_fractions[3] / ( 0.
    + mole_fractions[0] / ((3.16346812E-29 + -4.68461507E-28 * ln_T + -5.56206992E-27 * ln_T_2 + -4.38967779E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((2.26027830E-29 + -3.02592794E-28 * ln_T + -3.64857378E-27 * ln_T_2 + -2.89105021E-26 * ln_T_3)*T_32)
    + mole_fractions[2] / ((3.41281989E-29 + -4.01254458E-28 * ln_T + -4.94525371E-27 * ln_T_2 + -3.94027333E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((2.97799978E-29 + -4.09476935E-28 * ln_T + -4.91657063E-27 * ln_T_2 + -3.89155870E-26 * ln_T_3)*T_32)
    + mole_fractions[5] / ((2.97276501E-29 + -4.08757151E-28 * ln_T + -4.90792822E-27 * ln_T_2 + -3.88471806E-26 * ln_T_3)*T_32)
    + mole_fractions[6] / ((2.25922410E-29 + -3.02451664E-28 * ln_T + -3.64687207E-27 * ln_T_2 + -2.88970182E-26 * ln_T_3)*T_32)
    + mole_fractions[7] / ((2.25823193E-29 + -3.02318839E-28 * ln_T + -3.64527050E-27 * ln_T_2 + -2.88843277E-26 * ln_T_3)*T_32)
    + mole_fractions[8] / ((2.12702557E-29 + -2.87213061E-28 * ln_T + -3.45839739E-27 * ln_T_2 + -2.73939657E-26 * ln_T_3)*T_32)
    );
    Ddiag[4] = (1. - mass_fractions[4]) * mole_fractions[4] / ( 0.
    + mole_fractions[0] / ((1.42973456E-29 + -2.26454889E-28 * ln_T + -2.66308388E-27 * ln_T_2 + -2.09645889E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.22422501E-30 + -7.39246747E-29 * ln_T + -8.83692461E-28 * ln_T_2 + -6.98659487E-27 * ln_T_3)*T_32)
    + mole_fractions[2] / ((8.48774448E-30 + -1.04109554E-28 * ln_T + -1.27363857E-27 * ln_T_2 + -1.01292709E-26 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.97799978E-29 + -4.09476935E-28 * ln_T + -4.91657063E-27 * ln_T_2 + -3.89155870E-26 * ln_T_3)*T_32)
    + mole_fractions[5] / ((7.44309499E-30 + -1.08499161E-28 * ln_T + -1.29121163E-27 * ln_T_2 + -1.01966396E-26 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.19756719E-30 + -7.35474569E-29 * ln_T + -8.79183215E-28 * ln_T_2 + -6.95094415E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.17236387E-30 + -7.31908208E-29 * ln_T + -8.74920002E-28 * ln_T_2 + -6.91723860E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((5.04917663E-30 + -7.21332094E-29 * ln_T + -8.61029423E-28 * ln_T_2 + -6.80485874E-27 * ln_T_3)*T_32)
    );
    Ddiag[5] = (1. - mass_fractions[5]) * mole_fractions[5] / ( 0.
    + mole_fractions[0] / ((1.42498566E-29 + -2.25702714E-28 * ln_T + -2.65423839E-27 * ln_T_2 + -2.08949546E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.11997744E-30 + -7.24495339E-29 * ln_T + -8.66058690E-28 * ln_T_2 + -6.84717984E-27 * ln_T_3)*T_32)
    + mole_fractions[2] / ((8.35346892E-30 + -1.02462547E-28 * ln_T + -1.25348969E-27 * ln_T_2 + -9.96902648E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.97276501E-29 + -4.08757151E-28 * ln_T + -4.90792822E-27 * ln_T_2 + -3.88471806E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((7.44309499E-30 + -1.08499161E-28 * ln_T + -1.29121163E-27 * ln_T_2 + -1.01966396E-26 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.09277397E-30 + -7.20645950E-29 * ln_T + -8.61457146E-28 * ln_T_2 + -6.81079939E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.06704944E-30 + -7.17005835E-29 * ln_T + -8.57105769E-28 * ln_T_2 + -6.77639680E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((4.95302810E-30 + -7.07596187E-29 * ln_T + -8.44633341E-28 * ln_T_2 + -6.67527777E-27 * ln_T_3)*T_32)
    );
    Ddiag[6] = (1. - mass_fractions[6]) * mole_fractions[6] / ( 0.
    + mole_fractions[0] / ((1.12016983E-29 + -1.71438944E-28 * ln_T + -2.02583474E-27 * ln_T_2 + -1.59682601E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((3.46446413E-30 + -4.76616380E-29 * ln_T + -5.72224199E-28 * ln_T_2 + -4.52916710E-27 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.70303533E-30 + -6.84688808E-29 * ln_T + -8.40739566E-28 * ln_T_2 + -6.69266096E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.25922410E-29 + -3.02451664E-28 * ln_T + -3.64687207E-27 * ln_T_2 + -2.88970182E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.19756719E-30 + -7.35474569E-29 * ln_T + -8.79183215E-28 * ln_T_2 + -6.95094415E-27 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.09277397E-30 + -7.20645950E-29 * ln_T + -8.61457146E-28 * ln_T_2 + -6.81079939E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((3.41193820E-30 + -4.69390235E-29 * ln_T + -5.63548511E-28 * ln_T_2 + -4.46049884E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.40705780E-30 + -4.72984033E-29 * ln_T + -5.67064596E-28 * ln_T_2 + -4.48669817E-27 * ln_T_3)*T_32)
    );
    Ddiag[7] = (1. - mass_fractions[7]) * mole_fractions[7] / ( 0.
    + mole_fractions[0] / ((1.11921407E-29 + -1.71292668E-28 * ln_T + -2.02410624E-27 * ln_T_2 + -1.59546356E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((3.43910338E-30 + -4.73127428E-29 * ln_T + -5.68035373E-28 * ln_T_2 + -4.49601245E-27 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.67312088E-30 + -6.81097372E-29 * ln_T + -8.36329587E-28 * ln_T_2 + -6.65755557E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.25823193E-29 + -3.02318839E-28 * ln_T + -3.64527050E-27 * ln_T_2 + -2.88843277E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.17236387E-30 + -7.31908208E-29 * ln_T + -8.74920002E-28 * ln_T_2 + -6.91723860E-27 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.06704944E-30 + -7.17005835E-29 * ln_T + -8.57105769E-28 * ln_T_2 + -6.77639680E-27 * ln_T_3)*T_32)
    + mole_fractions[6] / ((3.41193820E-30 + -4.69390235E-29 * ln_T + -5.63548511E-28 * ln_T_2 + -4.46049884E-27 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.38380330E-30 + -4.69755733E-29 * ln_T + -5.63194158E-28 * ln_T_2 + -4.45607470E-27 * ln_T_3)*T_32)
    );
    Ddiag[8] = (1. - mass_fractions[8]) * mole_fractions[8] / ( 0.
    + mole_fractions[0] / ((1.06335867E-29 + -1.64545141E-28 * ln_T + -1.94134182E-27 * ln_T_2 + -1.52959726E-26 * ln_T_3)*T_32)
    + mole_fractions[1] / ((3.43160456E-30 + -4.76391732E-29 * ln_T + -5.71150114E-28 * ln_T_2 + -4.51902339E-27 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.55911099E-30 + -6.71959768E-29 * ln_T + -8.24133133E-28 * ln_T_2 + -6.55851747E-27 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.12702557E-29 + -2.87213061E-28 * ln_T + -3.45839739E-27 * ln_T_2 + -2.73939657E-26 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.04917663E-30 + -7.21332094E-29 * ln_T + -8.61029423E-28 * ln_T_2 + -6.80485874E-27 * ln_T_3)*T_32)
    + mole_fractions[5] / ((4.95302810E-30 + -7.07596187E-29 * ln_T + -8.44633341E-28 * ln_T_2 + -6.67527777E-27 * ln_T_3)*T_32)
    + mole_fractions[6] / ((3.40705780E-30 + -4.72984033E-29 * ln_T + -5.67064596E-28 * ln_T_2 + -4.48669817E-27 * ln_T_3)*T_32)
    + mole_fractions[7] / ((3.38380330E-30 + -4.69755733E-29 * ln_T + -5.63194158E-28 * ln_T_2 + -4.45607470E-27 * ln_T_3)*T_32)
    );
}
void fg_rates(dfloat * concentrations, dfloat Ts[], dfloat * wdot) {
    const dfloat rcp_P0_RT = 8.2057366080959685e-05 * Ts[1];
    const dfloat P0_RT = 1.2186596374704217e+04 * Ts[5];
    const dfloat T = Ts[1];
    dfloat sum_concentrations = 0.0;
    for (int i = 0; i < 9; i++) { sum_concentrations += concentrations[i]; }
    dfloat gibbs0_RT[9];
    fg_gibbs_RT(gibbs0_RT, Ts);
    {
    dfloat phi_f = concentrations[3]*concentrations[1];
    dfloat k_f = 1.000000e-06 * 3.547000e+15*fgexp(-4.060000e-01*Ts[0]-8.352941e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[4]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[3] + gibbs0_RT[1]) - (gibbs0_RT[1] + gibbs0_RT[1])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[4]*concentrations[0];
    dfloat k_f = 1.000000e-06 * 5.080000e+04*fgexp(2.670000e+00*Ts[0]-3.165251e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[3]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[4] + gibbs0_RT[0]) - (gibbs0_RT[0] + gibbs0_RT[0])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[0]*concentrations[5];
    dfloat k_f = 1.000000e-06 * 2.160000e+08*fgexp(1.510000e+00*Ts[0]-1.726043e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[2]*concentrations[3];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[0] + gibbs0_RT[5]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    }
    {
    dfloat phi_f = concentrations[4]*concentrations[2];
    dfloat k_f = 1.000000e-06 * 2.970000e+06*fgexp(2.020000e+00*Ts[0]-6.743142e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[5]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[4] + gibbs0_RT[2]) - (gibbs0_RT[2] + gibbs0_RT[2])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[0];
    dfloat alpha = sum_concentrations + 1.500000e+00*concentrations[0] + 1.100000e+01*concentrations[2];
    dfloat k_f = 1.000000e-06 * alpha * 4.577000e+19*fgexp(-1.400000e+00*Ts[0]-5.252605e+04*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[3]*concentrations[3];
    dfloat rcp_Kc = rcp_P0_RT * fgexp(-((gibbs0_RT[0]) - (gibbs0_RT[0] + gibbs0_RT[0])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    }
    {
    dfloat phi_f = concentrations[4]*concentrations[4];
    dfloat alpha = sum_concentrations + 1.500000e+00*concentrations[0] + 1.100000e+01*concentrations[2];
    dfloat k_f = 1.000000e-12 * alpha * 6.165000e+15*fgexp(-5.000000e-01*Ts[0]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[1];
    dfloat rcp_Kc = P0_RT * fgexp(-((gibbs0_RT[4] + gibbs0_RT[4]) - (gibbs0_RT[4])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[1] += qdot;
    }
    {
    dfloat phi_f = concentrations[4]*concentrations[3];
    dfloat alpha = sum_concentrations + 1.500000e+00*concentrations[0] + 1.100000e+01*concentrations[2];
    dfloat k_f = 1.000000e-12 * alpha * 4.714000e+18*fgexp(-1.000000e+00*Ts[0]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[5];
    dfloat rcp_Kc = P0_RT * fgexp(-((gibbs0_RT[4] + gibbs0_RT[3]) - (gibbs0_RT[3])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[3]*concentrations[5];
    dfloat alpha = sum_concentrations + 1.500000e+00*concentrations[0] + 1.100000e+01*concentrations[2];
    dfloat k_f = 1.000000e-12 * alpha * 3.800000e+22*fgexp(-2.000000e+00*Ts[0]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[2];
    dfloat rcp_Kc = P0_RT * fgexp(-((gibbs0_RT[3] + gibbs0_RT[5]) - (gibbs0_RT[5])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;
    }
    {
    dfloat phi_f = concentrations[3]*concentrations[1];
    dfloat alpha = sum_concentrations + concentrations[0] + 1.000000e+01*concentrations[2] + -2.200000e-01*concentrations[1];
    dfloat k_f = 1.000000e-06 * 1.475000e+12*fgexp(6.000000e-01*Ts[0]);
    dfloat redP = 1.0e-12 * alpha / k_f * 6.366000e+20*fgexp(-1.720000e+00*Ts[0]-2.640896e+02*Ts[5]);
    dfloat F = redP / (1 + redP);
    dfloat logPred = log10(redP);
    dfloat logFcent = log10((2.000000e-01*fgexp(T*(-1.000000e+30)))+ (8.000000e-01*fgexp(T*(-1.000000e-30))));
    dfloat troe_c = -.4 - .67 * logFcent;
    dfloat troe_n = .75 - 1.27 * logFcent;
    dfloat troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    dfloat F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[6];
    dfloat rcp_Kc = P0_RT * fgexp(-((gibbs0_RT[3] + gibbs0_RT[1]) - (gibbs0_RT[1])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[3];
    dfloat k_f = 1.000000e-06 * 1.660000e+13*fgexp(-4.141497e+02*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[0]*concentrations[1];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[3]) - (gibbs0_RT[3] + gibbs0_RT[3])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[3] -= qdot;
    wdot[0] += qdot;
    wdot[1] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[3];
    dfloat k_f = 1.000000e-06 * 7.079000e+13*fgexp(-1.484498e+02*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[5]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[3]) - (gibbs0_RT[3] + gibbs0_RT[3])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[4];
    dfloat k_f = 1.000000e-06 * 3.250000e+13;
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[1]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[4]) - (gibbs0_RT[4] + gibbs0_RT[4])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[4] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[5];
    dfloat k_f = 1.000000e-06 * 2.890000e+13*fgexp(+2.501001e+02*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[2]*concentrations[1];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[5]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;
    wdot[1] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[6];
    dfloat k_f = 1.000000e-06 * 4.200000e+14*fgexp(-6.029576e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[7]*concentrations[1];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[6]) - (gibbs0_RT[6] + gibbs0_RT[6])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[1] += qdot;
    }
    {
    dfloat phi_f = concentrations[6]*concentrations[6];
    dfloat k_f = 1.000000e-06 * 1.300000e+11*fgexp(+8.198956e+02*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[7]*concentrations[1];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[6]) - (gibbs0_RT[6] + gibbs0_RT[6])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[1] += qdot;
    }
    {
    dfloat phi_f = concentrations[7];
    dfloat alpha = sum_concentrations + 1.500000e+00*concentrations[0] + 1.100000e+01*concentrations[2];
    dfloat k_f = 1.000000e+00 * 2.951000e+14*fgexp(-2.437092e+04*Ts[5]);
    dfloat redP = 1.0e-12 * alpha / k_f * 1.202000e+17*fgexp(-2.289649e+04*Ts[5]);
    dfloat F = redP / (1 + redP);
    dfloat logPred = log10(redP);
    dfloat logFcent = log10((5.000000e-01*fgexp(T*(-1.000000e+30)))+ (5.000000e-01*fgexp(T*(-1.000000e-30))));
    dfloat troe_c = -.4 - .67 * logFcent;
    dfloat troe_n = .75 - 1.27 * logFcent;
    dfloat troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    dfloat F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[5]*concentrations[5];
    dfloat rcp_Kc = rcp_P0_RT * fgexp(-((gibbs0_RT[7]) - (gibbs0_RT[7] + gibbs0_RT[7])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[7]*concentrations[3];
    dfloat k_f = 1.000000e-06 * 2.410000e+13*fgexp(-1.997782e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[2]*concentrations[5];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[3]) - (gibbs0_RT[3] + gibbs0_RT[3])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[3] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    }
    {
    dfloat phi_f = concentrations[7]*concentrations[3];
    dfloat k_f = 1.000000e-06 * 4.820000e+13*fgexp(-4.000595e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[6]*concentrations[0];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[3]) - (gibbs0_RT[3] + gibbs0_RT[3])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[0] += qdot;
    }
    {
    dfloat phi_f = concentrations[7]*concentrations[4];
    dfloat k_f = 1.000000e-06 * 9.550000e+06*fgexp(2.000000e+00*Ts[0]-1.997782e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[5]*concentrations[6];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[4]) - (gibbs0_RT[4] + gibbs0_RT[4])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    }
    {
    dfloat phi_f = concentrations[7]*concentrations[5];
    dfloat k_f = 1.000000e-06 * 1.000000e+12;
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[6]*concentrations[2];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[5]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[2] += qdot;
    }
    {
    dfloat phi_f = concentrations[7]*concentrations[5];
    dfloat k_f = 1.000000e-06 * 5.800000e+14*fgexp(-4.809269e+03*Ts[5]);
    dfloat q_f = phi_f * k_f;
    dfloat phi_r = concentrations[6]*concentrations[2];
    dfloat rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[5]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    dfloat k_r = k_f * rcp_Kc;
    dfloat q_r = phi_r * k_r;
    dfloat qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[2] += qdot;
    }
}
const int inert_specie= 8;
