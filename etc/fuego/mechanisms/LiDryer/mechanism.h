#define n_species 9
const dfloat fg_molar_mass[9] = {0.00201594, 0.0319988, 0.01801534, 0.00100797, 0.0159994, 0.01700737, 0.03300677, 0.03401474, 0.0280134};
const dfloat fg_rcp_molar_mass[9] = {496.046509321, 31.2511719189, 55.5082501912, 992.093018641, 62.5023438379, 58.7980387326, 30.2968148656, 29.3990193663, 35.6972020533};
void fg_molar_heat_capacity_at_constant_pressure_R(const dfloat log_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
    if (T < 1000.0) {
        species[0] = +3.29812431e+00 +8.24944174e-04 * T -8.14301529e-07 * T_2 -9.47543433e-11 * T_3 +4.13487224e-13 * T_4;
        species[1] = +3.21293640e+00 +1.12748635e-03 * T -5.75615047e-07 * T_2 +1.31387723e-09 * T_3 -8.76855392e-13 * T_4;
        species[2] = +3.38684249e+00 +3.47498246e-03 * T -6.35469633e-06 * T_2 +6.96858127e-09 * T_3 -2.50658847e-12 * T_4;
        species[3] = +2.50000000e+00 +0.00000000e+00 * T +0.00000000e+00 * T_2 +0.00000000e+00 * T_3 +0.00000000e+00 * T_4;
        species[4] = +2.94642878e+00 -1.63816649e-03 * T +2.42103170e-06 * T_2 -1.60284319e-09 * T_3 +3.89069636e-13 * T_4;
        species[5] = +4.12530561e+00 -3.22544939e-03 * T +6.52764691e-06 * T_2 -5.79853643e-09 * T_3 +2.06237379e-12 * T_4;
        species[6] = +4.30179801e+00 -4.74912051e-03 * T +2.11582891e-05 * T_2 -2.42763894e-08 * T_3 +9.29225124e-12 * T_4;
        species[7] = +3.38875365e+00 +6.56922581e-03 * T -1.48501258e-07 * T_2 -4.62580552e-09 * T_3 +2.47151475e-12 * T_4;
        species[8] = +3.29867700e+00 +1.40824000e-03 * T -3.96322200e-06 * T_2 +5.64151500e-09 * T_3 -2.44485500e-12 * T_4;
    } else {
        species[0] = +2.99142337e+00 +7.00064411e-04 * T -5.63382869e-08 * T_2 -9.23157818e-12 * T_3 +1.58275179e-15 * T_4;
        species[1] = +3.69757819e+00 +6.13519689e-04 * T -1.25884199e-07 * T_2 +1.77528148e-11 * T_3 -1.13643531e-15 * T_4;
        species[2] = +2.67214561e+00 +3.05629289e-03 * T -8.73026011e-07 * T_2 +1.20099639e-10 * T_3 -6.39161787e-15 * T_4;
        species[3] = +2.50000000e+00 +0.00000000e+00 * T +0.00000000e+00 * T_2 +0.00000000e+00 * T_3 +0.00000000e+00 * T_4;
        species[4] = +2.54205966e+00 -2.75506191e-05 * T -3.10280335e-09 * T_2 +4.55106742e-12 * T_3 -4.36805150e-16 * T_4;
        species[5] = +2.86472886e+00 +1.05650448e-03 * T -2.59082758e-07 * T_2 +3.05218674e-11 * T_3 -1.33195876e-15 * T_4;
        species[6] = +4.01721090e+00 +2.23982013e-03 * T -6.33658150e-07 * T_2 +1.14246370e-10 * T_3 -1.07908535e-14 * T_4;
        species[7] = +4.57316685e+00 +4.33613639e-03 * T -1.47468882e-06 * T_2 +2.34890357e-10 * T_3 -1.43165356e-14 * T_4;
        species[8] = +2.92664000e+00 +1.48797700e-03 * T -5.68476100e-07 * T_2 +1.00970400e-10 * T_3 -6.75335100e-15 * T_4;
    }
}

void fg_enthalpy_RT(const dfloat log_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
    if (T < 1000.0) {
        species[0] = +3.29812431e+00 +4.12472087e-04 * T -2.71433843e-07 * T_2 -2.36885858e-11 * T_3 +8.26974448e-14 * T_4 -1.01252087e+03 * rcp_T;
        species[1] = +3.21293640e+00 +5.63743175e-04 * T -1.91871682e-07 * T_2 +3.28469308e-10 * T_3 -1.75371078e-13 * T_4 -1.00524902e+03 * rcp_T;
        species[2] = +3.38684249e+00 +1.73749123e-03 * T -2.11823211e-06 * T_2 +1.74214532e-09 * T_3 -5.01317694e-13 * T_4 -3.02081133e+04 * rcp_T;
        species[3] = +2.50000000e+00 +0.00000000e+00 * T +0.00000000e+00 * T_2 +0.00000000e+00 * T_3 +0.00000000e+00 * T_4 +2.54716270e+04 * rcp_T;
        species[4] = +2.94642878e+00 -8.19083245e-04 * T +8.07010567e-07 * T_2 -4.00710797e-10 * T_3 +7.78139272e-14 * T_4 +2.91476445e+04 * rcp_T;
        species[5] = +4.12530561e+00 -1.61272470e-03 * T +2.17588230e-06 * T_2 -1.44963411e-09 * T_3 +4.12474758e-13 * T_4 +3.34630913e+03 * rcp_T;
        species[6] = +4.30179801e+00 -2.37456025e-03 * T +7.05276303e-06 * T_2 -6.06909735e-09 * T_3 +1.85845025e-12 * T_4 +2.94808040e+02 * rcp_T;
        species[7] = +3.38875365e+00 +3.28461290e-03 * T -4.95004193e-08 * T_2 -1.15645138e-09 * T_3 +4.94302950e-13 * T_4 -1.76631465e+04 * rcp_T;
        species[8] = +3.29867700e+00 +7.04120000e-04 * T -1.32107400e-06 * T_2 +1.41037875e-09 * T_3 -4.88971000e-13 * T_4 -1.02090000e+03 * rcp_T;
    } else {
        species[0] = +2.99142337e+00 +3.50032206e-04 * T -1.87794290e-08 * T_2 -2.30789455e-12 * T_3 +3.16550358e-16 * T_4 -8.35033997e+02 * rcp_T;
        species[1] = +3.69757819e+00 +3.06759845e-04 * T -4.19613997e-08 * T_2 +4.43820370e-12 * T_3 -2.27287062e-16 * T_4 -1.23393018e+03 * rcp_T;
        species[2] = +2.67214561e+00 +1.52814644e-03 * T -2.91008670e-07 * T_2 +3.00249098e-11 * T_3 -1.27832357e-15 * T_4 -2.98992090e+04 * rcp_T;
        species[3] = +2.50000000e+00 +0.00000000e+00 * T +0.00000000e+00 * T_2 +0.00000000e+00 * T_3 +0.00000000e+00 * T_4 +2.54716270e+04 * rcp_T;
        species[4] = +2.54205966e+00 -1.37753096e-05 * T -1.03426778e-09 * T_2 +1.13776685e-12 * T_3 -8.73610300e-17 * T_4 +2.92308027e+04 * rcp_T;
        species[5] = +2.86472886e+00 +5.28252240e-04 * T -8.63609193e-08 * T_2 +7.63046685e-12 * T_3 -2.66391752e-16 * T_4 +3.68362875e+03 * rcp_T;
        species[6] = +4.01721090e+00 +1.11991006e-03 * T -2.11219383e-07 * T_2 +2.85615925e-11 * T_3 -2.15817070e-15 * T_4 +1.11856713e+02 * rcp_T;
        species[7] = +4.57316685e+00 +2.16806820e-03 * T -4.91562940e-07 * T_2 +5.87225893e-11 * T_3 -2.86330712e-15 * T_4 -1.80069609e+04 * rcp_T;
        species[8] = +2.92664000e+00 +7.43988500e-04 * T -1.89492033e-07 * T_2 +2.52426000e-11 * T_3 -1.35067020e-15 * T_4 -9.22797700e+02 * rcp_T;
    }
}

void fg_exp_Gibbs_RT(const dfloat log_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
    if (T < 1000.0) {
        species[0] = exp2(-1.46075884e+03 * rcp_T +9.51056079e+00 -3.29812431e+00 * log_T -5.95071434e-04 * T +1.95798130e-07 * T_2 +1.13918018e-11 * T_3 -2.98267984e-14 * T_4);
        species[1] = exp2(-1.45026778e+03 * rcp_T -4.07099858e+00 -3.21293640e+00 * log_T -8.13309483e-04 * T +1.38406162e-07 * T_2 -1.57960347e-10 * T_3 +6.32517463e-14 * T_4);
        species[2] = exp2(-4.35810953e+04 * rcp_T +1.14926478e+00 -3.38684249e+00 * log_T -2.50666998e-03 * T +1.52798148e-06 * T_2 -8.37794803e-10 * T_3 +1.80812138e-13 * T_4);
        species[3] = exp2(+3.67477900e+04 * rcp_T +4.27054699e+00 -2.50000000e+00 * log_T -0.00000000e+00 * T -0.00000000e+00 * T_2 -0.00000000e+00 * T_3 -0.00000000e+00 * T_4);
        species[4] = exp2(+4.20511622e+04 * rcp_T -2.53426696e-02 -2.94642878e+00 * log_T +1.18168734e-03 * T -5.82135071e-07 * T_2 +1.92701160e-10 * T_3 -2.80654417e-14 * T_4);
        species[5] = exp2(+4.82770359e+03 * rcp_T +6.94764215e+00 -4.12530561e+00 * log_T +2.32666992e-03 * T -1.56956730e-06 * T_2 +6.97126646e-10 * T_3 -1.48768822e-13 * T_4);
        species[6] = exp2(+4.25318097e+02 * rcp_T +8.44172171e-01 -4.30179801e+00 * log_T +3.42576630e-03 * T -5.08749313e-06 * T_2 +2.91861888e-09 * T_3 -6.70294239e-13 * T_4);
        species[7] = exp2(-2.54825339e+04 * rcp_T -4.90027175e+00 -3.38875365e+00 * log_T -4.73869475e-03 * T +3.57070047e-08 * T_2 +5.56135557e-10 * T_3 -1.78282104e-13 * T_4);
        species[8] = exp2(-1.47284737e+03 * rcp_T -9.40197145e-01 -3.29867700e+00 * log_T -1.01583043e-03 * T +9.52953454e-07 * T_2 -6.78248809e-10 * T_3 +1.76359009e-13 * T_4);
    } else {
        species[0] = exp2(-1.20469941e+03 * rcp_T +6.27072238e+00 -2.99142337e+00 * log_T -5.04989727e-04 * T +1.35464945e-08 * T_2 +1.10986267e-12 * T_3 -1.14171408e-16 * T_4);
        species[1] = exp2(-1.78018495e+03 * rcp_T +7.33484337e-01 -3.69757819e+00 * log_T -4.42560906e-04 * T +3.02687516e-08 * T_2 -2.13432482e-12 * T_3 +8.19764793e-17 * T_4);
        species[2] = exp2(-4.31354406e+04 * rcp_T -6.04586056e+00 -2.67214561e+00 * log_T -2.20464930e-03 * T +2.09918383e-07 * T_2 -1.44389295e-11 * T_3 +4.61057770e-16 * T_4);
        species[3] = exp2(+3.67477900e+04 * rcp_T +4.27054704e+00 -2.50000000e+00 * log_T -0.00000000e+00 * T -0.00000000e+00 * T_2 -0.00000000e+00 * T_3 -0.00000000e+00 * T_4);
        species[4] = exp2(+4.21711341e+04 * rcp_T -3.43108724e+00 -2.54205966e+00 * log_T +1.98735708e-05 * T +7.46066501e-10 * T_2 -5.47150200e-13 * T_3 +3.15088312e-17 * T_4);
        species[5] = exp2(+5.31435293e+03 * rcp_T -4.09279869e+00 -2.86472886e+00 * log_T -7.62106887e-04 * T +6.22962350e-08 * T_2 -3.66947889e-12 * T_3 +9.60805149e-17 * T_4);
        species[6] = exp2(+1.61375125e+02 * rcp_T +3.34862143e-01 -4.01721090e+00 * log_T -1.61568870e-03 * T +1.52362578e-07 * T_2 -1.37352226e-11 * T_3 +7.78395542e-16 * T_4);
        species[7] = exp2(-2.59785532e+04 * rcp_T +5.87469733e+00 -4.57316685e+00 * log_T -3.12786123e-03 * T +3.54587708e-07 * T_2 -2.82395961e-11 * T_3 +1.03271975e-15 * T_4);
        species[8] = exp2(-1.33131567e+03 * rcp_T -4.40582907e+00 -2.92664000e+00 * log_T -1.07334852e-03 * T +1.36689608e-07 * T_2 -1.21391246e-11 * T_3 +4.87151300e-16 * T_4);
    }
}

dfloat fg_viscosity(dfloat T, const dfloat mole_fractions[]) {
    dfloat T_14 =	sqrt(sqrt(T));
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((3.96729649E-08 + -5.49228488E-07 * ln_T + 3.29052645E-06 * ln_T_2 + -6.10092193E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[1] * pow((1.33781415E-07 + -2.19383368E-06 * ln_T + 1.46168053E-05 * ln_T_2 + -3.18633978E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[2] * pow((-4.00891580E-09 + 1.49878970E-07 * ln_T + -1.23459035E-06 * ln_T_2 + 3.03914265E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[3] * pow((6.74729606E-08 + -1.14565921E-06 * ln_T + 7.88485306E-06 * ln_T_2 + -1.79413191E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[4] * pow((1.41640059E-07 + -2.21121155E-06 * ln_T + 1.42544884E-05 * ln_T_2 + -2.96780713E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[5] * pow((1.46033613E-07 + -2.27980145E-06 * ln_T + 1.46966505E-05 * ln_T_2 + -3.05986599E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[6] * pow((1.35872151E-07 + -2.22811892E-06 * ln_T + 1.48452368E-05 * ln_T_2 + -3.23613591E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[7] * pow((1.37931200E-07 + -2.26188453E-06 * ln_T + 1.50702062E-05 * ln_T_2 + -3.28517733E-05 * ln_T_3)*T_14, 12.)
    + mole_fractions[8] * pow((1.12591716E-07 + -1.81926405E-06 * ln_T + 1.19929798E-05 * ln_T_2 + -2.57678846E-05 * ln_T_3)*T_14, 12.)
    , 1./6.);
}
dfloat fg_thermal_conductivity(dfloat T, const dfloat mole_fractions[]) {
    dfloat T_12 = sqrt(T);
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((1.72659622E-04 + -2.24461746E-03 * ln_T + 9.87640957E-03 * ln_T_2 + -4.39869825E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[1] * pow((2.06263285E-05 + -3.92558038E-04 * ln_T + 3.15126308E-03 * ln_T_2 + -7.51806299E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[2] * pow((-1.33297463E-04 + 3.47918009E-03 * ln_T + -2.70813974E-02 * ln_T_2 + 6.63710275E-02 * ln_T_3)*T_12, 4.)
    + mole_fractions[3] * pow((2.70732991E-04 + -5.96292949E-03 * ln_T + 4.74329801E-02 * ln_T_2 + -1.09806368E-01 * ln_T_3)*T_12, 4.)
    + mole_fractions[4] * pow((2.16721053E-05 + -4.03546198E-04 * ln_T + 2.94852134E-03 * ln_T_2 + -4.78408051E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[5] * pow((-5.92728440E-05 + 1.60790704E-03 * ln_T + -1.24204961E-02 * ln_T_2 + 3.30129237E-02 * ln_T_3)*T_12, 4.)
    + mole_fractions[6] * pow((1.66771525E-06 + 9.94426046E-05 * ln_T + -3.65902661E-04 * ln_T_2 + 2.36053379E-04 * ln_T_3)*T_12, 4.)
    + mole_fractions[7] * pow((-9.59488373E-05 + 2.05216263E-03 * ln_T + -1.28581296E-02 * ln_T_2 + 2.65389544E-02 * ln_T_3)*T_12, 4.)
    + mole_fractions[8] * pow((-7.13350511E-05 + 1.52694323E-03 * ln_T + -1.01926196E-02 * ln_T_2 + 2.32321507E-02 * ln_T_3)*T_12, 4.)
    , 1./4.);
}
void fg_mixture_diffusion_coefficients(const dfloat mole_fractions[n_species], const dfloat mass_fractions[n_species], dfloat T, dfloat* Ddiag) {
    dfloat T_12 = sqrt(T);
    dfloat ln_T = log(T);
    dfloat ln_T_2 = ln_T*ln_T; 
    dfloat ln_T_3 = ln_T_2*ln_T; 
    dfloat T_32 = T*T_12;
    Ddiag[0] = (1. - mass_fractions[0]) * mole_fractions[0] / ( 0.
    + mole_fractions[1] / ((8.62418622E-06 + -1.68915443E-04 * ln_T + 1.40824085E-03 * ln_T_2 + -2.56132655E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((2.17983369E-05 + -4.81332814E-04 * ln_T + 3.95481711E-03 * ln_T_2 + -9.22873381E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.80532870E-05 + -5.65916495E-04 * ln_T + 4.64888361E-03 * ln_T_2 + -9.09176434E-03 * ln_T_3)*T_32)
    + mole_fractions[4] / ((9.37178644E-06 + -1.76843510E-04 * ln_T + 1.50978537E-03 * ln_T_2 + -2.48211304E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((9.34065782E-06 + -1.76256120E-04 * ln_T + 1.50477058E-03 * ln_T_2 + -2.47386864E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((8.61637822E-06 + -1.68762514E-04 * ln_T + 1.40696588E-03 * ln_T_2 + -2.55900763E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((8.60902651E-06 + -1.68618521E-04 * ln_T + 1.40576542E-03 * ln_T_2 + -2.55682422E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((7.75349377E-06 + -1.50122099E-04 * ln_T + 1.26050651E-03 * ln_T_2 + -2.22447501E-03 * ln_T_3)*T_32)
    );
    Ddiag[1] = (1. - mass_fractions[1]) * mole_fractions[1] / ( 0.
    + mole_fractions[0] / ((8.62418622E-06 + -1.68915443E-04 * ln_T + 1.40824085E-03 * ln_T_2 + -2.56132655E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.25389680E-06 + -1.28235545E-04 * ln_T + 1.17612487E-03 * ln_T_2 + -3.05390468E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.81943576E-05 + -6.10831873E-04 * ln_T + 4.97587795E-03 * ln_T_2 + -1.12592633E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.56979620E-06 + -1.16307445E-04 * ln_T + 9.44652682E-04 * ln_T_2 + -1.99463722E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.45865288E-06 + -1.13986571E-04 * ln_T + 9.25802471E-04 * ln_T_2 + -1.95483493E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((4.03880015E-06 + -8.59650818E-05 * ln_T + 6.97879872E-04 * ln_T_2 + -1.52976288E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((4.00923511E-06 + -8.53357957E-05 * ln_T + 6.92771216E-04 * ln_T_2 + -1.51856463E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.89040364E-06 + -8.23029468E-05 * ln_T + 6.68013786E-04 * ln_T_2 + -1.44740787E-03 * ln_T_3)*T_32)
    );
    Ddiag[2] = (1. - mass_fractions[2]) * mole_fractions[2] / ( 0.
    + mole_fractions[0] / ((2.17983369E-05 + -4.81332814E-04 * ln_T + 3.95481711E-03 * ln_T_2 + -9.22873381E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.25389680E-06 + -1.28235545E-04 * ln_T + 1.17612487E-03 * ln_T_2 + -3.05390468E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.01460159E-05 + -5.41017848E-04 * ln_T + 5.55554866E-03 * ln_T_2 + -1.53159246E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((9.59148263E-06 + -2.24308215E-04 * ln_T + 1.95098418E-03 * ln_T_2 + -4.87968419E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((9.43974601E-06 + -2.20759674E-04 * ln_T + 1.92011974E-03 * ln_T_2 + -4.80248791E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.22492715E-06 + -1.27528425E-04 * ln_T + 1.16963903E-03 * ln_T_2 + -3.03706299E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.19752054E-06 + -1.26859493E-04 * ln_T + 1.16350386E-03 * ln_T_2 + -3.02113251E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((5.54983595E-06 + -1.33176211E-04 * ln_T + 1.19565030E-03 * ln_T_2 + -3.06215755E-03 * ln_T_3)*T_32)
    );
    Ddiag[3] = (1. - mass_fractions[3]) * mole_fractions[3] / ( 0.
    + mole_fractions[0] / ((2.80532870E-05 + -5.65916495E-04 * ln_T + 4.64888361E-03 * ln_T_2 + -9.09176434E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((2.81943576E-05 + -6.10831873E-04 * ln_T + 4.97587795E-03 * ln_T_2 + -1.12592633E-02 * ln_T_3)*T_32)
    + mole_fractions[2] / ((2.01460159E-05 + -5.41017848E-04 * ln_T + 5.55554866E-03 * ln_T_2 + -1.53159246E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((3.47673493E-05 + -7.40275463E-04 * ln_T + 6.00991816E-03 * ln_T_2 + -1.31823966E-02 * ln_T_3)*T_32)
    + mole_fractions[5] / ((3.47062348E-05 + -7.38974197E-04 * ln_T + 5.99935385E-03 * ln_T_2 + -1.31592244E-02 * ln_T_3)*T_32)
    + mole_fractions[6] / ((2.81812077E-05 + -6.10546980E-04 * ln_T + 4.97355719E-03 * ln_T_2 + -1.12540119E-02 * ln_T_3)*T_32)
    + mole_fractions[7] / ((2.81688315E-05 + -6.10278850E-04 * ln_T + 4.97137299E-03 * ln_T_2 + -1.12490696E-02 * ln_T_3)*T_32)
    + mole_fractions[8] / ((2.59510402E-05 + -5.59123123E-04 * ln_T + 4.54951116E-03 * ln_T_2 + -1.01939025E-02 * ln_T_3)*T_32)
    );
    Ddiag[4] = (1. - mass_fractions[4]) * mole_fractions[4] / ( 0.
    + mole_fractions[0] / ((9.37178644E-06 + -1.76843510E-04 * ln_T + 1.50978537E-03 * ln_T_2 + -2.48211304E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.56979620E-06 + -1.16307445E-04 * ln_T + 9.44652682E-04 * ln_T_2 + -1.99463722E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((9.59148263E-06 + -2.24308215E-04 * ln_T + 1.95098418E-03 * ln_T_2 + -4.87968419E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((3.47673493E-05 + -7.40275463E-04 * ln_T + 6.00991816E-03 * ln_T_2 + -1.31823966E-02 * ln_T_3)*T_32)
    + mole_fractions[5] / ((7.07160706E-06 + -1.44488755E-04 * ln_T + 1.18043538E-03 * ln_T_2 + -2.37790885E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.54137502E-06 + -1.15713959E-04 * ln_T + 9.39832372E-04 * ln_T_2 + -1.98445912E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.51450455E-06 + -1.15152855E-04 * ln_T + 9.35275068E-04 * ln_T_2 + -1.97483635E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((5.18959361E-06 + -1.07632457E-04 * ln_T + 8.75498870E-04 * ln_T_2 + -1.82232347E-03 * ln_T_3)*T_32)
    );
    Ddiag[5] = (1. - mass_fractions[5]) * mole_fractions[5] / ( 0.
    + mole_fractions[0] / ((9.34065782E-06 + -1.76256120E-04 * ln_T + 1.50477058E-03 * ln_T_2 + -2.47386864E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((5.45865288E-06 + -1.13986571E-04 * ln_T + 9.25802471E-04 * ln_T_2 + -1.95483493E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((9.43974601E-06 + -2.20759674E-04 * ln_T + 1.92011974E-03 * ln_T_2 + -4.80248791E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((3.47062348E-05 + -7.38974197E-04 * ln_T + 5.99935385E-03 * ln_T_2 + -1.31592244E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((7.07160706E-06 + -1.44488755E-04 * ln_T + 1.18043538E-03 * ln_T_2 + -2.37790885E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((5.42964996E-06 + -1.13380937E-04 * ln_T + 9.20883496E-04 * ln_T_2 + -1.94444850E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((5.40222381E-06 + -1.12808229E-04 * ln_T + 9.16231946E-04 * ln_T_2 + -1.93462674E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((5.09077120E-06 + -1.05582874E-04 * ln_T + 8.58827255E-04 * ln_T_2 + -1.78762202E-03 * ln_T_3)*T_32)
    );
    Ddiag[6] = (1. - mass_fractions[6]) * mole_fractions[6] / ( 0.
    + mole_fractions[0] / ((8.61637822E-06 + -1.68762514E-04 * ln_T + 1.40696588E-03 * ln_T_2 + -2.55900763E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((4.03880015E-06 + -8.59650818E-05 * ln_T + 6.97879872E-04 * ln_T_2 + -1.52976288E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.22492715E-06 + -1.27528425E-04 * ln_T + 1.16963903E-03 * ln_T_2 + -3.03706299E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.81812077E-05 + -6.10546980E-04 * ln_T + 4.97355719E-03 * ln_T_2 + -1.12540119E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.54137502E-06 + -1.15713959E-04 * ln_T + 9.39832372E-04 * ln_T_2 + -1.98445912E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.42964996E-06 + -1.13380937E-04 * ln_T + 9.20883496E-04 * ln_T_2 + -1.94444850E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((3.97756651E-06 + -8.46617357E-05 * ln_T + 6.87299076E-04 * ln_T_2 + -1.50656962E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.86257503E-06 + -8.17142220E-05 * ln_T + 6.63235388E-04 * ln_T_2 + -1.43705435E-03 * ln_T_3)*T_32)
    );
    Ddiag[7] = (1. - mass_fractions[7]) * mole_fractions[7] / ( 0.
    + mole_fractions[0] / ((8.60902651E-06 + -1.68618521E-04 * ln_T + 1.40576542E-03 * ln_T_2 + -2.55682422E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((4.00923511E-06 + -8.53357957E-05 * ln_T + 6.92771216E-04 * ln_T_2 + -1.51856463E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.19752054E-06 + -1.26859493E-04 * ln_T + 1.16350386E-03 * ln_T_2 + -3.02113251E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.81688315E-05 + -6.10278850E-04 * ln_T + 4.97137299E-03 * ln_T_2 + -1.12490696E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.51450455E-06 + -1.15152855E-04 * ln_T + 9.35275068E-04 * ln_T_2 + -1.97483635E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.40222381E-06 + -1.12808229E-04 * ln_T + 9.16231946E-04 * ln_T_2 + -1.93462674E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((3.97756651E-06 + -8.46617357E-05 * ln_T + 6.87299076E-04 * ln_T_2 + -1.50656962E-03 * ln_T_3)*T_32)
    + mole_fractions[8] / ((3.83621145E-06 + -8.11564905E-05 * ln_T + 6.58708547E-04 * ln_T_2 + -1.42724590E-03 * ln_T_3)*T_32)
    );
    Ddiag[8] = (1. - mass_fractions[8]) * mole_fractions[8] / ( 0.
    + mole_fractions[0] / ((7.75349377E-06 + -1.50122099E-04 * ln_T + 1.26050651E-03 * ln_T_2 + -2.22447501E-03 * ln_T_3)*T_32)
    + mole_fractions[1] / ((3.89040364E-06 + -8.23029468E-05 * ln_T + 6.68013786E-04 * ln_T_2 + -1.44740787E-03 * ln_T_3)*T_32)
    + mole_fractions[2] / ((5.54983595E-06 + -1.33176211E-04 * ln_T + 1.19565030E-03 * ln_T_2 + -3.06215755E-03 * ln_T_3)*T_32)
    + mole_fractions[3] / ((2.59510402E-05 + -5.59123123E-04 * ln_T + 4.54951116E-03 * ln_T_2 + -1.01939025E-02 * ln_T_3)*T_32)
    + mole_fractions[4] / ((5.18959361E-06 + -1.07632457E-04 * ln_T + 8.75498870E-04 * ln_T_2 + -1.82232347E-03 * ln_T_3)*T_32)
    + mole_fractions[5] / ((5.09077120E-06 + -1.05582874E-04 * ln_T + 8.58827255E-04 * ln_T_2 + -1.78762202E-03 * ln_T_3)*T_32)
    + mole_fractions[6] / ((3.86257503E-06 + -8.17142220E-05 * ln_T + 6.63235388E-04 * ln_T_2 + -1.43705435E-03 * ln_T_3)*T_32)
    + mole_fractions[7] / ((3.83621145E-06 + -8.11564905E-05 * ln_T + 6.58708547E-04 * ln_T_2 + -1.42724590E-03 * ln_T_3)*T_32)
    );
}
void fg_rates(const dfloat log_T, const dfloat T, const dfloat T2, const dfloat T4, const dfloat rcp_T, const dfloat rcp_T2, const dfloat P0_RT, const dfloat rcp_P0_RT,const dfloat exp_Gibbs0_RT[], const dfloat concentrations[], dfloat* molar_rates) {
    dfloat c, Pr, logFcent, logPr_c, f1;
    c = exp2(-12050.7 * rcp_T + -0.406 * log_T + 31.724);
    const dfloat cR0 = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[4]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3]) * concentrations[4]*concentrations[5]);
    c = exp2(-4566.49 * rcp_T + 2.67 * log_T + -4.29903);
    const dfloat cR1 = c * (concentrations[0]*concentrations[4] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[4]) * concentrations[3]*concentrations[5]);
    c = exp2(-2490.15 * rcp_T + 1.51 * log_T + 7.75489);
    const dfloat cR2 = c * (concentrations[0]*concentrations[5] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[5]) * concentrations[2]*concentrations[3]);
    c = exp2(-9728.3 * rcp_T + 2.02 * log_T + 1.57046);
    const dfloat cR3 = c * (concentrations[2]*concentrations[4] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[4]) * concentrations[5]*concentrations[5]);
    c = exp2(-75779.1 * rcp_T + -1.4 * log_T + 45.3795) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat cR4 = c * (concentrations[0] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0])* rcp_P0_RT * concentrations[3]*concentrations[3]);
    c = exp2(-0 * rcp_T + -0.5 * log_T + 12.5899) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat cR5 = c * (concentrations[4]*concentrations[4] - exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[4])* P0_RT * concentrations[1]);
    c = exp2(-0 * rcp_T + -1 * log_T + 22.1685) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat cR6 = c * (concentrations[3]*concentrations[4] - exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[4])* P0_RT * concentrations[5]);
    c = exp2(-0 * rcp_T + -2 * log_T + 35.1453) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat cR7 = c * (concentrations[3]*concentrations[5] - exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5])* P0_RT * concentrations[2]);
    Pr = exp2(-381.001 * rcp_T + -1.72 * log_T + 29.2458) * (2*concentrations[0] + 0.78*concentrations[1] + 11*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    logFcent = log2(0.2 * exp2(-1.4427e+30*T) + 0.8 * exp2(-1.4427e-30*T) + exp2(0*rcp_T));
    logPr_c = log2(Pr) - 0.67*logFcent - 1.32877;
    f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.49145);
    c = Pr / (exp2(0 * rcp_T - 0.6 * log_T - 20.4923) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const dfloat cR8 = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3])* P0_RT * concentrations[6]);
    c = exp2(-597.492 * rcp_T + 0 * log_T + 23.9847);
    const dfloat cR9 = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[0]*concentrations[1]);
    c = exp2(-214.168 * rcp_T + 0 * log_T + 26.077);
    const dfloat cR10 = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[5]*concentrations[5]);
    c = exp2(-0 * rcp_T + 0 * log_T + 24.9539);
    const dfloat cR11 = c * (concentrations[4]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[5]);
    c = exp2(360.818 * rcp_T + 0 * log_T + 24.7846);
    const dfloat cR12 = c * (concentrations[5]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[2]);
    c = exp2(-8698.84 * rcp_T + 0 * log_T + 28.6458);
    const dfloat cR13 = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);
    c = exp2(1182.86 * rcp_T + 0 * log_T + 16.9882);
    const dfloat cR14 = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);
    Pr = exp2(-33032.7 * rcp_T + 0 * log_T + 36.8066) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    logFcent = log2(0.5 * exp2(-1.4427e+30*T) + 0.5 * exp2(-1.4427e-30*T) + exp2(0*rcp_T));
    logPr_c = log2(Pr) - 0.67*logFcent - 1.32877;
    f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.49145);
    c = Pr / (exp2(35159.8 * rcp_T - 0 * log_T - 48.0682) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const dfloat cR15 = c * (concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[7])* rcp_P0_RT * concentrations[5]*concentrations[5]);
    c = exp2(-2882.19 * rcp_T + 0 * log_T + 24.5225);
    const dfloat cR16 = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[5]);
    c = exp2(-5771.64 * rcp_T + 0 * log_T + 25.5225);
    const dfloat cR17 = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[0]*concentrations[6]);
    c = exp2(-2882.19 * rcp_T + 2 * log_T + 3.2555);
    const dfloat cR18 = c * (concentrations[4]*concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[7]) * concentrations[5]*concentrations[6]);
    c = exp2(-0 * rcp_T + 0 * log_T + 19.9316);
    const dfloat cR19 = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);
    c = exp2(-6938.31 * rcp_T + 0 * log_T + 29.1115);
    const dfloat cR20 = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);
    molar_rates[0] = -cR1+-cR2+-cR4+cR9+cR17;
    molar_rates[1] = -cR0+cR5+-cR8+cR9+cR11+cR12+cR13+cR14;
    molar_rates[2] = cR2+-cR3+cR7+cR12+cR16+cR19+cR20;
    molar_rates[3] = -cR0+cR1+cR2+2*cR4+-cR6+-cR7+-cR8+-cR9+-cR10+-cR16+-cR17;
    molar_rates[4] = cR0+-cR1+-cR3+-2*cR5+-cR6+-cR11+-cR18;
    molar_rates[5] = cR0+cR1+-cR2+2*cR3+cR6+-cR7+2*cR10+cR11+-cR12+2*cR15+cR16+cR18+-cR19+-cR20;
    molar_rates[6] = cR8+-cR9+-cR10+-cR11+-cR12+-2*cR13+-2*cR14+cR17+cR18+cR19+cR20;
    molar_rates[7] = cR13+cR14+-cR15+-cR16+-cR17+-cR18+-cR19+-cR20;
}

