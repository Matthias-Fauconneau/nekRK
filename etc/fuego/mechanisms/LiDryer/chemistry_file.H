#define n_species 9
const float fg_molar_mass[9] = {0.00201594, 0.0319988, 0.01801534, 0.00100797, 0.0159994, 0.01700737, 0.03300677, 0.03401474, 0.0280134};
const float fg_rcp_molar_mass[9] = {496.046509321, 31.2511719189, 55.5082501912, 992.093018641, 62.5023438379, 58.7980387326, 30.2968148656, 29.3990193663, 35.6972020533};
void fg_molar_heat_capacity_at_constant_pressure_R(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {
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

void fg_enthalpy_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {
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

void fg_exp_Gibbs_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* species) {
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

float fg_viscosity(float T, const float mole_fractions[]) {
    float T_14 =	sqrt(sqrt(T));
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T; 
    float ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((1.49596873E-08 + -2.42222120E-08 * ln_T + -3.89210451E-07 * ln_T_2 + 2.39433308E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[1] * pow((3.05929516E-08 + -1.32625707E-08 * ln_T + -5.96041239E-07 * ln_T_2 + 3.11685631E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[2] * pow((3.41495491E-09 + -8.61954936E-09 * ln_T + -1.17453786E-07 * ln_T_2 + 4.43745536E-07 * ln_T_3)*T_14, 12.)
    + mole_fractions[3] * pow((1.29259894E-08 + 8.35264197E-09 * ln_T + -1.75669910E-07 * ln_T_2 + 6.15726371E-07 * ln_T_3)*T_14, 12.)
    + mole_fractions[4] * pow((3.88414490E-08 + -4.00784798E-08 * ln_T + -8.84162693E-07 * ln_T_2 + 5.11203081E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[5] * pow((4.00462777E-08 + -4.13216802E-08 * ln_T + -9.11588668E-07 * ln_T_2 + 5.27060166E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[6] * pow((3.10710584E-08 + -1.34698382E-08 * ln_T + -6.05356174E-07 * ln_T_2 + 3.16556656E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[7] * pow((3.15419189E-08 + -1.36739644E-08 * ln_T + -6.14529931E-07 * ln_T_2 + 3.21353855E-06 * ln_T_3)*T_14, 12.)
    + mole_fractions[8] * pow((2.73517122E-08 + -1.83558251E-08 * ln_T + -5.68509665E-07 * ln_T_2 + 3.10977871E-06 * ln_T_3)*T_14, 12.)
    , 1./6.);
}
float fg_thermal_conductivity(float T, const float mole_fractions[]) {
    float T_12 = sqrt(T);
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T; 
    float ln_T_3 = ln_T_2*ln_T; 
    return pow(0.
    + mole_fractions[0] * pow((7.49956429E-05 + -3.00587508E-04 * ln_T + -2.84357238E-03 * ln_T_2 + 2.28808225E-02 * ln_T_3)*T_12, 4.)
    + mole_fractions[1] * pow((1.59087902E-06 + 2.43351199E-05 * ln_T + 1.39419979E-04 * ln_T_2 + -3.51346544E-04 * ln_T_3)*T_12, 4.)
    + mole_fractions[2] * pow((3.43677584E-05 + -1.30340125E-04 * ln_T + -1.43534574E-03 * ln_T_2 + 6.31780241E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[3] * pow((-3.04017766E-05 + 3.86113384E-04 * ln_T + 3.23817508E-03 * ln_T_2 + -8.41359045E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[4] * pow((2.44628475E-06 + 5.60649971E-06 * ln_T + 7.68802066E-05 * ln_T_2 + 1.85200325E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[5] * pow((1.73570939E-05 + -5.28390823E-05 * ln_T + -5.41567611E-04 * ln_T_2 + 5.01011835E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[6] * pow((4.80353533E-06 + 3.16552121E-05 * ln_T + 1.07614738E-04 * ln_T_2 + -8.30265475E-04 * ln_T_3)*T_12, 4.)
    + mole_fractions[7] * pow((-4.98221992E-06 + 1.14542732E-04 * ln_T + 7.62505945E-04 * ln_T_2 + -5.00860724E-03 * ln_T_3)*T_12, 4.)
    + mole_fractions[8] * pow((-8.14370910E-07 + 3.42507739E-05 * ln_T + 2.44127293E-04 * ln_T_2 + -8.35094654E-04 * ln_T_3)*T_12, 4.)
    , 1./4.);
}
void fg_mixture_diffusion_coefficients(const float mole_fractions[n_species], const float mass_fractions[n_species], float T, float* Ddiag) {
    float T_12 = sqrt(T);
    float ln_T = log(T);
    float ln_T_2 = ln_T*ln_T; 
    float ln_T_3 = ln_T_2*ln_T; 
    float T_32 = T*T_12;
    Ddiag[0] = (1. - mass_fractions[0]) * mole_fractions[0] / ( 0.
    + mole_fractions[1] / ((-1.09818746E-07 + 1.47734688E-05 * ln_T + 1.32898503E-04 * ln_T_2 + 3.56617200E-04 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-2.88138478E-06 + 3.95373684E-05 * ln_T + 3.25344837E-04 * ln_T_2 + -8.92935299E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-9.99890802E-07 + 4.54635129E-05 * ln_T + 4.01555958E-04 * ln_T_2 + 6.32215802E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((1.11702259E-07 + 1.76342979E-05 * ln_T + 1.61454219E-04 * ln_T_2 + 5.98464666E-04 * ln_T_3)*T_32)
    + mole_fractions[5] / ((1.11331237E-07 + 1.75757252E-05 * ln_T + 1.60917945E-04 * ln_T_2 + 5.96476851E-04 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-1.09719321E-07 + 1.47600935E-05 * ln_T + 1.32778182E-04 * ln_T_2 + 3.56294333E-04 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-1.09625705E-07 + 1.47474998E-05 * ln_T + 1.32664893E-04 * ln_T_2 + 3.55990334E-04 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.70502511E-08 + 1.36690465E-05 * ln_T + 1.23700820E-04 * ln_T_2 + 3.75628034E-04 * ln_T_3)*T_32)
    );
    Ddiag[1] = (1. - mass_fractions[1]) * mole_fractions[1] / ( 0.
    + mole_fractions[0] / ((-1.09818746E-07 + 1.47734688E-05 * ln_T + 1.32898503E-04 * ln_T_2 + 3.56617200E-04 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-1.65374027E-06 + 1.84981540E-05 * ln_T + 1.46907301E-04 * ln_T_2 + -6.73973893E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-3.01463561E-06 + 4.72562396E-05 * ln_T + 3.94490818E-04 * ln_T_2 + -7.47631897E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.68979399E-07 + 8.74960591E-06 * ln_T + 7.52796157E-05 * ln_T_2 + -2.88223471E-06 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.61616546E-07 + 8.57501061E-06 * ln_T + 7.37774375E-05 * ln_T_2 + -2.82472073E-06 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.49417298E-07 + 6.48258589E-06 * ln_T + 5.49042830E-05 * ln_T_2 + -5.59685427E-05 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.46859475E-07 + 6.43513172E-06 * ln_T + 5.45023698E-05 * ln_T_2 + -5.55588389E-05 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.11059913E-07 + 6.18998909E-06 * ln_T + 5.26896356E-05 * ln_T_2 + -3.73432366E-05 * ln_T_3)*T_32)
    );
    Ddiag[2] = (1. - mass_fractions[2]) * mole_fractions[2] / ( 0.
    + mole_fractions[0] / ((-2.88138478E-06 + 3.95373684E-05 * ln_T + 3.25344837E-04 * ln_T_2 + -8.92935299E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-1.65374027E-06 + 1.84981540E-05 * ln_T + 1.46907301E-04 * ln_T_2 + -6.73973893E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-1.06561921E-05 + 1.16677757E-04 * ln_T + 9.18339329E-04 * ln_T_2 + -4.53594959E-03 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-2.20108308E-06 + 2.54691643E-05 * ln_T + 2.04129695E-04 * ln_T_2 + -8.52442397E-04 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-2.16626210E-06 + 2.50662438E-05 * ln_T + 2.00900376E-04 * ln_T_2 + -8.38956814E-04 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-1.64461851E-06 + 1.83961235E-05 * ln_T + 1.46097010E-04 * ln_T_2 + -6.70256209E-04 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-1.63599190E-06 + 1.82996293E-05 * ln_T + 1.45330679E-04 * ln_T_2 + -6.66740476E-04 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-1.55449010E-06 + 1.75474779E-05 * ln_T + 1.39772832E-04 * ln_T_2 + -6.23695135E-04 * ln_T_3)*T_32)
    );
    Ddiag[3] = (1. - mass_fractions[3]) * mole_fractions[3] / ( 0.
    + mole_fractions[0] / ((-9.99890802E-07 + 4.54635129E-05 * ln_T + 4.01555958E-04 * ln_T_2 + 6.32215802E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.01463561E-06 + 4.72562396E-05 * ln_T + 3.94490818E-04 * ln_T_2 + -7.47631897E-04 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-1.06561921E-05 + 1.16677757E-04 * ln_T + 9.18339329E-04 * ln_T_2 + -4.53594959E-03 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.02125816E-06 + 5.58409267E-05 * ln_T + 4.72808958E-04 * ln_T_2 + -4.90294721E-04 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.01594735E-06 + 5.57427688E-05 * ln_T + 4.71977848E-04 * ln_T_2 + -4.89432874E-04 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.01322958E-06 + 4.72341992E-05 * ln_T + 3.94306826E-04 * ln_T_2 + -7.47283200E-04 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.01190628E-06 + 4.72134556E-05 * ln_T + 3.94133661E-04 * ln_T_2 + -7.46955020E-04 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-2.60459823E-06 + 4.28804215E-05 * ln_T + 3.59584237E-04 * ln_T_2 + -5.83011389E-04 * ln_T_3)*T_32)
    );
    Ddiag[4] = (1. - mass_fractions[4]) * mole_fractions[4] / ( 0.
    + mole_fractions[0] / ((1.11702259E-07 + 1.76342979E-05 * ln_T + 1.61454219E-04 * ln_T_2 + 5.98464666E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.68979399E-07 + 8.74960591E-06 * ln_T + 7.52796157E-05 * ln_T_2 + -2.88223471E-06 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-2.20108308E-06 + 2.54691643E-05 * ln_T + 2.04129695E-04 * ln_T_2 + -8.52442397E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-3.02125816E-06 + 5.58409267E-05 * ln_T + 4.72808958E-04 * ln_T_2 + -4.90294721E-04 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.29950377E-07 + 1.12557539E-05 * ln_T + 9.85126080E-05 * ln_T_2 + 9.90554233E-05 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.67096596E-07 + 8.70495901E-06 * ln_T + 7.48954839E-05 * ln_T_2 + -2.86752744E-06 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.65316523E-07 + 8.66274814E-06 * ln_T + 7.45323112E-05 * ln_T_2 + -2.85362262E-06 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.11727535E-07 + 8.17342018E-06 * ln_T + 7.07095502E-05 * ln_T_2 + 2.08591165E-05 * ln_T_3)*T_32)
    );
    Ddiag[5] = (1. - mass_fractions[5]) * mole_fractions[5] / ( 0.
    + mole_fractions[0] / ((1.11331237E-07 + 1.75757252E-05 * ln_T + 1.60917945E-04 * ln_T_2 + 5.96476851E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.61616546E-07 + 8.57501061E-06 * ln_T + 7.37774375E-05 * ln_T_2 + -2.82472073E-06 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-2.16626210E-06 + 2.50662438E-05 * ln_T + 2.00900376E-04 * ln_T_2 + -8.38956814E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-3.01594735E-06 + 5.57427688E-05 * ln_T + 4.71977848E-04 * ln_T_2 + -4.89432874E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.29950377E-07 + 1.12557539E-05 * ln_T + 9.85126080E-05 * ln_T_2 + 9.90554233E-05 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.59695204E-07 + 8.52944986E-06 * ln_T + 7.33854431E-05 * ln_T_2 + -2.80971242E-06 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.57878318E-07 + 8.48636605E-06 * ln_T + 7.30147599E-05 * ln_T_2 + -2.79552006E-06 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.05791489E-07 + 8.01777850E-06 * ln_T + 6.93630694E-05 * ln_T_2 + 2.04619084E-05 * ln_T_3)*T_32)
    );
    Ddiag[6] = (1. - mass_fractions[6]) * mole_fractions[6] / ( 0.
    + mole_fractions[0] / ((-1.09719321E-07 + 1.47600935E-05 * ln_T + 1.32778182E-04 * ln_T_2 + 3.56294333E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.49417298E-07 + 6.48258589E-06 * ln_T + 5.49042830E-05 * ln_T_2 + -5.59685427E-05 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-1.64461851E-06 + 1.83961235E-05 * ln_T + 1.46097010E-04 * ln_T_2 + -6.70256209E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-3.01322958E-06 + 4.72341992E-05 * ln_T + 3.94306826E-04 * ln_T_2 + -7.47283200E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.67096596E-07 + 8.70495901E-06 * ln_T + 7.48954839E-05 * ln_T_2 + -2.86752744E-06 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.59695204E-07 + 8.52944986E-06 * ln_T + 7.33854431E-05 * ln_T_2 + -2.80971242E-06 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.44119662E-07 + 6.38430117E-06 * ln_T + 5.40718603E-05 * ln_T_2 + -5.51199844E-05 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.08834857E-07 + 6.14571121E-06 * ln_T + 5.23127391E-05 * ln_T_2 + -3.70761150E-05 * ln_T_3)*T_32)
    );
    Ddiag[7] = (1. - mass_fractions[7]) * mole_fractions[7] / ( 0.
    + mole_fractions[0] / ((-1.09625705E-07 + 1.47474998E-05 * ln_T + 1.32664893E-04 * ln_T_2 + 3.55990334E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.46859475E-07 + 6.43513172E-06 * ln_T + 5.45023698E-05 * ln_T_2 + -5.55588389E-05 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-1.63599190E-06 + 1.82996293E-05 * ln_T + 1.45330679E-04 * ln_T_2 + -6.66740476E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-3.01190628E-06 + 4.72134556E-05 * ln_T + 3.94133661E-04 * ln_T_2 + -7.46955020E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.65316523E-07 + 8.66274814E-06 * ln_T + 7.45323112E-05 * ln_T_2 + -2.85362262E-06 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.57878318E-07 + 8.48636605E-06 * ln_T + 7.30147599E-05 * ln_T_2 + -2.79552006E-06 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.44119662E-07 + 6.38430117E-06 * ln_T + 5.40718603E-05 * ln_T_2 + -5.51199844E-05 * ln_T_3)*T_32)
    + mole_fractions[8] / ((-3.06726939E-07 + 6.10376433E-06 * ln_T + 5.19556842E-05 * ln_T_2 + -3.68230561E-05 * ln_T_3)*T_32)
    );
    Ddiag[8] = (1. - mass_fractions[8]) * mole_fractions[8] / ( 0.
    + mole_fractions[0] / ((-3.70502511E-08 + 1.36690465E-05 * ln_T + 1.23700820E-04 * ln_T_2 + 3.75628034E-04 * ln_T_3)*T_32)
    + mole_fractions[1] / ((-3.11059913E-07 + 6.18998909E-06 * ln_T + 5.26896356E-05 * ln_T_2 + -3.73432366E-05 * ln_T_3)*T_32)
    + mole_fractions[2] / ((-1.55449010E-06 + 1.75474779E-05 * ln_T + 1.39772832E-04 * ln_T_2 + -6.23695135E-04 * ln_T_3)*T_32)
    + mole_fractions[3] / ((-2.60459823E-06 + 4.28804215E-05 * ln_T + 3.59584237E-04 * ln_T_2 + -5.83011389E-04 * ln_T_3)*T_32)
    + mole_fractions[4] / ((-3.11727535E-07 + 8.17342018E-06 * ln_T + 7.07095502E-05 * ln_T_2 + 2.08591165E-05 * ln_T_3)*T_32)
    + mole_fractions[5] / ((-3.05791489E-07 + 8.01777850E-06 * ln_T + 6.93630694E-05 * ln_T_2 + 2.04619084E-05 * ln_T_3)*T_32)
    + mole_fractions[6] / ((-3.08834857E-07 + 6.14571121E-06 * ln_T + 5.23127391E-05 * ln_T_2 + -3.70761150E-05 * ln_T_3)*T_32)
    + mole_fractions[7] / ((-3.06726939E-07 + 6.10376433E-06 * ln_T + 5.19556842E-05 * ln_T_2 + -3.68230561E-05 * ln_T_3)*T_32)
    );
}
void fg_rates(const float log_T, const float T, const float T2, const float T4, const float rcp_T, const float rcp_T2, const float P0_RT, const float rcp_P0_RT,const float exp_Gibbs0_RT[], const float concentrations[], float* molar_rates) {
    float c, Pr, logFcent, logPr_c, f1;
    c = exp2(-12050.7 * rcp_T + -0.406 * log_T + 31.724);
    const float cR0 = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[4]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3]) * concentrations[4]*concentrations[5]);
    c = exp2(-4566.49 * rcp_T + 2.67 * log_T + -4.29903);
    const float cR1 = c * (concentrations[0]*concentrations[4] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[4]) * concentrations[3]*concentrations[5]);
    c = exp2(-2490.15 * rcp_T + 1.51 * log_T + 7.75489);
    const float cR2 = c * (concentrations[0]*concentrations[5] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[5]) * concentrations[2]*concentrations[3]);
    c = exp2(-9728.3 * rcp_T + 2.02 * log_T + 1.57046);
    const float cR3 = c * (concentrations[2]*concentrations[4] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[4]) * concentrations[5]*concentrations[5]);
    c = exp2(-75779.1 * rcp_T + -1.4 * log_T + 45.3795) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const float cR4 = c * (concentrations[0] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0])* rcp_P0_RT * concentrations[3]*concentrations[3]);
    c = exp2(-0 * rcp_T + -0.5 * log_T + 12.5899) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const float cR5 = c * (concentrations[4]*concentrations[4] - exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[4])* P0_RT * concentrations[1]);
    c = exp2(-0 * rcp_T + -1 * log_T + 22.1685) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const float cR6 = c * (concentrations[3]*concentrations[4] - exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[4])* P0_RT * concentrations[5]);
    c = exp2(-0 * rcp_T + -2 * log_T + 35.1453) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const float cR7 = c * (concentrations[3]*concentrations[5] - exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5])* P0_RT * concentrations[2]);
    Pr = exp2(-381.001 * rcp_T + -1.72 * log_T + 29.2458) * (2*concentrations[0] + 0.78*concentrations[1] + 11*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    logFcent = log2(0.2 * exp2(-1.4427e+30*T) + 0.8 * exp2(-1.4427e-30*T) + exp2(0*rcp_T));
    logPr_c = log2(Pr) - 0.67*logFcent - 1.32877;
    f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.49145);
    c = Pr / (exp2(0 * rcp_T - 0.6 * log_T - 20.4923) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const float cR8 = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3])* P0_RT * concentrations[6]);
    c = exp2(-597.492 * rcp_T + 0 * log_T + 23.9847);
    const float cR9 = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[0]*concentrations[1]);
    c = exp2(-214.168 * rcp_T + 0 * log_T + 26.077);
    const float cR10 = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[5]*concentrations[5]);
    c = exp2(-0 * rcp_T + 0 * log_T + 24.9539);
    const float cR11 = c * (concentrations[4]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[5]);
    c = exp2(360.818 * rcp_T + 0 * log_T + 24.7846);
    const float cR12 = c * (concentrations[5]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[2]);
    c = exp2(-8698.84 * rcp_T + 0 * log_T + 28.6458);
    const float cR13 = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);
    c = exp2(1182.86 * rcp_T + 0 * log_T + 16.9882);
    const float cR14 = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);
    Pr = exp2(-33032.7 * rcp_T + 0 * log_T + 36.8066) * (2.5*concentrations[0] + concentrations[1] + 12*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    logFcent = log2(0.5 * exp2(-1.4427e+30*T) + 0.5 * exp2(-1.4427e-30*T) + exp2(0*rcp_T));
    logPr_c = log2(Pr) - 0.67*logFcent - 1.32877;
    f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.49145);
    c = Pr / (exp2(35159.8 * rcp_T - 0 * log_T - 48.0682) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const float cR15 = c * (concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[7])* rcp_P0_RT * concentrations[5]*concentrations[5]);
    c = exp2(-2882.19 * rcp_T + 0 * log_T + 24.5225);
    const float cR16 = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[5]);
    c = exp2(-5771.64 * rcp_T + 0 * log_T + 25.5225);
    const float cR17 = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[0]*concentrations[6]);
    c = exp2(-2882.19 * rcp_T + 2 * log_T + 3.2555);
    const float cR18 = c * (concentrations[4]*concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[7]) * concentrations[5]*concentrations[6]);
    c = exp2(-0 * rcp_T + 0 * log_T + 19.9316);
    const float cR19 = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);
    c = exp2(-6938.31 * rcp_T + 0 * log_T + 29.1115);
    const float cR20 = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);
    molar_rates[0] = -cR1+-cR2+-cR4+cR9+cR17;
    molar_rates[1] = -cR0+cR5+-cR8+cR9+cR11+cR12+cR13+cR14;
    molar_rates[2] = cR2+-cR3+cR7+cR12+cR16+cR19+cR20;
    molar_rates[3] = -cR0+cR1+cR2+2*cR4+-cR6+-cR7+-cR8+-cR9+-cR10+-cR16+-cR17;
    molar_rates[4] = cR0+-cR1+-cR3+-2*cR5+-cR6+-cR11+-cR18;
    molar_rates[5] = cR0+cR1+-cR2+2*cR3+cR6+-cR7+2*cR10+cR11+-cR12+2*cR15+cR16+cR18+-cR19+-cR20;
    molar_rates[6] = cR8+-cR9+-cR10+-cR11+-cR12+-2*cR13+-2*cR14+cR17+cR18+cR19+cR20;
    molar_rates[7] = cR13+cR14+-cR15+-cR16+-cR17+-cR18+-cR19+-cR20;
}

