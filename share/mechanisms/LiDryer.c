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

void fg_rates(const dfloat log_T, const dfloat T, const dfloat T2, const dfloat T4, const dfloat rcp_T, const dfloat rcp_T2, const dfloat P0_RT, const dfloat rcp_P0_RT, const dfloat exp_Gibbs0_RT[], const dfloat concentrations[], dfloat* molar_rates) {
    dfloat cR[21];
    {const dfloat c = exp2(-1.205075e+04 * rcp_T -4.060000e-01 * log_T + 3.172395e+01);
    cR[0] = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[4]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3]) * concentrations[4]*concentrations[5]);}
    {const dfloat c = exp2(-4.566492e+03 * rcp_T + 2.670000e+00 * log_T -4.299028e+00);
    cR[1] = c * (concentrations[0]*concentrations[4] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[4]) * concentrations[3]*concentrations[5]);}
    {const dfloat c = exp2(-2.490154e+03 * rcp_T + 1.510000e+00 * log_T + 7.754888e+00);
    cR[2] = c * (concentrations[0]*concentrations[5] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[5]) * concentrations[2]*concentrations[3]);}
    {const dfloat c = exp2(-9.728297e+03 * rcp_T + 2.020000e+00 * log_T + 1.570463e+00);
    cR[3] = c * (concentrations[2]*concentrations[4] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[4]) * concentrations[5]*concentrations[5]);}
    {const dfloat c = exp2(-7.577908e+04 * rcp_T -1.400000e+00 * log_T + 4.537947e+01) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    cR[4] = c * (concentrations[0] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[0])* rcp_P0_RT * concentrations[3]*concentrations[3]);}
    {const dfloat c = exp2(-0.000000e+00 * rcp_T -5.000000e-01 * log_T + 1.258989e+01) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    cR[5] = c * (concentrations[4]*concentrations[4] - exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[4])* P0_RT * concentrations[1]);}
    {const dfloat c = exp2(-0.000000e+00 * rcp_T -1.000000e+00 * log_T + 2.216852e+01) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    cR[6] = c * (concentrations[3]*concentrations[4] - exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[4])* P0_RT * concentrations[5]);}
    {const dfloat c = exp2(-0.000000e+00 * rcp_T -2.000000e+00 * log_T + 3.514528e+01) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    cR[7] = c * (concentrations[3]*concentrations[5] - exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5])* P0_RT * concentrations[2]);}
    {const dfloat Pr = exp2(-3.810008e+02 * rcp_T -1.720000e+00 * log_T + 2.924581e+01) * (2.000000e+00*concentrations[0] + 7.800000e-01*concentrations[1] + 1.100000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat logFcent = log2(2.000000e-01 * exp2(-1.442695e+30*T) + 8.000000e-01 * exp2(-1.442695e-30*T) + exp2(0.000000e+00*rcp_T));
    const dfloat logPr_c = log2(Pr) - 0.67*logFcent - 1.328771e+00;
    const dfloat f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.491446e+00);
    const dfloat c = Pr / (exp2(0.000000e+00 * rcp_T - 6.000000e-01 * log_T - 2.049228e+01) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    cR[8] = c * (concentrations[1]*concentrations[3] - exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3])* P0_RT * concentrations[6]);}
    {const dfloat c = exp2(-5.974917e+02 * rcp_T + 0.000000e+00 * log_T + 2.398468e+01);
    cR[9] = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[1]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[0]*concentrations[1]);}
    {const dfloat c = exp2(-2.141677e+02 * rcp_T + 0.000000e+00 * log_T + 2.607704e+01);
    cR[10] = c * (concentrations[3]*concentrations[6] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]) * concentrations[5]*concentrations[5]);}
    {const dfloat c = exp2(-0.000000e+00 * rcp_T + 0.000000e+00 * log_T + 2.495394e+01);
    cR[11] = c * (concentrations[4]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[5]);}
    {const dfloat c = exp2(3.608182e+02 * rcp_T + 0.000000e+00 * log_T + 2.478457e+01);
    cR[12] = c * (concentrations[5]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[2]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[2]);}
    {const dfloat c = exp2(-8.698840e+03 * rcp_T + 0.000000e+00 * log_T + 2.864581e+01);
    cR[13] = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);}
    {const dfloat c = exp2(1.182859e+03 * rcp_T + 0.000000e+00 * log_T + 1.698815e+01);
    cR[14] = c * (concentrations[6]*concentrations[6] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[6]*exp_Gibbs0_RT[6]) * concentrations[1]*concentrations[7]);}
    {const dfloat Pr = exp2(-3.303265e+04 * rcp_T + 0.000000e+00 * log_T + 3.680665e+01) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]);
    const dfloat logFcent = log2(5.000000e-01 * exp2(-1.442695e+30*T) + 5.000000e-01 * exp2(-1.442695e-30*T) + exp2(0.000000e+00*rcp_T));
    const dfloat logPr_c = log2(Pr) - 0.67*logFcent - 1.328771e+00;
    const dfloat f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.491446e+00);
    const dfloat c = Pr / (exp2(3.515981e+04 * rcp_T - 0.000000e+00 * log_T - 4.806820e+01) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    cR[15] = c * (concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[7])* rcp_P0_RT * concentrations[5]*concentrations[5]);}
    {const dfloat c = exp2(-2.882190e+03 * rcp_T + 0.000000e+00 * log_T + 2.452253e+01);
    cR[16] = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[5]);}
    {const dfloat c = exp2(-5.771639e+03 * rcp_T + 0.000000e+00 * log_T + 2.552253e+01);
    cR[17] = c * (concentrations[3]*concentrations[7] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]) * concentrations[0]*concentrations[6]);}
    {const dfloat c = exp2(-2.882190e+03 * rcp_T + 2.000000e+00 * log_T + 3.255501e+00);
    cR[18] = c * (concentrations[4]*concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[7]) * concentrations[5]*concentrations[6]);}
    {const dfloat c = exp2(-0.000000e+00 * rcp_T + 0.000000e+00 * log_T + 1.993157e+01);
    cR[19] = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);}
    {const dfloat c = exp2(-6.938309e+03 * rcp_T + 0.000000e+00 * log_T + 2.911148e+01);
    cR[20] = c * (concentrations[5]*concentrations[7] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]) * concentrations[2]*concentrations[6]);}
    molar_rates[0] = -cR[1]-cR[2]-cR[4]+cR[9]+cR[17];
    molar_rates[1] = -cR[0]+cR[5]-cR[8]+cR[9]+cR[11]+cR[12]+cR[13]+cR[14];
    molar_rates[2] = cR[2]-cR[3]+cR[7]+cR[12]+cR[16]+cR[19]+cR[20];
    molar_rates[3] = -cR[0]+cR[1]+cR[2]+2*cR[4]-cR[6]-cR[7]-cR[8]-cR[9]-cR[10]-cR[16]-cR[17];
    molar_rates[4] = cR[0]-cR[1]-cR[3]-2*cR[5]-cR[6]-cR[11]-cR[18];
    molar_rates[5] = cR[0]+cR[1]-cR[2]+2*cR[3]+cR[6]-cR[7]+2*cR[10]+cR[11]-cR[12]+2*cR[15]+cR[16]+cR[18]-cR[19]-cR[20];
    molar_rates[6] = cR[8]-cR[9]-cR[10]-cR[11]-cR[12]-2*cR[13]-2*cR[14]+cR[17]+cR[18]+cR[19]+cR[20];
    molar_rates[7] = cR[13]+cR[14]-cR[15]-cR[16]-cR[17]-cR[18]-cR[19]-cR[20];
}

