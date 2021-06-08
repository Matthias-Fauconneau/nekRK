#define n_species 9
const dfloat fg_molar_mass[9] = {2.01594, 31.9988, 18.01534, 1.00797, 15.9994, 17.00737, 33.00677, 34.01474, 28.0134};
const dfloat fg_rcp_molar_mass[9] = {0.496046509321, 0.0312511719189, 0.0555082501912, 0.992093018641, 0.0625023438379, 0.0587980387326, 0.0302968148656, 0.0293990193663, 0.0356972020533};
void fg_molar_heat_capacity_at_constant_pressure_R(const dfloat ln_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
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

void fg_enthalpy_RT(const dfloat ln_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
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

void fg_exp_Gibbs_RT(const dfloat ln_T, const dfloat T, const dfloat T_2, const dfloat T_3, const dfloat T_4, const dfloat rcp_T, dfloat* species) {
    if (T < 1000.0) {
        species[0] = exp2(-1.01252087e+03 * rcp_T +6.59221840e+00 -3.29812431e+00 * ln_T +4.12472087e-04 * T +1.35716922e-07 * T_2 +7.89619527e-12 * T_3 -2.06743612e-14 * T_4);
        species[1] = exp2(-1.00524902e+03 * rcp_T -2.82180119e+00 -3.21293640e+00 * ln_T +5.63743175e-04 * T +9.59358412e-08 * T_2 -1.09489769e-10 * T_3 +4.38427696e-14 * T_4);
        species[2] = exp2(-3.02081133e+04 * rcp_T +7.96609640e-01 -3.38684249e+00 * ln_T +1.73749123e-03 * T +1.05911606e-06 * T_2 -5.80715106e-10 * T_3 +1.25329424e-13 * T_4);
        species[3] = exp2(+2.54716270e+04 * rcp_T +2.96011761e+00 -2.50000000e+00 * ln_T +0.00000000e+00 * T -0.00000000e+00 * T_2 -0.00000000e+00 * T_3 -0.00000000e+00 * T_4);
        species[4] = exp2(+2.91476445e+04 * rcp_T -1.75662000e-02 -2.94642878e+00 * ln_T -8.19083245e-04 * T -4.03505283e-07 * T_2 +1.33570266e-10 * T_3 -1.94534818e-14 * T_4);
        species[5] = exp2(+3.34630913e+03 * rcp_T +4.81573857e+00 -4.12530561e+00 * ln_T -1.61272470e-03 * T -1.08794115e-06 * T_2 +4.83211369e-10 * T_3 -1.03118689e-13 * T_4);
        species[6] = exp2(+2.94808040e+02 * rcp_T +5.85135560e-01 -4.30179801e+00 * ln_T -2.37456025e-03 * T -3.52638152e-06 * T_2 +2.02303245e-09 * T_3 -4.64612562e-13 * T_4);
        species[7] = exp2(-1.76631465e+04 * rcp_T -3.39660955e+00 -3.38875365e+00 * ln_T +3.28461290e-03 * T +2.47502097e-08 * T_2 +3.85483793e-10 * T_3 -1.23575738e-13 * T_4);
        species[8] = exp2(-1.02090000e+03 * rcp_T -6.51695000e-01 -3.29867700e+00 * ln_T +7.04120000e-04 * T +6.60537000e-07 * T_2 -4.70126250e-10 * T_3 +1.22242750e-13 * T_4);
    } else {
        species[0] = exp2(-8.35033997e+02 * rcp_T +4.34653354e+00 -2.99142337e+00 * ln_T +3.50032206e-04 * T +9.38971448e-09 * T_2 +7.69298182e-13 * T_3 -7.91375895e-17 * T_4);
        species[1] = exp2(-1.23393018e+03 * rcp_T +5.08412600e-01 -3.69757819e+00 * ln_T +3.06759845e-04 * T +2.09806998e-08 * T_2 -1.47940123e-12 * T_3 +5.68217655e-17 * T_4);
        species[2] = exp2(-2.98992090e+04 * rcp_T -4.19067120e+00 -2.67214561e+00 * ln_T +1.52814644e-03 * T +1.45504335e-07 * T_2 -1.00083033e-11 * T_3 +3.19580894e-16 * T_4);
        species[3] = exp2(+2.54716270e+04 * rcp_T +2.96011764e+00 -2.50000000e+00 * ln_T +0.00000000e+00 * T -0.00000000e+00 * T_2 -0.00000000e+00 * T_3 -0.00000000e+00 * T_4);
        species[4] = exp2(+2.92308027e+04 * rcp_T -2.37824845e+00 -2.54205966e+00 * ln_T -1.37753096e-05 * T +5.17133892e-10 * T_2 -3.79255618e-13 * T_3 +2.18402575e-17 * T_4);
        species[5] = exp2(+3.68362875e+03 * rcp_T -2.83691187e+00 -2.86472886e+00 * ln_T +5.28252240e-04 * T +4.31804597e-08 * T_2 -2.54348895e-12 * T_3 +6.65979380e-17 * T_4);
        species[6] = exp2(+1.11856713e+02 * rcp_T +2.32108750e-01 -4.01721090e+00 * ln_T +1.11991006e-03 * T +1.05609692e-07 * T_2 -9.52053083e-12 * T_3 +5.39542675e-16 * T_4);
        species[7] = exp2(-1.80069609e+04 * rcp_T +4.07202989e+00 -4.57316685e+00 * ln_T +2.16806820e-03 * T +2.45781470e-07 * T_2 -1.95741964e-11 * T_3 +7.15826780e-16 * T_4);
        species[8] = exp2(-9.22797700e+02 * rcp_T -3.05388800e+00 -2.92664000e+00 * ln_T +7.43988500e-04 * T +9.47460167e-08 * T_2 -8.41420000e-12 * T_3 +3.37667550e-16 * T_4);
    }
}

void fg_rates(const dfloat log_T, const dfloat T, const dfloat T2, const dfloat T4, const dfloat rcp_T, const dfloat rcp_T2, const dfloat P0_RT, const dfloat rcp_P0_RT, const dfloat exp_Gibbs0_RT[], const dfloat concentrations[], dfloat* molar_rates) {
    dfloat cR[21];
    cR[0] = 3.547000e+15*fg_exp2(-5.857342e-01*log_T-8.352941e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[1]*concentrations[3] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[4]*exp_Gibbs0_RT[5]) * concentrations[4]*concentrations[5];
    cR[1] = 5.080000e+04*fg_exp2(3.851996e+00*log_T-3.165251e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[0]*concentrations[4] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[4]/(exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5]) * concentrations[3]*concentrations[5];
    cR[2] = 2.160000e+08*fg_exp2(2.178470e+00*log_T-1.726043e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[0]*concentrations[5] - exp_Gibbs0_RT[0]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[3]) * concentrations[2]*concentrations[3];
    cR[3] = 2.970000e+06*fg_exp2(2.914244e+00*log_T-6.743142e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[2]*concentrations[4] - exp_Gibbs0_RT[2]*exp_Gibbs0_RT[4]/(exp_Gibbs0_RT[5])+ P0_RT * concentrations[5];
    cR[4] = 4.577000e+19*fg_exp2(-2.019773e+00*log_T-5.252605e+04*rcp_T) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[0] - exp_Gibbs0_RT[0]/(exp_Gibbs0_RT[3]) * concentrations[3];
    cR[5] = 6.165000e+15*fg_exp2(-7.213475e-01*log_T) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[4] - exp_Gibbs0_RT[4]/(exp_Gibbs0_RT[1]) * concentrations[1];
    cR[6] = 4.714000e+18*fg_exp2(-1.442695e+00*log_T) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[4] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[4]/(exp_Gibbs0_RT[5])+ P0_RT * concentrations[5];
    cR[7] = 3.800000e+22*fg_exp2(-2.885390e+00*log_T) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[5] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[5]/(exp_Gibbs0_RT[2])+ P0_RT * concentrations[2];
    cR[8] = 1.475000e+12*fg_exp2(8.656170e-01*log_T) * (2.000000e+00*concentrations[0] + 7.800000e-01*concentrations[1] + 1.100000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[1]*concentrations[3] - exp_Gibbs0_RT[1]*exp_Gibbs0_RT[3]/(exp_Gibbs0_RT[6])+ P0_RT * concentrations[6];
    cR[9] = 1.660000e+13*fg_exp2(-4.141497e+02*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[6] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[1]) * concentrations[0]*concentrations[1];
    cR[10] = 7.079000e+13*fg_exp2(-1.484498e+02*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[6] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[5])+ P0_RT * concentrations[5];
    cR[11] = 3.250000e+13 * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[4]*concentrations[6] - exp_Gibbs0_RT[4]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[5]) * concentrations[1]*concentrations[5];
    cR[12] = 2.890000e+13*fg_exp2(+2.501001e+02*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[5]*concentrations[6] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[2]) * concentrations[1]*concentrations[2];
    cR[13] = 4.200000e+14*fg_exp2(-6.029576e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[6] - exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7])+ rcp_P0_RT * concentrations[1]*concentrations[7];
    cR[14] = 1.300000e+11*fg_exp2(+8.198956e+02*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[6] - exp_Gibbs0_RT[6]/(exp_Gibbs0_RT[1]*exp_Gibbs0_RT[7])+ rcp_P0_RT * concentrations[1]*concentrations[7];
    cR[15] = 2.951000e+14*fg_exp2(-2.437092e+04*rcp_T) * (2.500000e+00*concentrations[0] + concentrations[1] + 1.200000e+01*concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[7] - exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[5]) * concentrations[5];
    cR[16] = 2.410000e+13*fg_exp2(-1.997782e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[7] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[5]) * concentrations[2]*concentrations[5];
    cR[17] = 4.820000e+13*fg_exp2(-4.000595e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[3]*concentrations[7] - exp_Gibbs0_RT[3]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[0]*exp_Gibbs0_RT[6]) * concentrations[0]*concentrations[6];
    cR[18] = 9.550000e+06*fg_exp2(2.885390e+00*log_T-1.997782e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[4]*concentrations[7] - exp_Gibbs0_RT[4]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[5]*exp_Gibbs0_RT[6]) * concentrations[5]*concentrations[6];
    cR[19] = 1.000000e+12 * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[5]*concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]) * concentrations[2]*concentrations[6];
    cR[20] = 5.800000e+14*fg_exp2(-4.809269e+03*rcp_T) * (concentrations[0] + concentrations[1] + concentrations[2] + concentrations[3] + concentrations[4] + concentrations[5] + concentrations[6] + concentrations[7] + concentrations[8]) * concentrations[5]*concentrations[7] - exp_Gibbs0_RT[5]*exp_Gibbs0_RT[7]/(exp_Gibbs0_RT[2]*exp_Gibbs0_RT[6]) * concentrations[2]*concentrations[6];
    molar_rates[0] = cR[1]+cR[2]+cR[4]-cR[9]-cR[17];
    molar_rates[1] = cR[0]-cR[5]+cR[8]-cR[9]-cR[11]-cR[12]-cR[13]-cR[14];
    molar_rates[2] = -cR[2]+cR[3]-cR[7]-cR[12]-cR[16]-cR[19]-cR[20];
    molar_rates[3] = cR[0]-cR[1]-cR[2]-cR[4]+cR[6]+cR[7]+cR[8]+cR[9]+cR[10]+cR[16]+cR[17];
    molar_rates[4] = -cR[0]+cR[1]+cR[3]+cR[5]+cR[6]+cR[11]+cR[18];
    molar_rates[5] = -cR[0]-cR[1]+cR[2]-cR[3]-cR[6]+cR[7]-cR[10]-cR[11]+cR[12]-cR[15]-cR[16]-cR[18]+cR[19]+cR[20];
    molar_rates[6] = -cR[8]+cR[9]+cR[10]+cR[11]+cR[12]+cR[13]+cR[14]-cR[17]-cR[18]-cR[19]-cR[20];
    molar_rates[7] = -cR[13]-cR[14]+cR[15]+cR[16]+cR[17]+cR[18]+cR[19]+cR[20];
}

