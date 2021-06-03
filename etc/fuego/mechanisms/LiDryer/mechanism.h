#define n_species 9
const int inert_specie= 8;
const dfloat fg_molar_mass[9] = {
    2.0159400000000001e+00 
    ,
    3.1998799999999999e+01 
    ,
    1.8015339999999998e+01 
    ,
    1.0079700000000000e+00 
    ,
    1.5999400000000000e+01 
    ,
    1.7007369999999998e+01 
    ,
    3.3006769999999996e+01 
    ,
    3.4014739999999996e+01 
    ,
    2.8013400000000001e+01 
};
const dfloat fg_rcp_molar_mass[9] = {
    4.9604650932071387e-01 
    ,
    3.1251171918946961e-02 
    ,
    5.5508250191225926e-02 
    ,
    9.9209301864142774e-01 
    ,
    6.2502343837893923e-02 
    ,
    5.8798038732620038e-02 
    ,
    3.0296814865556371e-02 
    ,
    2.9399019366310019e-02 
    ,
    3.5697202053303061e-02 
};
void gibbs_RT(dfloat * species, dfloat tc[])
{
    dfloat T = tc[1];
    /*species with midpoint at T=1.000000e+03 kelvin */
    if (T < 1.000000e+03) {
        /*species 0: H2 */
        species[0] =
            -1.01252087e+03 * tc[5]
            +6.59221840e+00
            -3.29812431e+00 * tc[0]
            -4.12472087e-04 * tc[1]
            +1.35716922e-07 * tc[2]
            +7.89619527e-12 * tc[3]
            -2.06743612e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524902e+03 * tc[5]
            -2.82180119e+00
            -3.21293640e+00 * tc[0]
            -5.63743175e-04 * tc[1]
            +9.59358412e-08 * tc[2]
            -1.09489769e-10 * tc[3]
            +4.38427696e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02081133e+04 * tc[5]
            +7.96609640e-01
            -3.38684249e+00 * tc[0]
            -1.73749123e-03 * tc[1]
            +1.05911606e-06 * tc[2]
            -5.80715106e-10 * tc[3]
            +1.25329424e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * tc[5]
            +2.96011761e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.91476445e+04 * tc[5]
            -1.75662000e-02
            -2.94642878e+00 * tc[0]
            +8.19083245e-04 * tc[1]
            -4.03505283e-07 * tc[2]
            +1.33570266e-10 * tc[3]
            -1.94534818e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * tc[5]
            +4.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * tc[5]
            +5.85135560e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631465e+04 * tc[5]
            -3.39660955e+00
            -3.38875365e+00 * tc[0]
            -3.28461290e-03 * tc[1]
            +2.47502097e-08 * tc[2]
            +3.85483793e-10 * tc[3]
            -1.23575738e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * tc[5]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35033997e+02 * tc[5]
            +4.34653354e+00
            -2.99142337e+00 * tc[0]
            -3.50032206e-04 * tc[1]
            +9.38971448e-09 * tc[2]
            +7.69298182e-13 * tc[3]
            -7.91375895e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393018e+03 * tc[5]
            +5.08412600e-01
            -3.69757819e+00 * tc[0]
            -3.06759845e-04 * tc[1]
            +2.09806998e-08 * tc[2]
            -1.47940123e-12 * tc[3]
            +5.68217655e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.98992090e+04 * tc[5]
            -4.19067120e+00
            -2.67214561e+00 * tc[0]
            -1.52814644e-03 * tc[1]
            +1.45504335e-07 * tc[2]
            -1.00083033e-11 * tc[3]
            +3.19580894e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * tc[5]
            +2.96011764e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.92308027e+04 * tc[5]
            -2.37824845e+00
            -2.54205966e+00 * tc[0]
            +1.37753096e-05 * tc[1]
            +5.17133892e-10 * tc[2]
            -3.79255618e-13 * tc[3]
            +2.18402575e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * tc[5]
            -2.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * tc[5]
            +2.32108750e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069609e+04 * tc[5]
            +4.07202989e+00
            -4.57316685e+00 * tc[0]
            -2.16806820e-03 * tc[1]
            +2.45781470e-07 * tc[2]
            -1.95741964e-11 * tc[3]
            +7.15826780e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * tc[5]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
}
/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(dfloat * species, dfloat tc[])
{
    dfloat T = tc[1];
    /*species with midpoint at T=1.000000e+03 kelvin */
    if (T < 1.000000e+03) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
}
void fg_speciesEnthalpy_RT(dfloat * species, dfloat tc[])
{
    dfloat T = tc[1];
    /*species with midpoint at T=1.000000e+03 kelvin */
    if (T < 1.000000e+03) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * tc[5];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * tc[5];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * tc[5];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * tc[5];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * tc[5];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * tc[5];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * tc[5];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * tc[5];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * tc[5];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * tc[5];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * tc[5];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * tc[5];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * tc[5];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * tc[5];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * tc[5];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * tc[5];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * tc[5];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * tc[5];
    }
}
void fg_rates(dfloat * sc, dfloat tc[], dfloat * wdot)
{
    dfloat qdot;
    int id; /*loop counter */
    dfloat mixture;/*mixture concentration */
    dfloat gibbs0_RT[9];
    dfloat rcp_Kc;                      /*equilibrium constant */
    dfloat k_f;                     /*forward reaction rate */
    dfloat k_r;                     /*reverse reaction rate */
    dfloat q_f;                     /*forward progress rate */
    dfloat q_r;                     /*reverse progress rate */
    dfloat phi_f;                   /*forward phase space factor */
    dfloat phi_r;                   /*reverse phase space factor */
    dfloat alpha;                   /*enhancement */
    dfloat redP;                    /*reduced pressure */
    dfloat logPred;                 /*log of above */
    dfloat F;                       /*fallof rate enhancement */
    dfloat F_troe;                  /*TROE intermediate */
    dfloat logFcent;                /*TROE intermediate */
    dfloat troe;                    /*TROE intermediate */
    dfloat troe_c;                  /*TROE intermediate */
    dfloat troe_n;                  /*TROE intermediate */
    const dfloat refC = (1.0132500000000000e+05 / 8.3145100000000003e+00) * tc[5];
    const dfloat rcp_refC = 1/refC;
    const dfloat T = tc[1];
    mixture = 0.0;
    for (id = 0; id < 9; ++id) {
        mixture += sc[id];
    }
    gibbs_RT(gibbs0_RT, tc);

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[3]*sc[1];
    k_f = 1.000000e-06 * 3.547000e+15*fgexp(-4.060000e-01*tc[0]-8.352893e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[3] + gibbs0_RT[1]) - (gibbs0_RT[4] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[4]*sc[0];
    k_f = 1.000000e-06 * 5.080000e+04*fgexp(2.670000e+00*tc[0]-3.165233e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[4] + gibbs0_RT[0]) - (gibbs0_RT[3] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[0]*sc[5];
    k_f = 1.000000e-06 * 2.160000e+08*fgexp(1.510000e+00*tc[0]-1.726033e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[3];
    rcp_Kc = fgexp(-((gibbs0_RT[0] + gibbs0_RT[5]) - (gibbs0_RT[2] + gibbs0_RT[3])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;

    /*reaction 4: O + H2O <=> OH + OH */
    phi_f = sc[4]*sc[2];
    k_f = 1.000000e-06 * 2.970000e+06*fgexp(2.020000e+00*tc[0]-6.743103e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[4] + gibbs0_RT[2]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    /*reaction 5: H2 + M <=> H + H + M */
    phi_f = sc[0];
    alpha = mixture + 1.500000e+00*sc[0] + 1.100000e+01*sc[2];
    k_f = 1.000000e-06 * alpha * 4.577000e+19*fgexp(-1.400000e+00*tc[0]-5.252576e+04*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    rcp_Kc = rcp_refC * fgexp(-((gibbs0_RT[0]) - (gibbs0_RT[3] + gibbs0_RT[3])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;

    /*reaction 6: O + O + M <=> O2 + M */
    phi_f = sc[4]*sc[4];
    alpha = mixture + 1.500000e+00*sc[0] + 1.100000e+01*sc[2];
    k_f = 1.000000e-12 * alpha * 6.165000e+15*fgexp(-5.000000e-01*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[1];
    rcp_Kc = refC * fgexp(-((gibbs0_RT[4] + gibbs0_RT[4]) - (gibbs0_RT[1])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[1] += qdot;

    /*reaction 7: O + H + M <=> OH + M */
    phi_f = sc[4]*sc[3];
    alpha = mixture + 1.500000e+00*sc[0] + 1.100000e+01*sc[2];
    k_f = 1.000000e-12 * alpha * 4.714000e+18*fgexp(-1.000000e+00*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    rcp_Kc = refC * fgexp(-((gibbs0_RT[4] + gibbs0_RT[3]) - (gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[3]*sc[5];
    alpha = mixture + 1.500000e+00*sc[0] + 1.100000e+01*sc[2];
    k_f = 1.000000e-12 * alpha * 3.800000e+22*fgexp(-2.000000e+00*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[2];
    rcp_Kc = refC * fgexp(-((gibbs0_RT[3] + gibbs0_RT[5]) - (gibbs0_RT[2])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + sc[0] + 1.000000e+01*sc[2] + -2.200000e-01*sc[1];
    k_f = 1.000000e-06 * 1.475000e+12*fgexp(6.000000e-01*tc[0]);
    redP = 1.0e-12 * alpha / k_f * 6.366000e+20*fgexp(-1.720000e+00*tc[0]-2.640881e+02*tc[5]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((2.000000e-01*fgexp(T*(-1.000000e+30)))+ (8.000000e-01*fgexp(T*(-1.000000e-30))));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    rcp_Kc = refC * fgexp(-((gibbs0_RT[3] + gibbs0_RT[1]) - (gibbs0_RT[6])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    phi_f = sc[6]*sc[3];
    k_f = 1.000000e-06 * 1.660000e+13*fgexp(-4.141473e+02*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[1];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[3]) - (gibbs0_RT[0] + gibbs0_RT[1])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[3] -= qdot;
    wdot[0] += qdot;
    wdot[1] += qdot;

    /*reaction 11: HO2 + H <=> OH + OH */
    phi_f = sc[6]*sc[3];
    k_f = 1.000000e-06 * 7.079000e+13*fgexp(-1.484489e+02*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[3]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    /*reaction 12: HO2 + O <=> O2 + OH */
    phi_f = sc[6]*sc[4];
    k_f = 1.000000e-06 * 3.250000e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[4]) - (gibbs0_RT[1] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[4] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    phi_f = sc[6]*sc[5];
    k_f = 1.000000e-06 * 2.890000e+13*fgexp(+2.500987e+02*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[1];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[5]) - (gibbs0_RT[2] + gibbs0_RT[1])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[5] -= qdot;
    wdot[2] += qdot;
    wdot[1] += qdot;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1.000000e-06 * 4.200000e+14*fgexp(-6.029542e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[6]) - (gibbs0_RT[7] + gibbs0_RT[1])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[1] += qdot;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1.000000e-06 * 1.300000e+11*fgexp(+8.198909e+02*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    rcp_Kc = fgexp(-((gibbs0_RT[6] + gibbs0_RT[6]) - (gibbs0_RT[7] + gibbs0_RT[1])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[1] += qdot;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 1.500000e+00*sc[0] + 1.100000e+01*sc[2];
    k_f = 1.000000e+00 * 2.951000e+14*fgexp(-2.437078e+04*tc[5]);
    redP = 1.0e-12 * alpha / k_f * 1.202000e+17*fgexp(-2.289636e+04*tc[5]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((5.000000e-01*fgexp(T*(-1.000000e+30)))+ (5.000000e-01*fgexp(T*(-1.000000e-30))));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    rcp_Kc = rcp_refC * fgexp(-((gibbs0_RT[7]) - (gibbs0_RT[5] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[3];
    k_f = 1.000000e-06 * 2.410000e+13*fgexp(-1.997770e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[3]) - (gibbs0_RT[2] + gibbs0_RT[5])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[3] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[3];
    k_f = 1.000000e-06 * 4.820000e+13*fgexp(-4.000572e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[3]) - (gibbs0_RT[6] + gibbs0_RT[0])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[0] += qdot;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    phi_f = sc[7]*sc[4];
    k_f = 1.000000e-06 * 9.550000e+06*fgexp(2.000000e+00*tc[0]-1.997770e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[6];
    rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[4]) - (gibbs0_RT[5] + gibbs0_RT[6])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = 1.000000e-06 * 1.000000e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[5]) - (gibbs0_RT[6] + gibbs0_RT[2])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[2] += qdot;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = 1.000000e-06 * 5.800000e+14*fgexp(-4.809242e+03*tc[5]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    rcp_Kc = fgexp(-((gibbs0_RT[7] + gibbs0_RT[5]) - (gibbs0_RT[6] + gibbs0_RT[2])));
    k_r = k_f * rcp_Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[2] += qdot;

    return;
}
dfloat fg_mean_specific_heat_at_CP_R(dfloat T, const dfloat* mole_fractions)
{
    dfloat tc[6]; 
    tc[0] = log(T); 
    tc[1] = T; 
    tc[2] = T*T; 
    tc[3] = T*T*T; 
    tc[4] = T*T*T*T; 
    tc[5] = 1./T; 
    dfloat result = 0; 
    dfloat cpor[9]; 
    cp_R(cpor, tc);
    result += cpor[0]*mole_fractions[0]; 
    result += cpor[1]*mole_fractions[1]; 
    result += cpor[2]*mole_fractions[2]; 
    result += cpor[3]*mole_fractions[3]; 
    result += cpor[4]*mole_fractions[4]; 
    result += cpor[5]*mole_fractions[5]; 
    result += cpor[6]*mole_fractions[6]; 
    result += cpor[7]*mole_fractions[7]; 
    result += cpor[8]*mole_fractions[8]; 
    return result;
}
void fg_viscosity_and_thermal_conductivity(dfloat T, const dfloat mole_fractions[], /*->*/ dfloat& viscosity, dfloat& thermal_conductivity) {
    dfloat ln_T = log(T);
    {
        dfloat sum = 0.;
        {
            dfloat y = 0.;
            y += 2.13652322E+00*pow(ln_T,0);
            y += 6.86724779E+00*pow(ln_T,1);
            y += -9.61029535E-01*pow(ln_T,2);
            y += 4.90901894E-02*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[0] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += 1.64738554E+01*pow(ln_T,0);
            y += -1.37574613E+00*pow(ln_T,1);
            y += 3.72164321E-01*pow(ln_T,2);
            y += -2.00906343E-02*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[1] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[2] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += 9.59325894E+00*pow(ln_T,0);
            y += 3.03342641E+00*pow(ln_T,1);
            y += -2.98187398E-01*pow(ln_T,2);
            y += 1.23483557E-02*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[3] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += 1.02760990E+01*pow(ln_T,0);
            y += 2.10197867E+00*pow(ln_T,1);
            y += -1.85769950E-01*pow(ln_T,2);
            y += 7.83598622E-03*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[4] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += 1.60039812E+01*pow(ln_T,0);
            y += -3.44252546E-01*pow(ln_T,1);
            y += 1.58988102E-01*pow(ln_T,2);
            y += -7.40753810E-03*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[5] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[6] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[7] * pow(viscosity_i, 6.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat viscosity_i = exp(y);
            sum += mole_fractions[8] * pow(viscosity_i, 6.);
        }
        viscosity = pow(sum, 1./6.);
    }
    {
        dfloat sum = 0.;
        {
            dfloat y = 0.;
            y += 2.13652322E+00*pow(ln_T,0);
            y += 6.86724779E+00*pow(ln_T,1);
            y += -9.61029535E-01*pow(ln_T,2);
            y += 4.90901894E-02*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[0] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += 1.64738554E+01*pow(ln_T,0);
            y += -1.37574613E+00*pow(ln_T,1);
            y += 3.72164321E-01*pow(ln_T,2);
            y += -2.00906343E-02*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[1] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[2] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += 9.59325894E+00*pow(ln_T,0);
            y += 3.03342641E+00*pow(ln_T,1);
            y += -2.98187398E-01*pow(ln_T,2);
            y += 1.23483557E-02*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[3] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += 1.02760990E+01*pow(ln_T,0);
            y += 2.10197867E+00*pow(ln_T,1);
            y += -1.85769950E-01*pow(ln_T,2);
            y += 7.83598622E-03*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[4] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += 1.60039812E+01*pow(ln_T,0);
            y += -3.44252546E-01*pow(ln_T,1);
            y += 1.58988102E-01*pow(ln_T,2);
            y += -7.40753810E-03*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[5] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[6] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[7] * pow(thermal_conductivity_i, 4.);
        }
        {
            dfloat y = 0.;
            y += NAN*pow(ln_T,0);
            y += NAN*pow(ln_T,1);
            y += NAN*pow(ln_T,2);
            y += NAN*pow(ln_T,3);
            dfloat thermal_conductivity_i = exp(y);
            sum += mole_fractions[8] * pow(thermal_conductivity_i, 4.);
        }
        thermal_conductivity = pow(sum, 1./4.);
    }
}
void fg_mixture_diffusion_coefficients(const dfloat mole_fractions[n_species], const dfloat mass_fractions[n_species], dfloat T, dfloat* Ddiag) {
    dfloat logT[3] = {log(T), pow(log(T), 2), pow(log(T), 3)};
    dfloat T32 = pow(T, 3./2.);
    Ddiag[0] = (1 - mass_fractions[0]) * mole_fractions[0] / ( 0.
    + mole_fractions[1] / (T32 * (7.26955388E-02 + -1.48471209E+00 * logT[0] + 1.32244920E+01 * logT[1] + -2.52001979E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (1.05975266E-01 + -2.42913315E+00 * logT[0] + 2.26396802E+01 * logT[1] + -5.38526479E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (2.23339098E-01 + -4.66731562E+00 * logT[0] + 4.10653807E+01 * logT[1] + -8.25393662E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (8.34607323E-02 + -1.65746900E+00 * logT[0] + 1.50598121E+01 * logT[1] + -2.67687992E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (8.31835154E-02 + -1.65196368E+00 * logT[0] + 1.50097906E+01 * logT[1] + -2.66798860E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (7.26297231E-02 + -1.48336789E+00 * logT[0] + 1.32125190E+01 * logT[1] + -2.51773826E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (7.25677536E-02 + -1.48210225E+00 * logT[0] + 1.32012458E+01 * logT[1] + -2.51559006E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (6.68180001E-02 + -1.35322764E+00 * logT[0] + 1.21183884E+01 * logT[1] + -2.26260874E+01 * logT[2]))
    );
    Ddiag[1] = (1 - mass_fractions[1]) * mole_fractions[1] / ( 0.
    + mole_fractions[0] / (T32 * (7.26955388E-02 + -1.48471209E+00 * logT[0] + 1.32244920E+01 * logT[1] + -2.52001979E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (9.56460784E-03 + -2.75765521E-01 * logT[0] + 3.97162421E+00 * logT[1] + -1.10256643E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (1.73076755E-01 + -3.83052027E+00 * logT[0] + 3.39364658E+01 * logT[1] + -7.65481205E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (3.99107144E-02 + -8.56998080E-01 * logT[0] + 7.49540468E+00 * logT[1] + -1.59690109E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (3.91143102E-02 + -8.39896986E-01 * logT[0] + 7.34583652E+00 * logT[1] + -1.56503549E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (2.69299922E-02 + -5.87212945E-01 * logT[0] + 5.15272716E+00 * logT[1] + -1.13141106E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (2.67328579E-02 + -5.82914398E-01 * logT[0] + 5.11500790E+00 * logT[1] + -1.12312885E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (2.66246920E-02 + -5.77715829E-01 * logT[0] + 5.06015099E+00 * logT[1] + -1.10072756E+01 * logT[2]))
    );
    Ddiag[2] = (1 - mass_fractions[2]) * mole_fractions[2] / ( 0.
    + mole_fractions[0] / (T32 * (1.05975266E-01 + -2.42913315E+00 * logT[0] + 2.26396802E+01 * logT[1] + -5.38526479E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (9.56460784E-03 + -2.75765521E-01 * logT[0] + 3.97162421E+00 * logT[1] + -1.10256643E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (3.75733740E-02 + -1.23138619E+00 * logT[0] + 2.09037094E+01 * logT[1] + -6.00403283E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (3.88179487E-02 + -9.29942193E-01 * logT[0] + 9.51055914E+00 * logT[1] + -2.37642829E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (3.82038513E-02 + -9.15230569E-01 * logT[0] + 9.36010272E+00 * logT[1] + -2.33883335E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (1.70980389E-02 + -4.36422285E-01 * logT[0] + 5.09891767E+00 * logT[1] + -1.33793276E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (1.70083536E-02 + -4.34133094E-01 * logT[0] + 5.07217203E+00 * logT[1] + -1.33091482E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (1.26119610E-02 + -3.38573587E-01 * logT[0] + 4.34624453E+00 * logT[1] + -1.17236455E+01 * logT[2]))
    );
    Ddiag[3] = (1 - mass_fractions[3]) * mole_fractions[3] / ( 0.
    + mole_fractions[0] / (T32 * (2.23339098E-01 + -4.66731562E+00 * logT[0] + 4.10653807E+01 * logT[1] + -8.25393662E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (1.73076755E-01 + -3.83052027E+00 * logT[0] + 3.39364658E+01 * logT[1] + -7.65481205E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (3.75733740E-02 + -1.23138619E+00 * logT[0] + 2.09037094E+01 * logT[1] + -6.00403283E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (2.31489002E-01 + -5.04908585E+00 * logT[0] + 4.43106693E+01 * logT[1] + -9.73472147E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (2.31082087E-01 + -5.04021050E+00 * logT[0] + 4.42327794E+01 * logT[1] + -9.71760965E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (1.72996032E-01 + -3.82873371E+00 * logT[0] + 3.39206377E+01 * logT[1] + -7.65124183E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (1.72920058E-01 + -3.82705227E+00 * logT[0] + 3.39057410E+01 * logT[1] + -7.64788168E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (1.64303441E-01 + -3.61903454E+00 * logT[0] + 3.19399232E+01 * logT[1] + -7.14445930E+01 * logT[2]))
    );
    Ddiag[4] = (1 - mass_fractions[4]) * mole_fractions[4] / ( 0.
    + mole_fractions[0] / (T32 * (8.34607323E-02 + -1.65746900E+00 * logT[0] + 1.50598121E+01 * logT[1] + -2.67687992E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (3.99107144E-02 + -8.56998080E-01 * logT[0] + 7.49540468E+00 * logT[1] + -1.59690109E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (3.88179487E-02 + -9.29942193E-01 * logT[0] + 9.51055914E+00 * logT[1] + -2.37642829E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (2.31489002E-01 + -5.04908585E+00 * logT[0] + 4.43106693E+01 * logT[1] + -9.73472147E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (5.43146795E-02 + -1.14632422E+00 * logT[0] + 1.00506473E+01 * logT[1] + -2.06474117E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (3.97070607E-02 + -8.52625048E-01 * logT[0] + 7.45715763E+00 * logT[1] + -1.58875253E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (3.95145188E-02 + -8.48490618E-01 * logT[0] + 7.42099742E+00 * logT[1] + -1.58104857E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (3.81125608E-02 + -8.13979512E-01 * logT[0] + 7.12030284E+00 * logT[1] + -1.50021508E+01 * logT[2]))
    );
    Ddiag[5] = (1 - mass_fractions[5]) * mole_fractions[5] / ( 0.
    + mole_fractions[0] / (T32 * (8.31835154E-02 + -1.65196368E+00 * logT[0] + 1.50097906E+01 * logT[1] + -2.66798860E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (3.91143102E-02 + -8.39896986E-01 * logT[0] + 7.34583652E+00 * logT[1] + -1.56503549E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (3.82038513E-02 + -9.15230569E-01 * logT[0] + 9.36010272E+00 * logT[1] + -2.33883335E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (2.31082087E-01 + -5.04021050E+00 * logT[0] + 4.42327794E+01 * logT[1] + -9.71760965E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (5.43146795E-02 + -1.14632422E+00 * logT[0] + 1.00506473E+01 * logT[1] + -2.06474117E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (3.89064880E-02 + -8.35434444E-01 * logT[0] + 7.30680661E+00 * logT[1] + -1.55672014E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (3.87099642E-02 + -8.31214512E-01 * logT[0] + 7.26989860E+00 * logT[1] + -1.54885686E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (3.73868055E-02 + -7.98479374E-01 * logT[0] + 6.98471505E+00 * logT[1] + -1.47164736E+01 * logT[2]))
    );
    Ddiag[6] = (1 - mass_fractions[6]) * mole_fractions[6] / ( 0.
    + mole_fractions[0] / (T32 * (7.26297231E-02 + -1.48336789E+00 * logT[0] + 1.32125190E+01 * logT[1] + -2.51773826E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (2.69299922E-02 + -5.87212945E-01 * logT[0] + 5.15272716E+00 * logT[1] + -1.13141106E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (1.70980389E-02 + -4.36422285E-01 * logT[0] + 5.09891767E+00 * logT[1] + -1.33793276E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (1.72996032E-01 + -3.82873371E+00 * logT[0] + 3.39206377E+01 * logT[1] + -7.65124183E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (3.97070607E-02 + -8.52625048E-01 * logT[0] + 7.45715763E+00 * logT[1] + -1.58875253E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (3.89064880E-02 + -8.35434444E-01 * logT[0] + 7.30680661E+00 * logT[1] + -1.55672014E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (2.65216974E-02 + -5.78310008E-01 * logT[0] + 5.07460489E+00 * logT[1] + -1.11425735E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (2.64342417E-02 + -5.73583345E-01 * logT[0] + 5.02395501E+00 * logT[1] + -1.09285389E+01 * logT[2]))
    );
    Ddiag[7] = (1 - mass_fractions[7]) * mole_fractions[7] / ( 0.
    + mole_fractions[0] / (T32 * (7.25677536E-02 + -1.48210225E+00 * logT[0] + 1.32012458E+01 * logT[1] + -2.51559006E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (2.67328579E-02 + -5.82914398E-01 * logT[0] + 5.11500790E+00 * logT[1] + -1.12312885E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (1.70083536E-02 + -4.34133094E-01 * logT[0] + 5.07217203E+00 * logT[1] + -1.33091482E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (1.72920058E-01 + -3.82705227E+00 * logT[0] + 3.39057410E+01 * logT[1] + -7.64788168E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (3.95145188E-02 + -8.48490618E-01 * logT[0] + 7.42099742E+00 * logT[1] + -1.58104857E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (3.87099642E-02 + -8.31214512E-01 * logT[0] + 7.26989860E+00 * logT[1] + -1.54885686E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (2.65216974E-02 + -5.78310008E-01 * logT[0] + 5.07460489E+00 * logT[1] + -1.11425735E+01 * logT[2]))
    + mole_fractions[8] / (T32 * (2.62538177E-02 + -5.69668414E-01 * logT[0] + 4.98966455E+00 * logT[1] + -1.08539474E+01 * logT[2]))
    );
    Ddiag[8] = (1 - mass_fractions[8]) * mole_fractions[8] / ( 0.
    + mole_fractions[0] / (T32 * (6.68180001E-02 + -1.35322764E+00 * logT[0] + 1.21183884E+01 * logT[1] + -2.26260874E+01 * logT[2]))
    + mole_fractions[1] / (T32 * (2.66246920E-02 + -5.77715829E-01 * logT[0] + 5.06015099E+00 * logT[1] + -1.10072756E+01 * logT[2]))
    + mole_fractions[2] / (T32 * (1.26119610E-02 + -3.38573587E-01 * logT[0] + 4.34624453E+00 * logT[1] + -1.17236455E+01 * logT[2]))
    + mole_fractions[3] / (T32 * (1.64303441E-01 + -3.61903454E+00 * logT[0] + 3.19399232E+01 * logT[1] + -7.14445930E+01 * logT[2]))
    + mole_fractions[4] / (T32 * (3.81125608E-02 + -8.13979512E-01 * logT[0] + 7.12030284E+00 * logT[1] + -1.50021508E+01 * logT[2]))
    + mole_fractions[5] / (T32 * (3.73868055E-02 + -7.98479374E-01 * logT[0] + 6.98471505E+00 * logT[1] + -1.47164736E+01 * logT[2]))
    + mole_fractions[6] / (T32 * (2.64342417E-02 + -5.73583345E-01 * logT[0] + 5.02395501E+00 * logT[1] + -1.09285389E+01 * logT[2]))
    + mole_fractions[7] / (T32 * (2.62538177E-02 + -5.69668414E-01 * logT[0] + 4.98966455E+00 * logT[1] + -1.08539474E+01 * logT[2]))
    );
}

/* version */
/* $Id$ */

/* Generated automatically by CPickler on Thu Jun  3 09:11:33 2021 */

/* End of file  */
