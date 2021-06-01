void fg_rho_Di(dfloat pressure, dfloat T, const dfloat mole_fractions[], const dfloat mass_fractions[], const dfloat mean_rcp_molar_mass, const dfloat mean_molar_mass, /*->*/
                                                dfloat* rho_Di) {
    const dfloat K = 1.380649e-23; // J / K
    const dfloat NA = 6.02214076e23;
    const dfloat R = K*NA;
    dfloat concentration = pressure/R/T;
    dfloat density = concentration * mean_molar_mass;
    dfloat logT[3] = {log(T), pow(log(T), 2), pow(log(T), 3)};
    fg_mixture_diffusion_coefficients(
        mean_rcp_molar_mass,
        /*Xloc:*/mole_fractions,
        /*Yloc:*/mass_fractions,
        logT,
        rholoc,
        Tloc,
        /*Ddiag:*/ rho_Di // Is "Ddiag" <=> "rho_Di" ?
    );
    for (int i = 0; i < NUM_SPECIES; ++i) {
        Ddiag[i] *= density;
    }
}
