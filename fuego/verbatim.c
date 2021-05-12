#define NUM_SPECIES n_species
// No idea what this code does
void Pele_binary_diffusion_coefficients_Ddiag(
    const dfloat wbar,
    const dfloat Xloc[NUM_SPECIES],
    const dfloat Yloc[NUM_SPECIES],
    dfloat logT[3],
    dfloat rholoc,
    dfloat Tloc,
    dfloat* Ddiag,
    dfloat trans_fitdbin[NUM_SPECIES*NUM_SPECIES*4],
    const dfloat trans_wt[NUM_SPECIES] // What is this ?
) {
    for (int i = 0; i < NUM_SPECIES; ++i) {
        dfloat term1 = 0.0;
        dfloat term2 = 0.0;
        for (int j = 0; j < NUM_SPECIES; ++j) {
            if (i != j) {
                dfloat dbintemp =
                    trans_fitdbin[4 * (i + NUM_SPECIES * j)] +
                    trans_fitdbin[1 + 4 * (i + NUM_SPECIES * j)] * logT[0] +
                    trans_fitdbin[2 + 4 * (i + NUM_SPECIES * j)] * logT[1] +
                    trans_fitdbin[3 + 4 * (i + NUM_SPECIES * j)] * logT[2];
                term1 = term1 + Yloc[j];
                term2 = term2 + Xloc[j] / exp(dbintemp);
            }
        }
        Ddiag[i] = trans_wt[i] * term1 / term2 / wbar;
    }
    // Call CKRP ?
    const dfloat atmospheric_pressure = 101325.;
    const dfloat RU = 8.31447e7; // R . e7 ?
    const dfloat pscale = atmospheric_pressure * wbar / (RU * Tloc * rholoc);
    for (int i = 0; i < NUM_SPECIES; ++i) {
        Ddiag[i] = rholoc * pscale * Ddiag[i];
    }
}

void fg_rho_Di(dfloat pressure, dfloat T, const dfloat mass_fractions[], /*->*/ dfloat* rho_Di) {
    dfloat mean_rcp_molar_mass = 0.;
    for (int i = 0; i < NUM_SPECIES; ++i) {
        mean_rcp_molar_mass += mass_fractions[i]*fg_rcp_molar_mass[i];
    }
    dfloat mean_molar_mass = 1./mean_rcp_molar_mass;
    dfloat mole_fractions[n_species];
    for (int i = 0; i < NUM_SPECIES; ++i) {
        mole_fractions[i] = mass_fractions[i]*fg_rcp_molar_mass[i]*mean_molar_mass;
    }
    const dfloat K = 1.380649e-23; // J / K
    const dfloat NA = 6.02214076e23;
    const dfloat R = K*NA;
    dfloat concentration = pressure/R/T;
    dfloat density = concentration * mean_molar_mass;
    dfloat wbar = mean_molar_mass;
    dfloat logT[3] = {log(T), pow(log(T), 2), pow(log(T), 3)};
    dfloat rholoc = density;
    dfloat Tloc = T;
    Pele_binary_diffusion_coefficients_Ddiag(
        wbar,
        /*Xloc:*/mole_fractions,
        /*Yloc:*/mass_fractions,
        logT,
        rholoc,
        Tloc,
        /*Ddiag:*/ rho_Di, // Is "Ddiag" <=> "rho_Di" ?
        /*trans_fitdbin:*/fg_binary_diffusion_coefficients,
        /*trans_wt:*/fg_molar_mass
    );
}
