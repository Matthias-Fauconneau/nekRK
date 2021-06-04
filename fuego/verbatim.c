dfloat fg_mean_specific_heat_at_CP_R(dfloat T, const dfloat* mole_fractions) {
    dfloat Ts[] = { log(T), T, T*T, T*T*T, T*T*T*T, 1./T };
    dfloat Cp_R[n_species];
    molar_heat_capacity_at_constant_pressure_R(Ts, Cp_R);
    dfloat sum = 0;
    for (int i = 0; i < n_species; i++) {
        sum += Cp_R[i]*mole_fractions[i]
    }
}

void fg_rho_Di(dfloat pressure, dfloat T, const dfloat mole_fractions[], const dfloat mass_fractions[], const dfloat mean_molar_mass, dfloat* rho_Di) {
    fg_mixture_diffusion_coefficients(
        mole_fractions,
        mass_fractions,
        T,
        rho_Di
    );
    const dfloat K = 1.380649e-23; // J / K
    const dfloat NA = 6.02214076e23;
    const dfloat R = K*NA;
    dfloat concentration = pressure/R/T;
    dfloat density = concentration * mean_molar_mass;
    for (int i = 0; i < n_species; i++) {
        rho_Di[i] *= density;
    }
}
