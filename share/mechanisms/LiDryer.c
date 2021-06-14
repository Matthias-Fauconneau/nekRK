#define n_species 9
const float fg_molar_mass[9] = {0.002016, 0.031998, 0.018015, 0.001008, 0.015999, 0.017006999999999998, 0.033006, 0.034013999999999996, 0.0280134};
const float fg_rcp_molar_mass[9] = {496.031746031746, 31.251953247077942, 55.50929780738274, 992.063492063492, 62.503906494155885, 58.79931792791204, 30.297521662727988, 29.39965896395602, 35.69720205330306};
void fg_molar_heat_capacity_at_constant_pressure_R(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {
 if (T < 1000.0) {
    _[0] = 3.29812431 + 0.000824944174 * T -8.14301529e-07 * T_2 -9.47543433e-11 * T_3 + 4.13487224e-13 * T_4;
    _[1] = 3.2129364 + 0.00112748635 * T -5.75615047e-07 * T_2 + 1.31387723e-09 * T_3 -8.76855392e-13 * T_4;
    _[2] = 3.38684249 + 0.00347498246 * T -6.35469633e-06 * T_2 + 6.96858127e-09 * T_3 -2.50658847e-12 * T_4;
    _[3] = 2.5 + 0.0 * T + 0.0 * T_2 + 0.0 * T_3 + 0.0 * T_4;
    _[4] = 2.94642878 -0.00163816649 * T + 2.4210317e-06 * T_2 -1.60284319e-09 * T_3 + 3.89069636e-13 * T_4;
    _[5] = 4.12530561 -0.00322544939 * T + 6.52764691e-06 * T_2 -5.79853643e-09 * T_3 + 2.06237379e-12 * T_4;
    _[6] = 4.30179801 -0.00474912051 * T + 2.11582891e-05 * T_2 -2.42763894e-08 * T_3 + 9.29225124e-12 * T_4;
    _[7] = 3.38875365 + 0.00656922581 * T -1.48501258e-07 * T_2 -4.62580552e-09 * T_3 + 2.47151475e-12 * T_4;
    _[8] = 3.298677 + 0.00140824 * T -3.963222e-06 * T_2 + 5.641515e-09 * T_3 -2.444855e-12 * T_4;
 } else {
    _[0] = 2.99142337 + 0.000700064411 * T -5.63382869e-08 * T_2 -9.23157818e-12 * T_3 + 1.58275179e-15 * T_4;
    _[1] = 3.69757819 + 0.000613519689 * T -1.25884199e-07 * T_2 + 1.77528148e-11 * T_3 -1.13643531e-15 * T_4;
    _[2] = 2.67214561 + 0.00305629289 * T -8.73026011e-07 * T_2 + 1.20099639e-10 * T_3 -6.39161787e-15 * T_4;
    _[3] = 2.5 + 0.0 * T + 0.0 * T_2 + 0.0 * T_3 + 0.0 * T_4;
    _[4] = 2.54205966 -2.75506191e-05 * T -3.10280335e-09 * T_2 + 4.55106742e-12 * T_3 -4.3680515e-16 * T_4;
    _[5] = 2.86472886 + 0.00105650448 * T -2.59082758e-07 * T_2 + 3.05218674e-11 * T_3 -1.33195876e-15 * T_4;
    _[6] = 4.0172109 + 0.00223982013 * T -6.3365815e-07 * T_2 + 1.1424637e-10 * T_3 -1.07908535e-14 * T_4;
    _[7] = 4.57316685 + 0.00433613639 * T -1.47468882e-06 * T_2 + 2.34890357e-10 * T_3 -1.43165356e-14 * T_4;
    _[8] = 2.92664 + 0.001487977 * T -5.684761e-07 * T_2 + 1.009704e-10 * T_3 -6.753351e-15 * T_4;
 }
}
void fg_enthalpy_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {
 if (T < 1000.0) {
    _[0] = 3.29812431 + 0.000412472087 * T -2.71433843e-07 * T_2 -2.3688585825e-11 * T_3 + 8.26974448e-14 * T_4 -1012.52087 * rcp_T;
    _[1] = 3.2129364 + 0.000563743175 * T -1.9187168233333332e-07 * T_2 + 3.284693075e-10 * T_3 -1.753710784e-13 * T_4 -1005.24902 * rcp_T;
    _[2] = 3.38684249 + 0.00173749123 * T -2.11823211e-06 * T_2 + 1.7421453175e-09 * T_3 -5.01317694e-13 * T_4 -30208.1133 * rcp_T;
    _[3] = 2.5 + 0.0 * T + 0.0 * T_2 + 0.0 * T_3 + 0.0 * T_4 + 25471.627 * rcp_T;
    _[4] = 2.94642878 -0.000819083245 * T + 8.070105666666667e-07 * T_2 -4.007107975e-10 * T_3 + 7.78139272e-14 * T_4 + 29147.6445 * rcp_T;
    _[5] = 4.12530561 -0.001612724695 * T + 2.1758823033333334e-06 * T_2 -1.4496341075e-09 * T_3 + 4.1247475799999997e-13 * T_4 + 3346.30913 * rcp_T;
    _[6] = 4.30179801 -0.002374560255 * T + 7.0527630333333326e-06 * T_2 -6.06909735e-09 * T_3 + 1.8584502480000002e-12 * T_4 + 294.80804 * rcp_T;
    _[7] = 3.38875365 + 0.003284612905 * T -4.950041933333333e-08 * T_2 -1.15645138e-09 * T_3 + 4.9430295e-13 * T_4 -17663.1465 * rcp_T;
    _[8] = 3.298677 + 0.00070412 * T -1.3210739999999999e-06 * T_2 + 1.41037875e-09 * T_3 -4.889710000000001e-13 * T_4 -1020.9 * rcp_T;
 } else {
    _[0] = 2.99142337 + 0.0003500322055 * T -1.8779428966666666e-08 * T_2 -2.307894545e-12 * T_3 + 3.16550358e-16 * T_4 -835.033997 * rcp_T;
    _[1] = 3.69757819 + 0.0003067598445 * T -4.1961399666666664e-08 * T_2 + 4.4382037e-12 * T_3 -2.27287062e-16 * T_4 -1233.93018 * rcp_T;
    _[2] = 2.67214561 + 0.001528146445 * T -2.910086703333333e-07 * T_2 + 3.002490975e-11 * T_3 -1.278323574e-15 * T_4 -29899.209 * rcp_T;
    _[3] = 2.5 + 0.0 * T + 0.0 * T_2 + 0.0 * T_3 + 0.0 * T_4 + 25471.627 * rcp_T;
    _[4] = 2.54205966 -1.377530955e-05 * T -1.0342677833333334e-09 * T_2 + 1.137766855e-12 * T_3 -8.736102999999999e-17 * T_4 + 29230.8027 * rcp_T;
    _[5] = 2.86472886 + 0.00052825224 * T -8.636091933333334e-08 * T_2 + 7.63046685e-12 * T_3 -2.66391752e-16 * T_4 + 3683.62875 * rcp_T;
    _[6] = 4.0172109 + 0.001119910065 * T -2.1121938333333332e-07 * T_2 + 2.85615925e-11 * T_3 -2.1581707e-15 * T_4 + 111.856713 * rcp_T;
    _[7] = 4.57316685 + 0.002168068195 * T -4.9156294e-07 * T_2 + 5.872258925e-11 * T_3 -2.86330712e-15 * T_4 -18006.9609 * rcp_T;
    _[8] = 2.92664 + 0.0007439885 * T -1.8949203333333333e-07 * T_2 + 2.52426e-11 * T_3 -1.3506701999999999e-15 * T_4 -922.7977 * rcp_T;
 }
}
void fg_exp_Gibbs_RT(const float log_T, const float T, const float T_2, const float T_3, const float T_4, const float rcp_T, float* _) {
 if (T < 1000.0) {
    _[0] = exp2(-1460.7588379455788 * rcp_T + 9.510560794136978 -3.29812431 * log_T -0.0005950714344200211 * T + 1.9579812961276678e-07 * T_2 + 1.1391801765133363e-11 * T_3 -2.9826798376787193e-14 * T_4);
    _[1] = exp2(-1450.2677760124905 * rcp_T -4.070998583187576 -3.2129364 * log_T -0.000813309482907499 * T + 1.3840616229466128e-07 * T_2 -1.579603470048273e-10 * T_3 + 6.325174628075739e-14 * T_4);
    _[2] = exp2(-43581.09525252194 * rcp_T + 1.149264777152342 -3.38684249 * log_T -0.0025066699811090655 * T + 1.5279814802743832e-06 * T_2 -8.377948033550593e-10 * T_3 + 1.8081213776092263e-13 * T_4);
    _[3] = exp2(36747.789956273424 * rcp_T + 4.270546993509701 -2.5 * log_T -0.0 * T -0.0 * T_2 -0.0 * T_3 -0.0 * T_4);
    _[4] = exp2(42051.162173744466 * rcp_T -0.025342669627263194 -2.94642878 * log_T + 0.0011816873356367398 * T -5.821350712374962e-07 * T_2 + 1.9270116012797052e-10 * T_3 -2.8065441720883702e-14 * T_4);
    _[5] = exp2(4827.7035871324615 * rcp_T + 6.947642153156708 -4.12530561 * log_T + 0.0023266699197956662 * T -1.569567304288528e-06 * T_2 + 6.97126645997916e-10 * T_3 -1.4876882196461877e-13 * T_4);
    _[6] = exp2(425.3180973221952 * rcp_T + 0.8441721706597863 -4.30179801 * log_T + 0.003425766304180532 * T -5.087493126377502e-06 * T_2 + 2.9186188831724494e-09 * T_3 -6.70294239132116e-13 * T_4);
    _[7] = exp2(-25482.53386204525 * rcp_T -4.900271753621094 -3.38875365 * log_T -0.004738694749283392 * T + 3.570700474706209e-08 * T_2 + 5.56135556985066e-10 * T_3 -1.7828210366544628e-13 * T_4);
    _[8] = exp2(-1472.8473672435427 * rcp_T -0.9401971446721332 -3.298677 * log_T -0.001015830432190737 * T + 9.529534542236733e-07 * T_2 -6.78248809466725e-10 * T_3 + 1.763590092096293e-13 * T_4);
 } else {
    _[0] = exp2(-1204.6994064455896 * rcp_T + 6.270722383215552 -2.99142337 * log_T -0.0005049897270262766 * T + 1.3546494520468278e-08 * T_2 + 1.1098626716553966e-12 * T_3 -1.1417140791955647e-16 * T_4);
    _[1] = exp2(-1780.184951489226 * rcp_T + 0.7334843367454645 -3.69757819 * log_T -0.0004425609064040196 * T + 3.026875160392991e-08 * T_2 -2.134324822815016e-12 * T_3 + 8.197647930140558e-17 * T_4);
    _[2] = exp2(-43135.44055080266 * rcp_T -6.045860558236202 -2.67214561 * log_T -0.0022046492979535992 * T + 2.0991838277279561e-07 * T_2 -1.4438929466487895e-11 * T_3 + 4.610577702153138e-16 * T_4);
    _[3] = exp2(36747.789956273424 * rcp_T + 4.270547036790552 -2.5 * log_T -0.0 * T -0.0 * T_2 -0.0 * T_3 -0.0 * T_4);
    _[4] = exp2(42171.134096493726 * rcp_T -3.4310872448168634 -2.54205966 * log_T + 1.987357077449538e-05 * T + 7.460665009831105e-10 * T_2 -5.471501997987774e-13 * T_3 + 3.1508831186987985e-17 * T_4);
    _[5] = exp2(5314.352930101011 * rcp_T -4.092798686288036 -2.86472886 * log_T -0.0007621068869864866 * T + 6.229623502440591e-08 * T_2 -3.669478894720876e-12 * T_3 + 9.608051488603063e-17 * T_4);
    _[6] = exp2(161.37512513524004 * rcp_T + 0.3348621425719363 -4.0172109 * log_T -0.0016156886970171366 * T + 1.523625784373125e-07 * T_2 -1.3735222619880468e-11 * T_3 + 7.783955415704656e-16 * T_4);
    _[7] = exp2(-25978.553191911462 * rcp_T + 5.874697330097326 -4.57316685 * log_T -0.0031278612332355862 * T + 3.5458770791139956e-07 * T_2 -2.8239596099711515e-11 * T_3 + 1.0327197456415149e-15 * T_4);
    _[8] = exp2(-1331.3156655337414 * rcp_T -4.405829073030314 -2.92664 * log_T -0.0010733485194284185 * T + 1.3668960838898307e-07 * T_2 -1.2139124613047913e-11 * T_3 + 4.87151299854126e-16 * T_4);
 }
}
float fg_viscosity(float T, const float mole_fractions[]) {
 float T_14 =	sqrt(sqrt(T));
 float ln_T = log(T);
 float ln_T_2 = ln_T*ln_T;
 float ln_T_3 = ln_T_2*ln_T;
 return pow(
    mole_fractions[0] * pow((3.967355530742668e-08 -5.492366609241516e-07 * ln_T + 3.290575418242205e-06 * ln_T_2 -6.101012720850198e-06 * ln_T_3)*T_14, 12)+
    mole_fractions[1] * pow((1.337797421803962e-07 -2.193806254685158e-06 * ln_T + 1.4616622602222706e-05 * ln_T_2 -3.18629994963365e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[2] * pow((-4.008877968531974e-09 + 1.4987755605312997e-07 * ln_T -1.2345787016867175e-06 * ln_T_2 + 3.0391139720306797e-06 * ln_T_3)*T_14, 12)+
    mole_fractions[3] * pow((6.747396472184357e-08 -1.1456762578245252e-06 * ln_T + 7.88497039779098e-06 * ln_T_2 -1.7941586103841824e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[4] * pow((1.4163828841054896e-07 -2.2111839073867454e-06 * ln_T + 1.4254310198737872e-05 * ln_T_2 -2.9677700314938626e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[5] * pow((1.460320247050356e-07 -2.2797766523054733e-06 * ln_T + 1.4696490634380591e-05 * ln_T_2 -3.059832700757898e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[6] * pow((1.358705662256428e-07 -2.228092932123417e-06 * ln_T + 1.4845063661376895e-05 * ln_T_2 -3.2360981660264593e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[7] * pow((1.3792969990215125e-07 -2.2618599305129267e-06 * ln_T + 1.507004226685546e-05 * ln_T_2 -3.2851415968389787e-05 * ln_T_3)*T_14, 12)+
    mole_fractions[8] * pow((1.125917155902469e-07 -1.8192640475814333e-06 * ln_T + 1.1992979824329308e-05 * ln_T_2 -2.5767884623198545e-05 * ln_T_3)*T_14, 12), 1./6.);
}

float fg_thermal_conductivity(float T, const float mole_fractions[]) {
 float T_12 =	sqrt(T);
 float ln_T = log(T);
 float ln_T_2 = ln_T*ln_T;
 float ln_T_3 = ln_T_2*ln_T;
 return pow(
    T_12 * pow((0.00017265705282316272 -0.0022445840583933192*ln_T + 0.009876262597026984*ln_T_2 -0.004398632788062826*ln_T_3), 4) * mole_fractions[0]+
    T_12 * pow((2.0626586379486266e-05 -0.00039256294559316345*ln_T + 0.0031513024774719425*ln_T_2 -0.007518156975105011*ln_T_3), 4) * mole_fractions[1]+
    T_12 * pow((-0.00013329872125751008 + 0.003479212920528736*ln_T -0.027081652911999904*ln_T_2 + 0.06637165379943971*ln_T_3), 4) * mole_fractions[2]+
    T_12 * pow((0.0002707289626737129 -0.005962840756360254*ln_T + 0.04743227428711307*ln_T_2 -0.10980473426337048*ln_T_3), 4) * mole_fractions[3]+
    T_12 * pow((2.167237619717602e-05 -0.0004035512424930518*ln_T + 0.002948558201776504*ln_T_2 -0.004784140309819591*ln_T_3), 4) * mole_fractions[4]+
    T_12 * pow((-5.9273488766219073e-05 + 0.0016079245276583707*ln_T -0.01242063119977092*ln_T_2 + 0.03301328279311626*ln_T_3), 4) * mole_fractions[5]+
    T_12 * pow((1.667734703320205e-06 + 9.944376454316164e-05*ln_T -0.0003659069293791121*ln_T_2 + 0.00023605613197563577*ln_T_3), 4) * mole_fractions[6]+
    T_12 * pow((-9.594988102273465e-05 + 0.0020521849555358027*ln_T -0.01285826945296755*ln_T_2 + 0.026539243043749917*ln_T_3), 4) * mole_fractions[7]+
    T_12 * pow((-7.133505112581695e-05 + 0.0015269432255689602*ln_T -0.010192619586375469*ln_T_2 + 0.023232150668911666*ln_T_3), 4) * mole_fractions[8], 1./4.);
}

void fg_mixture_diffusion_coefficients(float T, const float mole_fractions[], const float mass_fractions[], float* _) {
 float ln_T = log(T);
 float ln_T_2 = ln_T*ln_T;
 float ln_T_3 = ln_T_2*ln_T;
 float T_32 = T*sqrt(T);
 _[0] = (1. - mass_fractions[0]) * mole_fractions[0] / (
    mole_fractions[0] / ((7.877690777909358e-06 -0.00013031015412669209*ln_T + 0.001208169975705796*ln_T_2 -0.0012661700540825235*ln_T_3)*T_32)+
    mole_fractions[1] / ((8.624071879765021e-06 -0.00016891320367856734*ln_T + 0.0014082221783546169*ln_T_2 -0.002561292591037956*ln_T_3)*T_32)+
    mole_fractions[2] / ((2.1798065820502298e-05 -0.00048132682963867733*ln_T + 0.003954767938153096*ln_T_2 -0.009228619064248676*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.8052869518166994e-05 -0.0005659080737233624*ln_T + 0.004648814427229764*ln_T_2 -0.00909162904085566*ln_T_3)*T_32)+
    mole_fractions[4] / ((9.371675695604753e-06 -0.00017684142019394345*ln_T + 0.0015097675313646263*ln_T_2 -0.0024820837114006055*ln_T_3)*T_32)+
    mole_fractions[5] / ((9.340544322994412e-06 -0.0001762539781692879*ln_T + 0.001504752298539638*ln_T_2 -0.002473838582633779*ln_T_3)*T_32)+
    mole_fractions[6] / ((8.616263166230654e-06 -0.0001687602602849943*ln_T + 0.0014069470958023526*ln_T_2 -0.002558973454509375*ln_T_3)*T_32)+
    mole_fractions[7] / ((8.608910807271654e-06 -0.00016861625516494404*ln_T + 0.0014057465312553913*ln_T_2 -0.002556789852286297*ln_T_3)*T_32)+
    mole_fractions[8] / ((7.753386136651787e-06 -0.00015012001517394324*ln_T + 0.001260489007680266*ln_T_2 -0.0022244441293002544*ln_T_3)*T_32));
    _[1] = (1. - mass_fractions[1]) * mole_fractions[1] / (
    mole_fractions[0] / ((8.624071879764066e-06 -0.00016891320367854693*ln_T + 0.0014082221783544716*ln_T_2 -0.002561292591037609*ln_T_3)*T_32)+
    mole_fractions[1] / ((4.070043194962521e-06 -8.66300840660579e-05*ln_T + 0.0007032784791925587*ln_T_2 -0.0015415967056044923*ln_T_3)*T_32)+
    mole_fractions[2] / ((5.253952172982132e-06 -0.0001282368964392761*ln_T + 0.0011761372700424147*ln_T_2 -0.003053936869071254*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.819396161258144e-05 -0.0006108232939658377*ln_T + 0.004975808066492487*ln_T_2 -0.011259105125226525*ln_T_3)*T_32)+
    mole_fractions[4] / ((5.569865822310673e-06 -0.00011630889857421319*ln_T + 0.0009446644908807677*ln_T_2 -0.0019946621508287895*ln_T_3)*T_32)+
    mole_fractions[5] / ((5.45871533159476e-06 -0.00011398787477156308*ln_T + 0.0009258130633826859*ln_T_2 -0.0019548573002363357*ln_T_3)*T_32)+
    mole_fractions[6] / ((4.0388489733899366e-06 -8.596612107899355e-05*ln_T + 0.0006978883091977138*ln_T_2 -0.0015297813741923412*ln_T_3)*T_32)+
    mole_fractions[7] / ((4.009282075250514e-06 -8.533679535720242e-05*ln_T + 0.0006927793307024408*ln_T_2 -0.001518582418657164*ln_T_3)*T_32)+
    mole_fractions[8] / ((3.8904263458668995e-06 -8.230342706964434e-05*ln_T + 0.0006680176844741743*ln_T_2 -0.001447416314766726*ln_T_3)*T_32));
    _[2] = (1. - mass_fractions[2]) * mole_fractions[2] / (
    mole_fractions[0] / ((2.1798065820502298e-05 -0.00048132682963867733*ln_T + 0.003954767938153096*ln_T_2 -0.009228619064248676*ln_T_3)*T_32)+
    mole_fractions[1] / ((5.253952172982132e-06 -0.0001282368964392761*ln_T + 0.0011761372700424147*ln_T_2 -0.003053936869071254*ln_T_3)*T_32)+
    mole_fractions[2] / ((-2.4936442456881388e-05 + 0.0005416098291005012*ln_T -0.003572136447540915*ln_T_2 + 0.007706933189849318*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.0145742105206148e-05 -0.000541010493787195*ln_T + 0.005555473146463371*ln_T_2 -0.015315716453509364*ln_T_3)*T_32)+
    mole_fractions[4] / ((9.591588708458134e-06 -0.00022431069595172284*ln_T + 0.001951005759173615*ln_T_2 -0.004879738156005005*ln_T_3)*T_32)+
    mole_fractions[5] / ((9.43984208546376e-06 -0.00022076192091071393*ln_T + 0.0019201392839320014*ln_T_2 -0.004802536786265476*ln_T_3)*T_32)+
    mole_fractions[6] / ((5.224980566221356e-06 -0.00012752972921163517*ln_T + 0.0011696509910805045*ln_T_2 -0.0030370940335783104*ln_T_3)*T_32)+
    mole_fractions[7] / ((5.1975721852789544e-06 -0.00012686075382399255*ln_T + 0.0011635154199473926*ln_T_2 -0.003021162523561163*ln_T_3)*T_32)+
    mole_fractions[8] / ((5.5498678214477205e-06 -0.00013317697631526514*ln_T + 0.0011956571686973264*ln_T_2 -0.003062175134485725*ln_T_3)*T_32));
    _[3] = (1. - mass_fractions[3]) * mole_fractions[3] / (
    mole_fractions[0] / ((2.8052869518166994e-05 -0.0005659080737233624*ln_T + 0.004648814427229764*ln_T_2 -0.00909162904085566*ln_T_3)*T_32)+
    mole_fractions[1] / ((2.8193961612580416e-05 -0.0006108232939658158*ln_T + 0.004975808066492328*ln_T_2 -0.011259105125226142*ln_T_3)*T_32)+
    mole_fractions[2] / ((2.0145742105206148e-05 -0.000541010493787195*ln_T + 0.005555473146463371*ln_T_2 -0.015315716453509364*ln_T_3)*T_32)+
    mole_fractions[3] / ((7.490062729482605e-05 -0.001650528584029963*ln_T + 0.013542731391787518*ln_T_2 -0.031507680103450036*ln_T_3)*T_32)+
    mole_fractions[4] / ((3.4766888330436316e-05 -0.0007402656486890132*ln_T + 0.006009838481386653*ln_T_2 -0.013182221813085888*ln_T_3)*T_32)+
    mole_fractions[5] / ((3.4705768325682896e-05 -0.0007389642656162268*ln_T + 0.005999273217303768*ln_T_2 -0.013159047537257077*ln_T_3)*T_32)+
    mole_fractions[6] / ((2.8180810494342415e-05 -0.0006105383744688015*ln_T + 0.004973487092905404*ln_T_2 -0.011253853297732248*ln_T_3)*T_32)+
    mole_fractions[7] / ((2.816843323426154e-05 -0.0006102702206393846*ln_T + 0.00497130269358626*ln_T_2 -0.011248910506281885*ln_T_3)*T_32)+
    mole_fractions[8] / ((2.5950667401864423e-05 -0.0005591150913811907*ln_T + 0.004549445811591998*ln_T_2 -0.010193756040796769*ln_T_3)*T_32));
    _[4] = (1. - mass_fractions[4]) * mole_fractions[4] / (
    mole_fractions[0] / ((9.371675695604214e-06 -0.0001768414201939315*ln_T + 0.0015097675313645394*ln_T_2 -0.0024820837114003943*ln_T_3)*T_32)+
    mole_fractions[1] / ((5.569865822310673e-06 -0.00011630889857421319*ln_T + 0.0009446644908807677*ln_T_2 -0.0019946621508287895*ln_T_3)*T_32)+
    mole_fractions[2] / ((9.591588708458134e-06 -0.00022431069595172284*ln_T + 0.001951005759173615*ln_T_2 -0.004879738156005005*ln_T_3)*T_32)+
    mole_fractions[3] / ((3.476688833043835e-05 -0.0007402656486890561*ln_T + 0.006009838481386958*ln_T_2 -0.013182221813086606*ln_T_3)*T_32)+
    mole_fractions[4] / ((7.178862148012006e-06 -0.0001466802164039969*ln_T + 0.0011983390493111078*ln_T_2 -0.0024139746041368644*ln_T_3)*T_32)+
    mole_fractions[5] / ((7.071689894186152e-06 -0.00014449044746017698*ln_T + 0.0011804492090948923*ln_T_2 -0.0023779367065328047*ln_T_3)*T_32)+
    mole_fractions[6] / ((5.541442774940739e-06 -0.00011571537380374643*ln_T + 0.0009398438642391872*ln_T_2 -0.001984483381966367*ln_T_3)*T_32)+
    mole_fractions[7] / ((5.514570624715073e-06 -0.00011515423457799546*ln_T + 0.0009352862739988427*ln_T_2 -0.0019748600153222692*ln_T_3)*T_32)+
    mole_fractions[8] / ((5.189634899189305e-06 -0.0001076333130435482*ln_T + 0.000875505836327247*ln_T_2 -0.0018223379695537362*ln_T_3)*T_32));
    _[5] = (1. - mass_fractions[5]) * mole_fractions[5] / (
    mole_fractions[0] / ((9.340544322994412e-06 -0.0001762539781692879*ln_T + 0.001504752298539638*ln_T_2 -0.002473838582633779*ln_T_3)*T_32)+
    mole_fractions[1] / ((5.45871533159476e-06 -0.00011398787477156308*ln_T + 0.0009258130633826859*ln_T_2 -0.0019548573002363357*ln_T_3)*T_32)+
    mole_fractions[2] / ((9.43984208546376e-06 -0.00022076192091071393*ln_T + 0.0019201392839320014*ln_T_2 -0.004802536786265476*ln_T_3)*T_32)+
    mole_fractions[3] / ((3.4705768325682896e-05 -0.0007389642656162268*ln_T + 0.005999273217303768*ln_T_2 -0.013159047537257077*ln_T_3)*T_32)+
    mole_fractions[4] / ((7.071689894186152e-06 -0.00014449044746017698*ln_T + 0.0011804492090948923*ln_T_2 -0.0023779367065328047*ln_T_3)*T_32)+
    mole_fractions[5] / ((6.9628682436852225e-06 -0.00014226697765175724*ln_T + 0.0011622840416188854*ln_T_2 -0.0023413441804091264*ln_T_3)*T_32)+
    mole_fractions[6] / ((5.429710473016042e-06 -0.00011338220072802005*ln_T + 0.0009208937599674065*ln_T_2 -0.0019444701750446754*ln_T_3)*T_32)+
    mole_fractions[7] / ((5.402282570155895e-06 -0.00011280945637947845*ln_T + 0.0009162419125588868*ln_T_2 -0.0019346477840828354*ln_T_3)*T_32)+
    mole_fractions[8] / ((5.090805654457106e-06 -0.00010558358907591895*ln_T + 0.0008588330679642043*ln_T_2 -0.0017876341245480194*ln_T_3)*T_32));
    _[6] = (1. - mass_fractions[6]) * mole_fractions[6] / (
    mole_fractions[0] / ((8.616263166230654e-06 -0.0001687602602849943*ln_T + 0.0014069470958023526*ln_T_2 -0.002558973454509375*ln_T_3)*T_32)+
    mole_fractions[1] / ((4.0388489733899366e-06 -8.596612107899355e-05*ln_T + 0.0006978883091977138*ln_T_2 -0.0015297813741923412*ln_T_3)*T_32)+
    mole_fractions[2] / ((5.2249805662217574e-06 -0.0001275297292116438*ln_T + 0.0011696509910805663*ln_T_2 -0.003037094033578458*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.8180810494342415e-05 -0.0006105383744688015*ln_T + 0.004973487092905404*ln_T_2 -0.011253853297732248*ln_T_3)*T_32)+
    mole_fractions[4] / ((5.541442774940739e-06 -0.00011571537380374643*ln_T + 0.0009398438642391872*ln_T_2 -0.001984483381966367*ln_T_3)*T_32)+
    mole_fractions[5] / ((5.429710473016042e-06 -0.00011338220072802005*ln_T + 0.0009208937599674065*ln_T_2 -0.0019444701750446754*ln_T_3)*T_32)+
    mole_fractions[6] / ((4.007411939249679e-06 -8.529698987327842e-05*ln_T + 0.0006924561826817481*ln_T_2 -0.0015178740734726416*ln_T_3)*T_32)+
    mole_fractions[7] / ((3.977611366124012e-06 -8.466269042449412e-05*ln_T + 0.0006873068265833038*ln_T_2 -0.001506586609641165*ln_T_3)*T_32)+
    mole_fractions[8] / ((3.86259571337696e-06 -8.171465961132411e-05*ln_T + 0.0006632389396733325*ln_T_2 -0.0014370620481812637*ln_T_3)*T_32));
    _[7] = (1. - mass_fractions[7]) * mole_fractions[7] / (
    mole_fractions[0] / ((8.608910807271654e-06 -0.00016861625516494404*ln_T + 0.0014057465312553913*ln_T_2 -0.002556789852286297*ln_T_3)*T_32)+
    mole_fractions[1] / ((4.009282075250514e-06 -8.533679535720242e-05*ln_T + 0.0006927793307024408*ln_T_2 -0.001518582418657164*ln_T_3)*T_32)+
    mole_fractions[2] / ((5.1975721852789544e-06 -0.00012686075382399255*ln_T + 0.0011635154199473926*ln_T_2 -0.003021162523561163*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.816843323426154e-05 -0.0006102702206393846*ln_T + 0.00497130269358626*ln_T_2 -0.011248910506281885*ln_T_3)*T_32)+
    mole_fractions[4] / ((5.514570624715073e-06 -0.00011515423457799546*ln_T + 0.0009352862739988427*ln_T_2 -0.0019748600153222692*ln_T_3)*T_32)+
    mole_fractions[5] / ((5.402282570155895e-06 -0.00011280945637947845*ln_T + 0.0009162419125588868*ln_T_2 -0.0019346477840828354*ln_T_3)*T_32)+
    mole_fractions[6] / ((3.977611366124012e-06 -8.466269042449412e-05*ln_T + 0.0006873068265833038*ln_T_2 -0.001506586609641165*ln_T_3)*T_32)+
    mole_fractions[7] / ((3.947585833011331e-06 -8.40236027457935e-05*ln_T + 0.0006821185987800866*ln_T_2 -0.0014952139384646593*ln_T_3)*T_32)+
    mole_fractions[8] / ((3.836230291913693e-06 -8.11568892412802e-05*ln_T + 0.000658711783462094*ln_T_2 -0.0014272529070283508*ln_T_3)*T_32));
    _[8] = (1. - mass_fractions[8]) * mole_fractions[8] / (
    mole_fractions[0] / ((7.753386136651787e-06 -0.00015012001517394324*ln_T + 0.001260489007680266*ln_T_2 -0.0022244441293002544*ln_T_3)*T_32)+
    mole_fractions[1] / ((3.8904263458668995e-06 -8.230342706964434e-05*ln_T + 0.0006680176844741743*ln_T_2 -0.001447416314766726*ln_T_3)*T_32)+
    mole_fractions[2] / ((5.5498678214477205e-06 -0.00013317697631526514*ln_T + 0.0011956571686973264*ln_T_2 -0.003062175134485725*ln_T_3)*T_32)+
    mole_fractions[3] / ((2.5950667401864423e-05 -0.0005591150913811907*ln_T + 0.004549445811591998*ln_T_2 -0.010193756040796769*ln_T_3)*T_32)+
    mole_fractions[4] / ((5.189634899189305e-06 -0.0001076333130435482*ln_T + 0.000875505836327247*ln_T_2 -0.0018223379695537362*ln_T_3)*T_32)+
    mole_fractions[5] / ((5.090805654457106e-06 -0.00010558358907591895*ln_T + 0.0008588330679642043*ln_T_2 -0.0017876341245480194*ln_T_3)*T_32)+
    mole_fractions[6] / ((3.86259571337696e-06 -8.171465961132411e-05*ln_T + 0.0006632389396733325*ln_T_2 -0.0014370620481812637*ln_T_3)*T_32)+
    mole_fractions[7] / ((3.836230291913693e-06 -8.11568892412802e-05*ln_T + 0.000658711783462094*ln_T_2 -0.0014272529070283508*ln_T_3)*T_32)+
    mole_fractions[8] / ((3.711965635388281e-06 -7.804341552948402e-05*ln_T + 0.0006334890212170214*ln_T_2 -0.0013561200989773644*ln_T_3)*T_32));
}

void fg_rates(const float log_T, const float T, const float T2, const float T4, const float rcp_T, const float rcp_T2, const float P0_RT, const float rcp_P0_RT, const float exp_Gibbs0_RT[], const float concentrations[], float* molar_rates) {
 float c, Pr, logFcent, logPr_c, f1;
    c = exp2(-12050.746610262851 * rcp_T -0.406 * log_T + 31.723952184259364);
    const float cR0 = c * concentrations[1]*concentrations[3];
    c = exp2(-4566.491727125329 * rcp_T + 2.67 * log_T -4.299027692777283);
    const float cR1 = c * concentrations[0]*concentrations[4];
    c = exp2(-2490.1536763179456 * rcp_T + 1.51 * log_T + 7.754887502163468);
    const float cR2 = c * concentrations[0]*concentrations[5];
    c = exp2(-9728.297161125503 * rcp_T + 2.02 * log_T + 1.570462931026041);
    const float cR3 = c * concentrations[2]*concentrations[4];
    c = exp2(-75779.07893121493 * rcp_T -1.4 * log_T + 45.37946752547428) * (2.5*concentrations[0]+concentrations[1]+12.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
    const float cR4 = c * concentrations[0];
    c = exp2(-0.0 * rcp_T -0.5 * log_T + 12.589885179290201) * (2.5*concentrations[0]+concentrations[1]+12.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
    const float cR5 = c * concentrations[4]*concentrations[4];
    c = exp2(-0.0 * rcp_T -1.0 * log_T + 22.168520327912255) * (2.5*concentrations[0]+concentrations[1]+12.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
    const float cR6 = c * concentrations[3]*concentrations[4];
    c = exp2(-0.0 * rcp_T -2.0 * log_T + 35.14528036742985) * (2.5*concentrations[0]+concentrations[1]+12.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
    const float cR7 = c * concentrations[3]*concentrations[5];
    
        Pr = exp2(-381.0007723999002 * rcp_T -1.72 * log_T + 29.245811916072732) * (2.0*concentrations[0]+0.78*concentrations[1]+11.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
        logFcent = log2(0.19999999999999996 * exp2(-1.4426950408889633e+30*T) + 0.8 * exp2(-1.4426950408889633e-30*T) + exp2(0.0*rcp_T));
        logPr_c = log2(Pr) - 0.67*logFcent - 1.328771237954945;
        f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.4914460711655217);
        c = Pr / (exp2(0.0 * rcp_T - 0.6 * log_T - 20.492283523798655) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const float cR8 = c * concentrations[1]*concentrations[3];
    c = exp2(-597.4916838512156 * rcp_T + 0.0 * log_T + 23.984679905783736);
    const float cR9 = c * concentrations[3]*concentrations[6];
    c = exp2(-214.16773600985246 * rcp_T + 0.0 * log_T + 26.07704223964188);
    const float cR10 = c * concentrations[3]*concentrations[6];
    c = exp2(-0.0 * rcp_T + 0.0 * log_T + 24.95393638235263);
    const float cR11 = c * concentrations[4]*concentrations[6];
    c = exp2(360.8181857521921 * rcp_T + 0.0 * log_T + 24.78456615693749);
    const float cR12 = c * concentrations[5]*concentrations[6];
    c = exp2(-8698.840043627297 * rcp_T + 0.0 * log_T + 28.645814086990296);
    const float cR13 = c * concentrations[6]*concentrations[6];
    c = exp2(1182.859295867297 * rcp_T + 0.0 * log_T + 16.98815209769054);
    const float cR14 = c * concentrations[6]*concentrations[6];
    
        Pr = exp2(-33032.65080829928 * rcp_T + 0.0 * log_T + 36.80664593981008) * (2.5*concentrations[0]+concentrations[1]+12.0*concentrations[2]+concentrations[3]+concentrations[4]+concentrations[5]+concentrations[6]+concentrations[7]+concentrations[8]);
        logFcent = log2(0.5 * exp2(-1.4426950408889633e+30*T) + 0.5 * exp2(-1.4426950408889633e-30*T) + exp2(0.0*rcp_T));
        logPr_c = log2(Pr) - 0.67*logFcent - 1.328771237954945;
        f1 = logPr_c / (-0.14*logPr_c-1.27*logFcent-2.4914460711655217);
        c = Pr / (exp2(35159.80832188866 * rcp_T - 0.0 * log_T - 48.06819724919299) * Pr + 1.) * exp2(logFcent/(f1*f1+1.));
    const float cR15 = c * concentrations[7];
    c = exp2(-2882.189532064794 * rcp_T + 0.0 * log_T + 24.522529810666775);
    const float cR16 = c * concentrations[3]*concentrations[7];
    c = exp2(-5771.63898738416 * rcp_T + 0.0 * log_T + 25.522529810666775);
    const float cR17 = c * concentrations[3]*concentrations[7];
    c = exp2(-2882.189532064794 * rcp_T + 2.0 * log_T + 3.255500733148386);
    const float cR18 = c * concentrations[4]*concentrations[7];
    c = exp2(-0.0 * rcp_T + 0.0 * log_T + 19.931568569324174);
    const float cR19 = c * concentrations[5]*concentrations[7];
    c = exp2(-6938.308654393763 * rcp_T + 0.0 * log_T + 29.11147765933911);
    const float cR20 = c * concentrations[5]*concentrations[7];
    molar_rates[0] = -cR1-cR2-cR4+cR9+cR17;
    molar_rates[1] = -cR0+cR5-cR8+cR9+cR11+cR12+cR13+cR14;
    molar_rates[2] = cR2-cR3+cR7+cR12+cR16+cR19+cR20;
    molar_rates[3] = -cR0+cR1+cR2+2*cR4-cR6-cR7-cR8-cR9-cR10-cR16-cR17;
    molar_rates[4] = cR0-cR1-cR3-2*cR5-cR6-cR11-cR18;
    molar_rates[5] = cR0+cR1-cR2+2*cR3+cR6-cR7+2*cR10+cR11-cR12+2*cR15+cR16+cR18-cR19-cR20;
    molar_rates[6] = cR8-cR9-cR10-cR11-cR12-2*cR13-2*cR14+cR17+cR18+cR19+cR20;
    molar_rates[7] = cR13+cR14-cR15-cR16-cR17-cR18-cR19-cR20;
}

