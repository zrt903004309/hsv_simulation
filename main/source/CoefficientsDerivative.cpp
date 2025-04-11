#include <CoefficientsDerivative.h>

void getCoefficientsDerivative(
	const double& Mach,
	const double& Alpha,
	const double& Beta,
	const double& Delta_e,
	const double& Delta_a,
	const double& Delta_r,
	const double& C,
	const double& B,
	const double& V,
	const double& p,
	const double& q,
	const double& r,
	vector<vector<double>>& ans) {

	// 所有角度值是弧度制
	double RE = Delta_a;		// Right Elevon 右襟翼
	double LE = Delta_e;		// Left Elevon 左襟翼
	double RUD = Delta_r;

	double CLbv, CL_dRE, CL_dLE, CDbv, CD_dRE, CD_dLE, CD_dRUD, CYbv, CY_dRE, CY_dLE, CY_dRUD;
	double Clbv, Cl_dRE, Cl_dLE, Cl_dRUD, Clr, Clp, Cl_state, Cnbv, Cn_dRE, Cn_dLE, Cn_dRUD, Cnp, Cnr, Cn_state, Cmbv, Cm_dRE, Cm_dLE, Cm_dRUD, Cm_q, Cm_state;

	if (Mach < 0 || Mach > 30) {
		std::cout << Delta_a << " " << Mach << "坠毁了" << endl;
		exit(1);
	}
	else if (Mach <= 1.25) {
		CLbv = -5.2491E-04 + Alpha * 1.5746E-02 + (Alpha * Mach) * 6.0213E-03
			- 3.4437E-04 * pow(Alpha, 2) + pow(Alpha * Mach, 2) * 1.4471E-04
			- 5.1952E-05 * pow(Alpha, 3) + 3.4771E-05 * pow(Alpha, 4)
			+ 2.7717E-03 * pow(Mach, 4) - 2.3034E-06 * pow(Alpha, 5);

		CL_dRE = 1.406E-04 * Alpha - 8.584E-04 * Mach + 8.879E-05 * (Alpha * Mach)
			- 1.703E-06 * 2 * pow(Mach, 2) * RE - 3.727E-05 * 3 * pow(RE, 2) 
			+ 1.781E-07 * 2 * RE + 7.912E-07 * 2 * pow(Alpha * Mach, 2) * RE 
			+ 2.465E-08 * 2 * pow(Alpha, 2) * RE - 5.942E-09 * 3 * pow(Alpha * Mach, 3) * pow(RE, 2)
			+ 2.672E-08 * 4 * pow(RE, 3) - 1.610E-11 * 4 * pow(Alpha * Mach, 4) * pow(RE, 3)
			+ 7.624E-08 * RE * pow(RE, 4) + 1.388E-13 * 5 * pow(Alpha * Mach, 5) * pow(RE, 4);
		
		CL_dLE = 1.406E-04 * Alpha - 8.584E-04 * Mach + 8.879E-05 * (Alpha * Mach)
			- 1.703E-06 * 2 * pow(Mach, 2) * LE - 3.727E-05 * 3 * pow(LE, 2) 
			+ 1.781E-07 * 2 * LE + 7.912E-07 * 2 * pow(Alpha * Mach, 2) * LE 
			+ 2.465E-08 * 2 * pow(Alpha, 2) * LE - 5.942E-09 * 3 * pow(Alpha * Mach, 3) * pow(LE, 2)
			+ 2.672E-08 * 4 * pow(RE, 3) - 1.610E-11 * 4 * pow(Alpha * Mach, 4) * pow(RE, 3)
			+ 7.624E-08 * 5 * pow(LE, 4) + 1.388E-13 * 5 * pow(Alpha * Mach, 5) * pow(RE, 4);

		CDbv = +1.1457E-02 + CLbv * (-2.4645E-02) + Mach * (0)
			+ (CLbv * Mach) * (4.9698E-02) + pow(CLbv, 2) * (-1.9112e+0)
			+ pow(Mach, 2) * (0) + pow(CLbv * Mach, 2) * (3.5404e+0)
			+ pow(CLbv, 3) * (4.4334e+01) + pow(Mach, 3) * (0)
			+ pow(CLbv * Mach, 3) * (-7.0367e+01)
			+ pow(CLbv, 4) * (-2.3841e+02) + pow(Mach, 4) * (0)
			+ pow(CLbv * Mach, 4) * (4.1750e+02) + pow(CLbv, 5) * (4.1734e+02)
			+ pow(Mach, 5) * (5.4910E-02)
			+ pow(CLbv * Mach, 5) * (-7.9055e+02);

		CD_dLE = 3.38E-07 * Alpha - 2.79E-04 * Mach - 1.53E-04 * Alpha * Mach
			+ 9.39E-08 * 2 * LE - 5.69E-07 * pow(Alpha * Mach, 2) * 2 * LE
			+ 4.14E-07 * pow(Alpha, 2) * 2 * LE - 1.68E-05 * pow(Mach, 2) * 2 * LE
			- 1.84E-06 * 3 * pow(LE, 2) + 5.76E-08 * 4 * pow(LE, 3) + 5.71E-09 * 5 * pow(LE, 4)
			- 8.93E-15 * pow(Alpha * Mach, 5) * 5 * pow(LE, 4)
			- 7.58E-12 * pow(Alpha * Mach, 4) * 4 * pow(LE, 3)
			- 3.94E-10 * pow(Alpha * Mach, 3) * 3 * pow(LE, 2);

		CD_dRE = 3.38E-07 * Alpha - 2.79E-04 * Mach	- 1.53E-04 * Alpha * Mach
			+ 9.39E-08 * 2 * RE	- 5.69E-07 * pow(Alpha * Mach, 2) * 2 * RE 
			+ 4.14E-07 * pow(Alpha, 2) * 2 * RE - 1.68E-05 * pow(Mach, 2) * 2 * RE
			- 1.84E-06 * 3 * pow(RE, 2) + 5.76E-08 * 4 * pow(RE, 3)	+ 5.71E-09 * 5 * pow(RE, 4) 
			- 8.93E-15 * pow(Alpha * Mach , 5) * 5 * pow(RE, 4)
			- 7.58E-12 * pow(Alpha * Mach , 4) * 4 * pow(RE, 3)
			- 3.94E-10 * pow(Alpha * Mach , 3) * 3 * pow(RE, 2);

		CD_dRUD = 5.03E-06 * 2 * RUD - 1.30E-07 * 2 * pow(Alpha * Mach, 2) * RUD 
			- 3.50E-08 * pow(Alpha, 2) * 2 * RUD + 4.53E-06 * pow(Mach, 2) * 2 * RUD
			+ 7.54E-09 * 4 * pow(RUD, 3) + 3.12E-12 * pow(Alpha * Mach, 4) * 4 * pow(RUD, 3);

		CYbv = -4.750E-01 - 5.000E-02 * Mach;

		CY_dLE = -(- 2.13E-07 * Alpha + 1.990E-05 * Mach + 6.17E-08 * (Alpha * Mach)
			+ 1.37E-07 * 2 * LE - 3.40E-07 * 3 * pow(LE, 2) + 3.53E-09 * pow(Alpha * Mach, 2) * 2 * LE
			- 2.66E-09 * pow(Alpha, 2) * 2 * LE	+ 3.92E-08 * pow(Mach, 2) * 2 * LE 
			+ 5.42E-11 * pow(Alpha * Mach, 3) * 3 * pow(LE, 2) - 4.73E-10 * 4 * pow(LE, 3)
			+ 7.35E-14 * pow(Alpha * Mach, 4) * 4 * pow(LE, 3) + 6.53E-10 * 5 * pow(LE, 4)
			- 1.11E-15 * pow(Alpha * Mach, 5) * 5 * pow(LE, 4));

		CY_dRE = -2.13E-07 * Alpha + 1.990E-05 * Mach + 6.17E-08 * (Alpha * Mach)
			+ 1.37E-07 * 2 * RE - 3.40E-07 * 3 * pow(RE, 2) + 3.53E-09 * pow(Alpha * Mach, 2) * 2 * RE
			- 2.66E-09 * pow(Alpha, 2) * 2 * RE + 3.92E-08 * pow(Mach, 2) * 2 * RE
			+ 5.42E-11 * pow(Alpha * Mach, 3) * 3 * pow(RE, 2) - 4.73E-10 * 4 * pow(RE, 3)
			+ 7.35E-14 * pow(Alpha * Mach, 4) * 4 * pow(RE, 3) + 6.53E-10 * 5 * pow(RE, 4)
			- 1.11E-15 * pow(Alpha * Mach, 5) * 5 * pow(RE, 4);

		CY_dRUD = 2.440E-03;

		Clbv = -9.380E-02 - 1.250E-02 * Mach;

		Cl_dLE = -(3.690E-05 * Alpha + 1.926E-04 * Mach - 8.500E-06 * (Alpha * Mach)
			+ 3.762E-06 * 2 * LE - 5.302E-08 * pow(Alpha * Mach, 2) * 2 * LE 
			+ 2.100E-06 * pow(Mach, 2) * 2 * LE + 8.400E-06 * 3 * pow(LE, 2)
			+ 1.153E-09 * pow(Alpha * Mach, 3) * 3 * pow(LE, 2) - 3.576E-08 * pow(Alpha, 2) * 2 * LE
			- 1.137E-08 * 4 * pow(LE, 3) + 1.011E-12 * pow(Alpha * Mach, 4) * 4 * pow(LE, 3) 
			- 1.676E-08 * 5 * pow(LE, 4) - 2.984E-14 * pow(Alpha * Mach, 5) * 5 * pow(LE, 4));

		Cl_dRE = 3.690E-05 * Alpha + 1.926E-04 * Mach - 8.500E-06 * (Alpha * Mach)
			+ 3.762E-06 * 2 * RE - 5.302E-08 * pow(Alpha * Mach, 2) * 2 * RE
			+ 2.100E-06 * pow(Mach, 2) * 2 * RE + 8.400E-06 * 3 * pow(RE, 2)
			+ 1.153E-09 * pow(Alpha * Mach, 3) * 3 * pow(RE, 2) - 3.576E-08 * pow(Alpha, 2) * 2 * RE
			- 1.137E-08 * 4 * pow(RE, 3) + 1.011E-12 * pow(Alpha * Mach, 4) * 4 * pow(RE, 3)
			- 1.676E-08 * 5 * pow(RE, 4) - 2.984E-14 * pow(Alpha * Mach, 5) * 5 * pow(RE, 4);

		Cl_dRUD = 7.0E-04;
		Clr = +2.6250E-01 + 2.50E-02 * (Mach);
		Clp = -1.33750E-01 - 1.250E-02 * (Mach);

		Cnbv = +1.062E-01 + 6.250E-02 * Mach;

		Cn_dLE = -(-2.7E-7 * Alpha - 1.008E-05 * Mach + 3.564E-07 * (Alpha * Mach)
			+ 1.1E-7 * 3 * pow(LE, 2) + 1.11E-07 * 3 * pow(LE, 2)
			- 9.32E-12 * pow(Alpha * Mach, 3) * 3 * pow(LE, 2) + 2.89E-25 * 4 * pow(LE, 3)
			+ 1.82E-28 * pow(Alpha * Mach, 4) * 4 * pow(LE, 3)
			- 2.2046E-10 * 5 * pow(LE, 4) + 2.22E-16 * pow(Alpha * Mach, 5) * 5 * pow(LE, 4));

		Cn_dRE = -2.7E-7 * Alpha - 1.008E-05 * Mach + 3.564E-07 * (Alpha * Mach)
			+ 1.1E-7 * 3 * pow(RE, 2) + 1.11E-07 * 3 * pow(RE, 2)
			- 9.32E-12 * pow(Alpha * Mach, 3) * 3 * pow(RE, 2) + 2.89E-25 * 4 * pow(RE, 3)
			+ 1.82E-28 * pow(Alpha * Mach, 4) * 4 * pow(RE, 3)
			- 2.2046E-10 * 5 * pow(RE, 4) + 2.22E-16 * pow(Alpha * Mach, 5) * 5 * pow(RE, 4);

		Cn_dRUD = -3.000E-03;
		Cnp = +1.790E-01 + 2.0E-02 * Mach;
		Cnr = -1.2787 - 1.375E-001 * Mach;

		Cmbv = +(-1.8316E-003) + CLbv * (-1.0306E-001) + Mach * (0)
			+ (CLbv * Mach) * (-1.8335E-001) + pow(CLbv, 2) * (-1.1839e+000)
			+ pow(Mach, 2) * (-2.8113E-03)
			+ pow(CLbv * Mach, 2) * (-1.3362e+00) + pow(CLbv, 3) * (9.0641e+00)
			+ pow(Mach, 3) * (0) + pow(CLbv * Mach, 3) * (2.6964e+001)
			+ pow(CLbv, 4) * (-6.3590e+01) + pow(Mach, 4) * (0)
			+ pow(CLbv * Mach, 4) * (-8.0921e+01)
			+ pow(CLbv, 5) * (1.6885e+02) + pow(Mach, 5) * (0)
			+ pow(CLbv * Mach, 5) * (-4.2209e+00);

		Cm_dRE = 4.55E-05 * Alpha + 6.665E-04 * Mach - 2.770E-05 * (Alpha * Mach)
			- 1.600E-06 * 2 * RE - 1.000E-07 * pow(Alpha * Mach, 2) * 2 * RE
			+ 2.300E-06 * pow(Mach, 2) * 2 * RE + 1.920E-05 * 3 * pow(RE, 2)
			+ 1.90E-09 * pow(Alpha * Mach, 3) * 3 * pow(RE, 2) - 4.69E-10 * 4 * pow(RE, 3)
			+ 1.29E-12 * pow(Alpha * Mach, 4) * 4 * pow(RE, 3) - 3.87E-08 * 5 * pow(RE, 4)
			- 4.67E-14 * pow(Alpha * Mach, 5) * 5 * pow(RE, 4);

		Cm_dLE = 4.55E-05 * Alpha + 6.665E-04 * Mach - 2.770E-05 * (Alpha * Mach)
			- 1.600E-06 * 2 * LE - 1.000E-07 * pow(Alpha * Mach, 2) * 2 * LE
			+ 2.300E-06 * pow(Mach, 2) * 2 * LE + 1.920E-05 * 3 * pow(LE, 2)
			+ 1.90E-09 * pow(Alpha * Mach, 3) * 3 * pow(LE, 2) - 4.69E-10 * 4 * pow(LE, 3)
			+ 1.29E-12 * pow(Alpha * Mach, 4) * 4 * pow(LE, 3) - 3.87E-08 * 5 * pow(LE, 4)
			- 4.67E-14 * pow(Alpha * Mach, 5) * 5 * pow(LE, 4);

		Cm_dRUD = - 1.0E-07 + 5.8E-06 * 2 * RUD + 6.482E-09 * pow(Alpha * Mach, 2) * 2 * RUD;

		Cm_q = -1.0313 - 3.125E-01 * Mach;
	}
	else if (Mach < 4) {
		CLbv = +1.9920E-001 + Mach * (2.3402E-001) + Alpha * (3.8202E-002)
			+ (Alpha * Mach) * (-2.4626E-003) + pow(Mach, 2) * (-6.4872E-001)
			+ pow(Alpha, 2) * (-6.9523E-003)
			+ pow(Alpha * Mach * Mach, 2) * (4.5735E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (2.1241E-007)
			+ pow(Alpha * Mach, 2) * (-1.0521E-004)
			+ pow(Alpha * Mach, 4) * (-9.5825E-009)
			+ pow(Mach, 3) * (3.9121E-001)
			+ pow(Alpha, 3) * (1.0295E-003) + pow(Mach, 4) * (-9.1356E-002)
			+ pow(Alpha, 4) * (-5.7398E-005) + pow(Mach, 5) * (7.4089E-003)
			+ pow(Alpha, 5) * (1.0934E-006);

		CL_dRE = Alpha * (-3.3093E-5) + Mach * (-1.4287E-4) + (Alpha * Mach) * 6.1071E-6
			+ 2 * RE * 2.7242E-4 + pow(Alpha * Mach, 2) * 2 * RE * (-9.1890E-8)
			+ pow(Alpha, 2) * 2 * RE * (3.4060E-7) + pow(Mach, 2) * 2 * RE * (-6.3863E-6)
			+ 3 * pow(RE, 2) * (3.8067E-6) + pow(Alpha * Mach, 3) *3 * pow(RE, 2) * (2.3165E-011)
			+ 4 * pow(RE, 3) * (-3.7716E-007) + pow(Alpha * Mach, 4) * 4 * pow(RE, 3) * (7.9006E-14)
			+ 5 * pow(RE, 4) * (-8.5345E-009) + pow(Alpha * Mach, 5) * 5 * pow(RE, 4) * (-2.5698E-017);

		CL_dLE = Alpha * (-3.3093E-5) + Mach * (-1.4287E-4) + (Alpha * Mach) * 6.1071E-6
			+ 2 * LE * 2.7242E-4 + pow(Alpha * Mach, 2) * 2 * LE * (-9.1890E-8)
			+ pow(Alpha, 2) * 2 * LE * (3.4060E-7) + pow(Mach, 2) * 2 * LE * (-6.3863E-6)
			+ 3 * pow(LE, 2) * (3.8067E-6) + pow(Alpha * Mach, 3) * 3 * pow(LE, 2) * (2.3165E-011)
			+ 4 * pow(LE, 3) * (-3.7716E-007) + pow(Alpha * Mach, 4) * 4 * pow(LE, 3) * (7.9006E-14)
			+ 5 * pow(LE, 4) * (-8.5345E-009) + pow(Alpha * Mach, 5) * 5 * pow(LE, 4) * (-2.5698E-017);

		CDbv = +(-8.2073E-002) + CLbv * (-9.1273E-002)
			+ Mach * (2.1845E-001)
			+ (CLbv * Mach) * (3.2202E-002) + pow(CLbv, 2) * (1.6325e+000)
			+ pow(Mach, 2) * (-1.3680E-001)
			+ pow(CLbv * Mach, 2) * (5.7526E-002)
			+ pow(CLbv, 3) * (-1.1575e+000) + pow(Mach, 3) * (3.8791E-002)
			+ pow(CLbv * Mach, 3) * (-2.4002E-001)
			+ pow(CLbv, 4) * (-8.5306e+000)
			+ pow(Mach, 4) * (-5.2527E-003)
			+ pow(CLbv * Mach, 4) * (3.5543E-001)
			+ pow(CLbv, 5) * (1.7259e+001) + pow(Mach, 5) * (2.7435E-004)
			+ pow(CLbv * Mach, 5) * (-1.4983E-001);

		CD_dRE = Alpha * (-3.6923E-005) + Mach * (1.3641E-007) + (Alpha * Mach) * (5.1142E-006)
			+ 2 * RE * (1.2125E-005) + 2 * RE * pow(Alpha * Mach, 2) * (3.5662E-009)
			+ 2 * RE * pow(Alpha, 2) * (-1.3848E-008) + 2 * RE * pow(Mach, 2) * (-3.3763E-007)
			+ 3 * pow(RE, 2) * (3.9119E-008) + 3 * pow(RE, 2) * pow(Alpha * Mach, 3) * (-9.7714E-013)
			+ 4 * pow(RE, 3) * (4.5137E-009) + 4 * pow(RE, 3) * pow(Alpha * Mach, 4) * (-6.6207E-016)
			+ 5 * pow(RE, 4) * (-1.1538E-010) + 5 * pow(RE, 4) * pow(Alpha * Mach, 5) * (4.1917E-019);

		CD_dLE = Alpha * (-3.6923E-005) + Mach * (1.3641E-007) + (Alpha * Mach) * (5.1142E-006)
			+ 2 * LE * (1.2125E-005) + 2 * LE * pow(Alpha * Mach, 2) * (3.5662E-009)
			+ 2 * LE * pow(Alpha, 2) * (-1.3848E-008) + 2 * LE * pow(Mach, 2) * (-3.3763E-007)
			+ 3 * pow(LE, 2) * (3.9119E-008) + 3 * pow(LE, 2) * pow(Alpha * Mach, 3) * (-9.7714E-013)
			+ 4 * pow(LE, 3) * (4.5137E-009) + 4 * pow(LE, 3) * pow(Alpha * Mach, 4) * (-6.6207E-016)
			+ 5 * pow(LE, 4) * (-1.1538E-010) + 5 * pow(LE, 4) * pow(Alpha * Mach, 5) * (4.1917E-019);

		CD_dRUD = Alpha * (2.6425E-021) + Mach * (1.8193E-020) + (Alpha * Mach) * (1.0319E-021)
			+ 2 * RUD * (8.7608E-006) + 2 * RUD * pow(Alpha * Mach, 2) * (5.4045E-010)
			+ 2 * RUD * pow(Alpha, 2) * (-2.8939E-008) + 2 * RUD * pow(Mach, 2) * (-2.9646E-007)
			+ 3 * pow(RUD, 2) * (-8.8556E-022) + 3 * pow(RUD, 2) * pow(Alpha * Mach, 3) * (-5.2022E-027)
			+ 4 * pow(RUD, 3) * (4.6754E-010) + 4 * pow(RUD, 3) * pow(Alpha * Mach, 4) * (2.6560E-016)
			+ 5 * pow(RUD, 4) * (1.5611E-024) + 5 * pow(RUD, 4) * pow(Alpha * Mach, 5) * (5.4442E-033);

		CYbv = +(0) + Mach * (0) + Alpha * (-1.1185E-002)
			+ (Alpha * Mach) * (3.0432E-003) + pow(Mach, 2) * (-3.7586E-001)
			+ pow(Alpha, 2) * (3.4004E-003)
			+ pow(Alpha * Mach * Mach, 2) * (-2.4047E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (3.6104E-007)
			+ pow(Alpha * Mach, 2) * (-8.7176E-005)
			+ pow(Alpha * Mach, 4) * (-5.3622E-010) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8160E-004) + pow(Mach, 4) * (9.4289E-002)
			+ pow(Alpha, 4) * (4.4848E-005) + pow(Mach, 5) * (-1.8384E-002)
			+ pow(Alpha, 5) * (-1.3021E-006);

		CY_dLE = -(2.27E-07	+ 4.11E-09 * (Alpha * Mach) - 5.04E-08 * 2 * LE
			+ 4.50E-14 * pow(Alpha * Mach, 2) * 2 * LE);

		CY_dRE = 2.27E-07 + 4.11E-09 * (Alpha * Mach) - 5.04E-08 * 2 * RE
			+ 4.50E-14 * pow(Alpha * Mach, 2) * 2 * RE;

		CY_dRUD = Alpha * (2.0067E-005)	+ Mach * (-5.7185E-004) + (Alpha * Mach) * (-1.5307E-005)
			+ 2 * RUD * (1.9243E-019) + 2 * RUD * pow(Alpha * Mach, 2) * (2.8011E-022)
			+ 2 * RUD * pow(Alpha, 2) * (-2.0404E-021) + 2 * RUD * pow(Mach, 2) * (-1.7950E-020)
			+ 3 * pow(RUD, 2) * (3.2768E-005) + 3 * pow(RUD, 2) * pow(Alpha * Mach, 3) * (1.2674E-012)
			+ 4 * pow(RUD, 3) * (7.7275E-023) + 4 * pow(RUD, 3) * pow(Alpha * Mach, 4) * (-3.2592E-029)
			+ 5 * pow(RUD, 4) * (-6.5825E-008) + 5 * pow(RUD, 4) * pow(Alpha * Mach, 5) * (1.2684E-017);

		Clbv = +(0) + Mach * (0) + Alpha * (5.9211E-004)
			+ (Alpha * Mach) * (-3.1579E-004) + pow(Mach, 2) * (-8.7296E-002)
			+ pow(Alpha, 2) * (-5.7398E-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1037E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (-6.8068E-008)
			+ pow(Alpha * Mach, 2) * (2.0549E-005)
			+ pow(Alpha * Mach, 4) * (3.6561E-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-2.8226E-016) + pow(Mach, 4) * (2.0334E-002)
			+ pow(Alpha, 4) * (1.9013E-007) + pow(Mach, 5) * (-3.7733E-003)
			+ pow(Alpha, 5) * (-9.6648E-019);

		Cl_dLE = -(1.170E-04 + 2.794E-08 * (Alpha * Mach) - 1.160E-06 * 2 * LE
			- 4.641E-11 * pow(Alpha * Mach, 2) * 2 * LE);

		Cl_dRE = 1.170E-04 + 2.794E-08 * (Alpha * Mach) - 1.160E-06 * 2 * RE
			- 4.641E-11 * pow(Alpha * Mach, 2) * 2 * RE;

		Cl_dRUD = 0.00011441 - 2.6824E-06 * Alpha - 3.5496E-06 * Mach
			+ 5.5547E-08 * (Alpha * Mach);

		Clr = +3.82E-01 - 1.06E-01 * Mach
			+ 1.94E-03 * Alpha - 8.15E-05 * (Alpha * Mach)
			+ 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2)
			- 1.02E-03 * pow(Mach, 3) - 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4)
			+ 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		Clp = +(0) + Mach * (0) + Alpha * (-1.2668E-005)
			+ (Alpha * Mach) * (1.7282E-005) + pow(Mach, 2) * (-1.0966E-001)
			+ pow(Alpha, 2) * (1.0751E-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.0989E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (6.1850E-009)
			+ pow(Alpha * Mach, 2) * (8.6481E-006)
			+ pow(Alpha * Mach, 4) * (-4.3707E-010)
			+ pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-1.1567E-005) + pow(Mach, 4) * (2.6725E-002)
			+ pow(Alpha, 4) * (1.5082E-006) + pow(Mach, 5) * (-5.0800E-003)
			+ pow(Alpha, 5) * (-6.1276E-008);

		Cnbv = +(0) + Mach * (0) + Alpha * (-2.3745E-003)
			+ (Alpha * Mach) * (8.5307E-004)
			+ pow(Mach, 2) * (1.4474E-001)
			+ pow(Alpha, 2) * (5.3105E-004)
			+ pow(Alpha * Mach * Mach, 2) * (-8.3462E-007)
			+ pow(Alpha * Alpha * Mach, 2) * (1.3335E-007)
			+ pow(Alpha * Mach, 2) * (-2.7081E-005)
			+ pow(Alpha * Mach, 4) * (-1.3450E-009)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.1046E-005)
			+ pow(Mach, 4) * (-3.9519E-002) + pow(Alpha, 4) * (-1.5141E-006)
			+ pow(Mach, 5) * (7.7646E-003) + pow(Alpha, 5) * (1.7278E-007);

		Cn_dRE = - 1.30E-05 - 8.93E-08 * (Alpha * Mach) + 1.97E-06 * 2 * RE
			+ 1.41E-11 * pow(Alpha * Mach, 2) * 2 * RE;

		Cn_dLE = -(- 1.30E-05 - 8.93E-08 * (Alpha * Mach) + 1.97E-06 * 2 * LE
			+ 1.41E-11 * pow(Alpha * Mach, 2) * 2 * LE);

		Cn_dRUD = - 5.28E-04 + 1.39E-05 * Alpha	+ 1.65E-05 * Mach - 3.13E-07 * Alpha * Mach;

		Cnp = +(1.7000E-001) + Alpha * (-6.4056E-018)
			+ Mach * (1.1333E-002) + (Alpha * Mach) * (2.3467E-018)
			+ pow(Alpha, 2) * (2.0917E-019)
			+ pow(Mach, 2) * (-5.3333E-003)
			+ pow(Alpha * Mach, 2) * (-5.0665E-020);

		Cnr = +(0) + Mach * (0) + Alpha * (-1.3332E-003)
			+ (Alpha * Mach) * (6.6899E-004)
			+ pow(Mach, 2) * (-1.0842e+000)
			+ pow(Alpha, 2) * (1.6434E-003)
			+ pow(Alpha * Mach * Mach, 2) * (-4.4258E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (1.2017E-007)
			+ pow(Alpha * Mach, 2) * (1.0819E-005)
			+ pow(Alpha * Mach, 4) * (-2.8899E-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8118E-004) + pow(Mach, 4) * (2.7379E-001)
			+ pow(Alpha, 4) * (6.7994E-005) + pow(Mach, 5) * (-5.2435E-002)
			+ pow(Alpha, 5) * (-2.5848E-006);

		Cmbv = +(-5.7643E-001) + Mach * (1.0553e+000) + CLbv * (-3.7951E-001)
			+ (CLbv * Mach) * (1.0483E-001) + pow(Mach, 2) * (-7.4344E-001)
			+ pow(CLbv, 2) * (-1.5412E-001)
			+ pow(CLbv * Mach * Mach, 2) * (-2.1133E-003)
			+ pow(CLbv * CLbv * Mach, 2) * (-1.7858E-001)
			+ pow(CLbv * Mach, 2) * (5.7805E-002)
			+ pow(CLbv * Mach, 4) * (-3.8875E-003)
			+ pow(Mach, 3) * (2.5341E-001)
			+ pow(CLbv, 3) * (-4.9731E-001) + pow(Mach, 4) * (-4.1938E-002)
			+ pow(CLbv, 4) * (7.1784e+000) + pow(Mach, 5) * (2.7017E-003)
			+ pow(CLbv, 5) * (-1.0331e+001);

		Cm_dRE = 2.89E-04 + 4.48E-06 * Alpha - 5.87E-06 * Mach + 9.72E-08 * Alpha * Mach;

		Cm_dLE = 2.89E-04 + 4.48E-06 * Alpha - 5.87E-06 * Mach + 9.72E-08 * Alpha * Mach;

		Cm_dRUD = 1.43E-07 * 4 * pow(RUD, 3) - 4.77E-22 * 5 * pow(RUD, 4) 
			- 3.38E-10 * 6 * pow(RUD, 5) + 2.63E-24 * 7 * pow(RUD, 6);

		Cm_q = +(0) + Mach * (0) + Alpha * (-1.0828E-002)
			+ (Alpha * Mach) * (4.2311E-003)
			+ pow(Mach, 2) * (-6.1171E-001)
			+ pow(Alpha, 2) * (4.6974E-003)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1593E-005)
			+ pow(Alpha * Alpha * Mach, 2) * (2.5378E-007)
			+ pow(Alpha * Mach, 2) * (-7.0964E-005)
			+ pow(Alpha * Mach, 4) * (4.1284E-008)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-1.1414E-003)
			+ pow(Mach, 4) * (1.5903E-001)
			+ pow(Alpha, 4) * (1.1176E-004) + pow(Mach, 5) * (-3.0665E-002)
			+ pow(Alpha, 5) * (-3.8123E-006);
	}
	else {
		// 1 lift force 升力 z轴
		CLbv = -8.19E-02 + 4.7E-02 * Mach + 1.86E-02 * Alpha - 4.73E-04 * Alpha * Mach - 9.19E-03 * pow(Mach, 2) - 1.52E-04 * pow(Alpha, 2) 
			+ 5.99E-07 * pow(Alpha * Mach, 2) + 7.74E-04 * pow(Mach, 3)	+ 4.08E-06 * pow(Alpha, 3) - 2.93E-05 * pow(Mach, 4)
			- 3.91E-07 * pow(Alpha, 4) + 4.12E-07 * pow(Mach, 5) + 1.3E-08 * pow(Alpha, 5);

		// 对右升降副翼的偏导
		CL_dRE = - 4.14E-04	- 3.51E-06 * Alpha + 8.72E-06 * Mach - 1.7E-07 * Alpha * Mach;

		// 对左升降副翼的偏导
		CL_dLE = - 4.14E-04 - 3.51E-06 * Alpha + 8.72E-06 * Mach - 1.7E-07 * Alpha * Mach;


		// 2 drag force 牵引力 x轴
		CDbv = 8.717E-02 - 3.307E-02 * Mach + 3.179E-03 * Alpha	- 1.25E-04 * Alpha * Mach + 5.036E-03 * pow(Mach, 2)
			- 1.1E-03 * pow(Alpha, 2) + 1.405E-07 * pow(Alpha * Mach, 2) - 3.658E-04 * pow(Mach, 3) + 3.175E-04 * pow(Alpha, 3)
			+ 1.274E-05 * pow(Mach, 4) - 2.985E-05 * pow(Alpha, 4) - 1.705E-07 * pow(Mach, 5) + 9.766E-07 * pow(Alpha, 5);

		// 对右升降副翼的偏导
		CD_dRE = - 3.6417E-05 - 5.3015E-07 * Alpha * Mach + 6.9629E-06 * 2 * RE + 2.1026E-12 * 2 * RE * pow(Alpha * Mach, 2);

		// 对左升降副翼的偏导
		CD_dLE = - 3.6417E-05 - 5.3015E-07 * Alpha * Mach + 6.9629E-06 * 2 * LE + 2.1026E-12 * 2 * LE * pow(Alpha * Mach, 2);

		// 对方向舵的偏导
		CD_dRUD = - 1.83E-06 + 9.13E-09 * Alpha * Mach + 1.97E-06 * 2 * RUD	- 1.77E-11 * pow(Alpha * Mach, 2) * 2 * RUD;


		// 3 side force 侧向力 y轴
		// 和侧滑角beta相关的系数
		CYbv = -2.9253E-01 * Mach + 2.8803E-03 * Alpha - 2.8943E-04 * (Alpha * Mach) + 5.4822E-02 * pow(Mach, 2)	+ pow(Alpha, 2) * (7.3535E-04)
			+ pow(Alpha * Mach * Mach, 2) * (-4.6490E-09) + pow(Alpha * Alpha * Mach, 2) * (-2.0675E-08) + pow(Alpha * Mach, 2) * (4.6205E-06)
			+ pow(Alpha * Mach, 4) * (2.6144E-11) + pow(Mach, 3) * (-4.3203E-03) + pow(Alpha, 3) * (-3.7405E-04) + pow(Mach, 4) * (1.5495E-04)
			+ pow(Alpha, 4) * (2.8183E-05) + pow(Mach, 5) * (-2.0829E-06) + pow(Alpha, 5) * (-5.2083E-07);

		// 对右升降副翼的偏导
		CY_dRE = 2.27E-07 + 4.11E-09 * Alpha * Mach	- 5.04E-08 * 2 * RE + 4.50E-14 * pow(Alpha * Mach, 2) * 2 * RE;

		// 对左升降副翼的偏导
		CY_dLE = -(2.27E-07	+ 4.11E-09 * Alpha * Mach - 5.04E-08 * 2 * LE + 4.50E-14 * pow(Alpha * Mach, 2) * 2 * LE);

		// 对方向舵的偏导
		CY_dRUD = 3.84E-04 - 1.17E-05 * Alpha - 1.07E-05 * Mach	+ 2.60E-07 * (Alpha * Mach);


		// 4 rolling moment 滚转力矩 绕x轴 改变滚转角速度p
		// 和侧滑角beta相关的系数
		Clbv = -1.402E-01 + 3.326E-02 * Mach - 7.590E-04 * Alpha + 8.596E-06 * (Alpha * Mach) - 3.794E-03 * pow(Mach, 2) + 2.354E-06 * pow(Alpha, 2)
			- 1.044E-08 * pow(Alpha * Mach, 2) + 2.219E-04 * pow(Mach, 3) - 8.964E-18 * pow(Alpha, 3) - 6.462E-06 * pow(Mach, 4) 
			+ 3.803E-19 * pow(Alpha, 4) + 7.419E-08 * pow(Mach, 5) - 3.353E-21 * pow(Alpha, 5);

		// 对右升降副翼的偏导
		Cl_dRE = 1.170E-04 + 2.794E-08 * Alpha * Mach - 1.160E-06 * 2 * RE	- 4.641E-11 * pow(Alpha * Mach, 2) * 2 * RE;

		// 对左升降副翼的偏导
		Cl_dLE = -(1.170E-04 + 2.794E-08 * Alpha * Mach - 1.160E-06 * 2 * LE - 4.641E-11 * pow(Alpha * Mach, 2) * 2 * LE);

		// 对方向舵的偏导
		Cl_dRUD = 0.00011441 - 2.6824E-06 * Alpha - 3.5496E-06 * Mach + 5.5547E-08 * Alpha * Mach;

		// 和偏航角速度r相关的系数
		Clr = +3.82E-01 - 1.06E-01 * Mach + 1.94E-03 * Alpha - 8.15E-05 * (Alpha * Mach) + 1.45E-02 * pow(Mach, 2)	- 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2) - 1.02E-03 * pow(Mach, 3)	- 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4) 
			+ 3.19E-08 * pow(Alpha, 4) - 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		// 和滚转角速度p相关的系数
		Clp = -2.99E-01 + 7.47E-02 * Mach + 1.38E-03 * Alpha - 8.78E-05 * (Alpha * Mach) - 9.13E-03 * pow(Mach, 2) - 2.04E-04 * pow(Alpha, 2)
			- 1.52E-07 * pow(Alpha * Mach, 2) + 5.73E-04 * pow(Mach, 3)	- 3.86E-05 * pow(Alpha, 3) - 1.79E-05 * pow(Mach, 4) + 4.21E-06 * pow(Alpha, 4)
			+ 2.20E-07 * pow(Mach, 5) - 1.15E-07 * pow(Alpha, 5);

		// 与舵角无关的气动力矩系数
		Cl_state = Clbv * Beta + Clp * (p * B / 2 / V) + Clr * (r * B / 2 / V);

		// 5 yawing moment 偏航力矩 绕z轴 改变偏航加速度r
		// 和侧滑角beta相关的系数
		Cnbv = Alpha * (6.9980E-4) + Mach * (5.9115E-2)	+ (Alpha * Mach) * (-7.5250E-5) + pow(Alpha, 2) * (2.5160E-4) + pow(Mach, 2) * (-1.4824E-2)
			+ pow(Alpha * Mach, 2) * (-2.1924E-7) + pow(Alpha, 3) * (-1.0777E-4) + pow(Mach, 3) * (1.2692E-3) + pow(Alpha * Mach, 3) * (1.0707E-8)
			+ pow(Alpha, 4) * (9.4989E-6) + pow(Mach, 4) * (-4.7098E-5)	+ pow(Alpha * Mach, 4) * (-5.5472E-11) + pow(Alpha, 5) * (-2.5953E-7)
			+ pow(Mach, 5) * (6.4284E-7) + pow(Alpha * Mach, 5) * (8.5863E-14);

		// 对左升降副翼的偏导
		Cn_dLE = -(- 1.30E-05 - 8.93E-08 * Alpha * Mach + 1.97E-06 * 2 * LE + 1.41E-11 * pow(Alpha * Mach, 2) * 2 * LE);

		// 对右升降副翼的偏导
		Cn_dRE = - 1.30E-05 - 8.93E-08 * Alpha * Mach + 1.97E-06 * 2 * RE + 1.41E-11 * pow(Alpha * Mach, 2) * 2 * RE;

		// 对方向舵的偏导
		Cn_dRUD = - 5.28E-04 + 1.39E-05 * Alpha + 1.65E-05 * Mach - 3.13E-07 * Alpha * Mach;

		// 和滚转角速度p相关的系数
		Cnp = 3.68E-01 - 9.79E-02 * Mach + 7.61E-16 * Alpha + 1.24E-02 * pow(Mach, 2) - 4.64E-16 * pow(Alpha, 2) - 8.05E-04 * pow(Mach, 3) 
			+ 1.01E-16 * pow(Alpha, 3) + 2.57E-05 * pow(Mach, 4) - 9.18E-18 * pow(Alpha, 4) - 3.20E-07 * pow(Mach, 5) + 2.96E-19 * pow(Alpha, 5);

		// 和偏航角速度r相关的系数
		Cnr = -2.41 + 5.96E-01 * Mach - 2.74E-03 * Alpha + 2.09E-04 * (Alpha * Mach) - 7.57E-02 * pow(Mach, 2) + 1.15E-03 * pow(Alpha, 2) 
			- 6.53E-08 * pow(Alpha * Mach, 2) + 4.90E-03 * pow(Mach, 3) - 3.87E-04 * pow(Alpha, 3) - 1.57E-04 * pow(Mach, 4)
			+ 3.60E-05 * pow(Alpha, 4) + 1.96E-06 * pow(Mach, 5) - 1.18E-06 * pow(Alpha, 5);

		// 与舵角无关的气动力矩系数
		Cn_state = Cnbv * Beta + Cnp * (p * B / 2 / V) + Cnr * (p * B / 2 / V);

		// 6 pitching moment 俯仰力矩 绕y轴 改变俯仰角速度q
		Cmbv = -2.192E-02 + 7.739E-03 * Mach - 2.260E-03 * Alpha + 1.808E-04 * (Alpha * Mach) - 8.849E-04 * pow(Mach, 2) + 2.616E-04 * pow(Alpha, 2)
			- 2.880E-07 * pow(Alpha * Mach, 2) + 4.617E-05 * pow(Mach, 3) - 7.887E-05 * pow(Alpha, 3) - 1.143E-06 * pow(Mach, 4)
			+ 8.288E-06 * pow(Alpha, 4) + 1.082E-08 * pow(Mach, 5) - 2.789E-07 * pow(Alpha, 5);

		// 对左升降副翼的偏导
		Cm_dLE = 2.89E-04 + 4.48E-06 * Alpha - 5.87E-06 * Mach + 9.72E-08 * Alpha * Mach;

		// 对右升降副翼的偏导
		Cm_dRE = 2.89E-04 + 4.48E-06 * Alpha - 5.87E-06 * Mach + 9.72E-08 * Alpha * Mach;

		// 对方向舵的偏导
		Cm_dRUD = 1.43E-07 * 4 * pow(RUD, 3) - 4.77E-22 * 5 * pow(RUD, 4) - 3.38E-10 * 6 * pow(RUD, 5) + 2.63E-24 * 7 * pow(RUD, 6);
		
		// 和俯仰角速度q相关的系数
		Cm_q = -1.36 + 3.86E-01 * Mach + 7.85E-04 * Alpha + 1.40E-04 * (Alpha * Mach) - 5.42E-02 * pow(Mach, 2)	+ 2.36E-03 * pow(Alpha, 2) 
			- 1.95E-06 * pow(Alpha * Mach, 2) + 3.80E-03 * pow(Mach, 3) - 1.48E-03 * pow(Alpha, 3) - 1.30E-04 * pow(Mach, 4)
			+ 1.69E-04 * pow(Alpha, 4) + 1.71E-06 * pow(Mach, 5) - 5.93E-06 * pow(Alpha, 5);

		// 与舵角无关的气动力矩系数
		Cm_state = Cmbv + Cm_q * (q * C / 2 / V);
	}

	ans[0][0] = CDbv; ans[0][1] = CD_dLE; ans[0][2] = CD_dRE; ans[0][3] = CD_dRUD;	// 沿着x轴的牵引力
	ans[1][0] = CYbv;  ans[1][1] = CY_dLE; ans[1][2] = CY_dRE; ans[1][3] = CY_dRUD;	// 沿着y轴的侧向力
	ans[2][0] = CLbv; ans[2][1] = CL_dLE; ans[2][2] = CL_dRE;						// 沿着z轴的升力
	ans[3][0] = Clbv; ans[3][1] = Cl_dLE; ans[3][2] = Cl_dRE; ans[3][3] = Cl_dRUD; ans[3][4] = Clp; ans[3][5] = Clr; ans[3][6] = Cl_state; 	// 引起相对于机体坐标系x轴的转动，产生滚转角速度p
	ans[4][0] = Cmbv; ans[4][1] = Cm_dLE; ans[4][2] = Cm_dRE; ans[4][3] = Cm_dRUD; ans[4][4] = Cm_q; ans[4][5] = Cm_state;					// 引起相对于机体坐标系y轴的转动，产生俯仰角速度q
	ans[5][0] = Cnbv; ans[5][1] = Cn_dLE; ans[5][2] = Cn_dRE; ans[5][3] = Cn_dRUD; ans[5][4] = Cnp; ans[5][5] = Cnr; ans[5][6] = Cn_state;	// 引起相对于机体坐标系z轴的转动，产生偏航角速度r

	return;
}