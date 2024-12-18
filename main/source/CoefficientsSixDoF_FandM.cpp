#include "CoefficientsSixDoF_FandM.h"

void CoefficientsSixDoF_FandM(
	const double& Mach, 
	const double& Alpha, 
	const double& Beta, 
	const double& Delta_e, 
	const double& Delta_a, 
	const double& Delta_r, 
	const double& C, 
	const double& B, 
	const double& Vel, 
	const double& omega_x, 
	const double& omega_y, 
	const double& omega_z, 
	vector<double>& ans) {

	double RE = Delta_a - Delta_e;
	double LE = -Delta_a - Delta_e;
	double RUD = Delta_r;

	double CL0, CL_e, CL_a, CD0, CD_e, CD_a, CD_r, CNb, CN_e, CN_a, CN_r;
	double mxb, mx_e, mx_a, mx_r, mxy, mxx, myb, my_e, my_a, my_r, myx, myy, mz0, mz_e, mz_a, mz_r, mzz;

	if (Mach < 0 || Mach > 30) {
		exit(1);
	}
	else if (Mach <= 1.25) {
		CL0 = -5.2491e-004 + Alpha * 1.5746e-002 + (Alpha * Mach) * 6.0213e-03
			- 3.4437e-004 * pow(Alpha, 2) + pow(Alpha * Mach, 2) * 1.4471E-04
			- 5.1952E-05 * pow(Alpha, 3) + 3.4771E-05 * pow(Alpha, 4)
			+ 2.7717E-03 * pow(Mach, 4) - 2.3034E-06 * pow(Alpha, 5);

		CL_e = -5.119E-04 + 1.000E-03 * Alpha - 1.406E-04 * (Alpha * RE)
			+ 1.313E-03 * (Alpha * Mach) - 8.584E-04 * (Mach * RE)
			+ 8.879E-05 * (Alpha * Mach) * RE - 1.604E-04 * pow(Mach, 2)
			- 3.477E-04 * pow(Alpha, 2) - 9.788E-05 * pow(Alpha * Mach, 2)
			- 1.703E-06 * pow(RE * Mach, 2) + 2.532E-05 * pow(Alpha, 3)
			- 3.727E-05 * pow(RE, 3) + 1.781E-07 * pow(RE, 2)
			+ 7.912E-07 * pow(Alpha * Mach * RE, 2) + 2.465E-08 * pow(Alpha * RE, 2)
			- 9.788E-05 * pow(Alpha * Mach, 2) - 5.942E-09 * pow(Alpha * Mach * RE, 3)
			- 7.377E-08 * pow(Alpha, 4) + 2.672E-08 * pow(RE, 4)
			- 1.610E-11 * pow(Alpha * Mach * RE, 4) - 3.273E-08 * pow(Alpha, 5)
			+ 7.624E-08 * pow(RE, 5) + 1.388E-13 * pow(Alpha * Mach * RE, 5);

		CL_a = -5.119E-04 + 1.000E-03 * Alpha - 1.406E-04 * (Alpha * LE)
			+ 1.313E-03 * (Alpha * Mach) - 8.584E-04 * (Mach * LE)
			+ 8.879E-05 * (Alpha * Mach) * LE - 1.604E-04 * pow(Mach, 2)
			- 3.477E-04 * pow(Alpha, 2) - 9.788E-05 * pow(Alpha * Mach, 2)
			- 1.703E-06 * pow(LE * Mach, 2) + 2.532E-05 * pow(Alpha, 3)
			- 3.727E-05 * pow(LE, 3) + 1.781E-07 * pow(LE, 2)
			+ 7.912E-07 * pow(Alpha * Mach * LE, 2) + 2.465E-08 * pow(Alpha * LE, 2)
			- 9.788E-05 * pow(Alpha * Mach, 2) - 5.942E-09 * pow(Alpha * Mach * LE, 3)
			- 7.377E-08 * pow(Alpha, 4) + 2.672E-08 * pow(LE, 4)
			- 1.610E-11 * pow(Alpha * Mach * LE, 4) - 3.273E-08 * pow(Alpha, 5)
			+ 7.624E-08 * pow(LE, 5) + 1.388E-13 * pow(Alpha * Mach * LE, 5);

		CD0 = +1.1457e-002 + CL0 * (-2.4645e-002) + Mach * (0)
			+ (CL0 * Mach) * (4.9698e-002) + pow(CL0, 2) * (-1.9112e+000)
			+ pow(Mach, 2) * (0) + pow(CL0 * Mach, 2) * (3.5404e+000)
			+ pow(CL0, 3) * (4.4334e+001) + pow(Mach, 3) * (0)
			+ pow(CL0 * Mach, 3) * (-7.0367e+001)
			+ pow(CL0, 4) * (-2.3841e+002) + pow(Mach, 4) * (0)
			+ pow(CL0 * Mach, 4) * (4.1750e+002) + pow(CL0, 5) * (4.1734e+002)
			+ pow(Mach, 5) * (5.4910e-002)
			+ pow(CL0 * Mach, 5) * (-7.9055e+002);

		CD_e = -5.184e-04 + 1.100e-03 * Alpha + 3.38e-07 * (Alpha * LE)
			- 1.36e-03 * (Alpha * Mach) - 2.79e-04 * (Mach * LE)
			- 1.53e-04 * (Alpha * Mach) * LE + 1.29e-03 * pow(Mach, 2)
			- 1.02e-04 * pow(Alpha, 2) + 9.39E-08 * pow(LE, 2)
			- 5.69E-07 * pow(Alpha * Mach * LE, 2) + 4.14E-07 * pow(Alpha * LE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * LE, 2)
			- 1.84E-06 * pow(LE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(LE, 4)
			+ 5.71E-09 * pow(LE, 5) - 8.93E-15 * pow(Alpha * Mach * LE, 5)
			- 7.58E-12 * pow(Alpha * Mach * LE, 4) - 3.94E-10 * pow(Alpha * Mach * LE, 3);

		CD_a = -5.184e-04 + 1.100e-03 * Alpha + 3.38e-07 * (Alpha * RE)
			- 1.36e-03 * (Alpha * Mach) - 2.79e-04 * (Mach * RE)
			- 1.53e-04 * (Alpha * Mach) * RE + 1.29e-03 * pow(Mach, 2)
			- 1.02e-04 * pow(Alpha, 2) + 9.39E-08 * pow(RE, 2)
			- 5.69E-07 * pow(Alpha * Mach * RE, 2) + 4.14E-07 * pow(Alpha * RE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * RE, 2)
			- 1.84E-06 * pow(RE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(RE, 4)
			+ 5.71E-09 * pow(RE, 5) - 8.93E-15 * pow(Alpha * Mach * RE, 5)
			- 7.58E-12 * pow(Alpha * Mach * RE, 4) - 3.94E-10 * pow(Alpha * Mach * RE, 3);

		CD_r = +2.47E-04 - 1.93E-04 * Alpha + 7.27E-05 * (Alpha * Mach)
			+ 4.73E-05 * pow(Mach, 2) + 1.50E-05 * pow(Alpha, 2) + 5.03E-06 * pow(RUD, 2)
			- 1.30E-07 * pow(Alpha * Mach * RUD, 2) - 3.50E-08 * pow(Alpha * RUD, 2)
			- 1.68E-06 * pow(Alpha * Mach, 2) + 4.53E-06 * pow(Mach * RUD, 2)
			- 1.98E-11 * pow(Alpha, 3) - 2.63E-08 * pow(Alpha, 4) + 7.54E-09 * pow(RUD, 4)
			+ 3.12E-12 * pow(Alpha * Mach * RUD, 4);

		CNb = -4.750E-01 - 5.000E-02 * Mach;

		CN_e = -(-1.845E-04 * Mach - 2.13E-07 * (Alpha * RE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * RE)
			+ 6.17E-08 * (Alpha * Mach) * RE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(RE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(RE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * RE, 2) - 2.66E-09 * pow(Alpha * RE, 2)
			+ 3.92E-08 * pow(Mach * RE, 2) + 5.42E-11 * pow(Alpha * Mach * RE, 3)
			- 4.73E-10 * pow(RE, 4) + 7.35E-14 * pow(Alpha * Mach * RE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(RE, 5)
			- 1.11E-15 * pow(Alpha * Mach * RE, 5));

		CN_a = -1.845E-04 * Mach - 2.13E-07 * (Alpha * LE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * LE)
			+ 6.17E-08 * (Alpha * Mach) * LE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(LE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(LE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * LE, 2) - 2.66E-09 * pow(Alpha * LE, 2)
			+ 3.92E-08 * pow(Mach * LE, 2) + 5.42E-11 * pow(Alpha * Mach * LE, 3)
			- 4.73E-10 * pow(LE, 4) + 7.35E-14 * pow(Alpha * Mach * LE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(LE, 5)
			- 1.11E-15 * pow(Alpha * Mach * LE, 5);

		CN_r = +2.440E-03 * RUD;

		mxb = -9.380E-02 - 1.250E-02 * Mach;

		mx_e = -(5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * RE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * RE)
			- 8.500E-06 * (Alpha * Mach) * RE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(RE, 2)
			- 5.302E-08 * pow(Alpha * Mach * RE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * RE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(RE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * RE, 3) - 3.576E-08 * pow(Alpha * RE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(RE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * RE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(RE, 5) - 2.984E-14 * pow(Alpha * Mach * RE, 5));

		mx_a = 5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * LE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * LE)
			- 8.500E-06 * (Alpha * Mach) * LE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(LE, 2)
			- 5.302E-08 * pow(Alpha * Mach * LE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * LE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(LE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * LE, 3) - 3.576E-08 * pow(Alpha * LE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(LE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * LE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(LE, 5) - 2.984E-14 * pow(Alpha * Mach * LE, 5);

		mx_r = +7.000000E-04 * RUD;
		mxy = +2.625000E-01 + 2.50E-02 * (Mach);
		mxx = -1.337500E-01 - 1.250000E-02 * (Mach);

		myb = +1.062E-01 + 6.250E-02 * Mach;

		my_e = -(-2.7e-7 * (Alpha * RE) - 1.008E-05 * (Mach * RE)
			+ 3.564E-07 * (Alpha * Mach) * RE + 1.1e-7 * pow(RE, 3) + 1.11e-07 * pow(RE, 3)
			- 9.32E-12 * pow(Alpha * Mach * RE, 3) - 1.9910e-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(RE, 4) + 1.82E-28 * pow(Alpha * Mach * RE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046e-010 * pow(RE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * RE, 5));

		my_a = -2.7e-7 * (Alpha * LE) - 1.008E-05 * (Mach * LE)
			+ 3.564E-07 * (Alpha * Mach) * LE + 1.1e-7 * pow(LE, 3) + 1.11e-07 * pow(LE, 3)
			- 9.32E-12 * pow(Alpha * Mach * LE, 3) - 1.9910e-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(LE, 4) + 1.82E-28 * pow(Alpha * Mach * LE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046e-010 * pow(LE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * LE, 5);

		my_r = -3.000E-03 * RUD;
		myx = +1.790E-01 + 2.000E-02 * Mach;
		myy = -1.2787 - 1.375e-001 * Mach;

		mz0 = +(-1.8316e-003) + CL0 * (-1.0306e-001) + Mach * (0)
			+ (CL0 * Mach) * (-1.8335e-001) + pow(CL0, 2) * (-1.1839e+000)
			+ pow(Mach, 2) * (-2.8113e-03)
			+ pow(CL0 * Mach, 2) * (-1.3362e+00) + pow(CL0, 3) * (9.0641e+00)
			+ pow(Mach, 3) * (0) + pow(CL0 * Mach, 3) * (2.6964e+001)
			+ pow(CL0, 4) * (-6.3590e+01) + pow(Mach, 4) * (0)
			+ pow(CL0 * Mach, 4) * (-8.0921e+01)
			+ pow(CL0, 5) * (1.6885e+02) + pow(Mach, 5) * (0)
			+ pow(CL0 * Mach, 5) * (-4.2209e+00);

		mz_e = +2.880000E-04 - 5.351000E-04 * Alpha + 4.550000E-05 * (Alpha * RE)
			+ 3.379000E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * RE)
			- 2.770E-05 * (Alpha * Mach) * RE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(RE, 2)
			- 1.000E-07 * pow(Alpha * Mach * RE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * RE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(RE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * RE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(RE, 4) + 1.29E-12 * pow(Alpha * Mach * RE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(RE, 5)
			- 4.67E-14 * pow(Alpha * Mach * RE, 5);

		mz_a = +2.880000E-04 - 5.351000E-04 * Alpha + 4.550000E-05 * (Alpha * LE)
			+ 3.379000E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * LE)
			- 2.770E-05 * (Alpha * Mach) * LE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(LE, 2)
			- 1.000E-07 * pow(Alpha * Mach * LE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * LE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(LE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * LE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(LE, 4) + 1.29E-12 * pow(Alpha * Mach * LE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(LE, 5)
			- 4.67E-14 * pow(Alpha * Mach * LE, 5);

		mz_r = -1.841E-04 + 3.5E-06 * Alpha + 2.762E-04 * Mach - 1.0E-07 * RUD
			- 4.0E-07 * pow(Alpha, 2) + 5.8E-06 * pow(RUD, 2)
			+ 6.482E-09 * pow(Alpha * Mach * RUD, 2);

		mzz = -1.0313 - 3.125000E-01 * Mach;
	}
	else if (Mach < 4) {
		CL0 = +1.9920e-001 + Mach * (2.3402e-001) + Alpha * (3.8202e-002)
			+ (Alpha * Mach) * (-2.4626e-003) + pow(Mach, 2) * (-6.4872e-001)
			+ pow(Alpha, 2) * (-6.9523e-003)
			+ pow(Alpha * Mach * Mach, 2) * (4.5735e-006)
			+ pow(Alpha * Alpha * Mach, 2) * (2.1241e-007)
			+ pow(Alpha * Mach, 2) * (-1.0521e-004)
			+ pow(Alpha * Mach, 4) * (-9.5825e-009)
			+ pow(Mach, 3) * (3.9121e-001)
			+ pow(Alpha, 3) * (1.0295e-003) + pow(Mach, 4) * (-9.1356e-002)
			+ pow(Alpha, 4) * (-5.7398e-005) + pow(Mach, 5) * (7.4089e-003)
			+ pow(Alpha, 5) * (1.0934e-006);

		CL_e = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.3093e-005) + (Alpha * Mach) * (0)
			+ (Mach * RE) * (-1.4287e-004)
			+ ((Alpha * Mach) * RE) * (6.1071e-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (2.7242e-004)
			+ pow(Alpha * Mach * RE, 2) * (-9.1890e-008)
			+ pow(Alpha * RE, 2) * (3.4060e-007)
			+ pow(Alpha * Mach, 2) * (-6.5093e-006)
			+ pow(Mach * RE, 2) * (-6.3863e-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092e-004)
			+ pow(RE, 3) * (3.8067e-006)
			+ pow(Alpha * Mach * RE, 3) * (2.3165e-011)
			+ pow(Mach, 4) * (-1.0680e-003)
			+ pow(Alpha, 4) * (-2.1893e-005) + pow(RE, 4) * (-3.7716e-007)
			+ pow(Alpha * Mach * RE, 4) * (7.9006e-014)
			+ pow(Mach, 5) * (2.6056e-004)
			+ pow(Alpha, 5) * (9.2099e-007) + pow(RE, 5) * (-8.5345e-009)
			+ pow(Alpha * Mach * RE, 5) * (-2.5698e-017);

		CL_a = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.3093e-005) + (Alpha * Mach) * (0)
			+ (Mach * LE) * (-1.4287e-004)
			+ ((Alpha * Mach) * LE) * (6.1071e-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (2.7242e-004)
			+ pow(Alpha * Mach * LE, 2) * (-9.1890e-008)
			+ pow(Alpha * LE, 2) * (3.4060e-007)
			+ pow(Alpha * Mach, 2) * (-6.5093e-006)
			+ pow(Mach * LE, 2) * (-6.3863e-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092e-004)
			+ pow(LE, 3) * (3.8067e-006)
			+ pow(Alpha * Mach * LE, 3) * (2.3165e-011)
			+ pow(Mach, 4) * (-1.0680e-003)
			+ pow(Alpha, 4) * (-2.1893e-005) + pow(LE, 4) * (-3.7716e-007)
			+ pow(Alpha * Mach * LE, 4) * (7.9006e-014)
			+ pow(Mach, 5) * (2.6056e-004)
			+ pow(Alpha, 5) * (9.2099e-007) + pow(LE, 5) * (-8.5345e-009)
			+ pow(Alpha * Mach * LE, 5) * (-2.5698e-017);

		CD0 = +(-8.2073e-002) + CL0 * (-9.1273e-002)
			+ Mach * (2.1845e-001)
			+ (CL0 * Mach) * (3.2202e-002) + pow(CL0, 2) * (1.6325e+000)
			+ pow(Mach, 2) * (-1.3680e-001)
			+ pow(CL0 * Mach, 2) * (5.7526e-002)
			+ pow(CL0, 3) * (-1.1575e+000) + pow(Mach, 3) * (3.8791e-002)
			+ pow(CL0 * Mach, 3) * (-2.4002e-001)
			+ pow(CL0, 4) * (-8.5306e+000)
			+ pow(Mach, 4) * (-5.2527e-003)
			+ pow(CL0 * Mach, 4) * (3.5543e-001)
			+ pow(CL0, 5) * (1.7259e+001) + pow(Mach, 5) * (2.7435e-004)
			+ pow(CL0 * Mach, 5) * (-1.4983e-001);

		CD_e = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.6923e-005) + (Alpha * Mach) * (1.5100e-005)
			+ (Mach * RE) * (1.3641e-007)
			+ ((Alpha * Mach) * RE) * (5.1142e-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (1.2125e-005)
			+ pow(Alpha * Mach * RE, 2) * (3.5662e-009)
			+ pow(Alpha * RE, 2) * (-1.3848e-008)
			+ pow(Alpha * Mach, 2) * (-4.7972e-007)
			+ pow(Mach * RE, 2) * (-3.3763e-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045e-008)
			+ pow(RE, 3) * (3.9119e-008)
			+ pow(Alpha * Mach * RE, 3) * (-9.7714e-013)
			+ pow(Mach, 4) * (9.6475e-007)
			+ pow(Alpha, 4) * (1.5015e-008) + pow(RE, 4) * (4.5137e-009)
			+ pow(Alpha * Mach * RE, 4) * (-6.6207e-016)
			+ pow(Mach, 5) * (-3.2682e-007)
			+ pow(Alpha, 5) * (-3.5360e-010) + pow(RE, 5) * (-1.1538e-010)
			+ pow(Alpha * Mach * RE, 5) * (4.1917e-019);

		CD_a = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.6923e-005) + (Alpha * Mach) * (1.5100e-005)
			+ (Mach * LE) * (1.3641e-007)
			+ ((Alpha * Mach) * LE) * (5.1142e-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (1.2125e-005)
			+ pow(Alpha * Mach * LE, 2) * (3.5662e-009)
			+ pow(Alpha * LE, 2) * (-1.3848e-008)
			+ pow(Alpha * Mach, 2) * (-4.7972e-007)
			+ pow(Mach * LE, 2) * (-3.3763e-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045e-008)
			+ pow(LE, 3) * (3.9119e-008)
			+ pow(Alpha * Mach * LE, 3) * (-9.7714e-013)
			+ pow(Mach, 4) * (9.6475e-007)
			+ pow(Alpha, 4) * (1.5015e-008) + pow(LE, 4) * (4.5137e-009)
			+ pow(Alpha * Mach * LE, 4) * (-6.6207e-016)
			+ pow(Mach, 5) * (-3.2682e-007)
			+ pow(Alpha, 5) * (-3.5360e-010) + pow(LE, 5) * (-1.1538e-010)
			+ pow(Alpha * Mach * LE, 5) * (4.1917e-019);

		CD_r = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.6425e-021)
			+ (Alpha * Mach) * (-9.8380e-006)
			+ (Mach * RUD) * (1.8193e-020)
			+ ((Alpha * Mach) * RUD) * (1.0319e-021)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RUD, 2) * (8.7608e-006)
			+ pow(Alpha * Mach * RUD, 2) * (5.4045e-010)
			+ pow(Alpha * RUD, 2) * (-2.8939e-008)
			+ pow(Alpha * Mach, 2) * (2.1842e-007)
			+ pow(Mach * RUD, 2) * (-2.9646e-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.0067e-007)
			+ pow(RUD, 3) * (-8.8556e-022)
			+ pow(Alpha * Mach * RUD, 3) * (-5.2022e-027)
			+ pow(Mach, 4) * (1.3388e-006) + pow(Alpha, 4) * (1.6460e-007)
			+ pow(RUD, 4) * (4.6754e-010)
			+ pow(Alpha * Mach * RUD, 4) * (2.6560e-016)
			+ pow(Mach, 5) * (-2.5185e-007)
			+ pow(Alpha, 5) * (-7.2766e-009) + pow(RUD, 5) * (1.5611e-024)
			+ pow(Alpha * Mach * RUD, 5) * (5.4442e-033);

		CNb = +(0) + Mach * (0) + Alpha * (-1.1185e-002)
			+ (Alpha * Mach) * (3.0432e-003) + pow(Mach, 2) * (-3.7586e-001)
			+ pow(Alpha, 2) * (3.4004e-003)
			+ pow(Alpha * Mach * Mach, 2) * (-2.4047e-006)
			+ pow(Alpha * Alpha * Mach, 2) * (3.6104e-007)
			+ pow(Alpha * Mach, 2) * (-8.7176e-005)
			+ pow(Alpha * Mach, 4) * (-5.3622e-010) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8160e-004) + pow(Mach, 4) * (9.4289e-002)
			+ pow(Alpha, 4) * (4.4848e-005) + pow(Mach, 5) * (-1.8384e-002)
			+ pow(Alpha, 5) * (-1.3021e-006);

		CN_e = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(RE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * RE, 2));

		CN_a = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(LE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * LE, 2));

		CN_r = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.0067e-005)
			+ (Alpha * Mach) * (0) + (Mach * RUD) * (-5.7185e-004)
			+ ((Alpha * Mach) * RUD) * (-1.5307e-005) + pow(Mach, 2) * (0)
			+ pow(Alpha, 2) * (0) + pow(RUD, 2) * (1.9243e-019)
			+ pow(Alpha * Mach * RUD, 2) * (2.8011e-022)
			+ pow(Alpha * RUD, 2) * (-2.0404e-021)
			+ pow(Alpha * Mach, 2) * (-1.2673e-020)
			+ pow(Mach * RUD, 2) * (-1.7950e-020)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.9873e-019)
			+ pow(RUD, 3) * (3.2768e-005)
			+ pow(Alpha * Mach * RUD, 3) * (1.2674e-012)
			+ pow(Mach, 4) * (-3.8438e-020)
			+ pow(Alpha, 4) * (1.9239e-019) + pow(RUD, 4) * (7.7275e-023)
			+ pow(Alpha * Mach * RUD, 4) * (-3.2592e-029)
			+ pow(Mach, 5) * (3.1048e-020)
			+ pow(Alpha, 5) * (-9.0794e-021) + pow(RUD, 5) * (-6.5825e-008)
			+ pow(Alpha * Mach * RUD, 5) * (1.2684e-017);

		mxb = +(0) + Mach * (0) + Alpha * (5.9211e-004)
			+ (Alpha * Mach) * (-3.1579e-004) + pow(Mach, 2) * (-8.7296e-002)
			+ pow(Alpha, 2) * (-5.7398e-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1037e-006)
			+ pow(Alpha * Alpha * Mach, 2) * (-6.8068e-008)
			+ pow(Alpha * Mach, 2) * (2.0549e-005)
			+ pow(Alpha * Mach, 4) * (3.6561e-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-2.8226e-016) + pow(Mach, 4) * (2.0334e-002)
			+ pow(Alpha, 4) * (1.9013e-007) + pow(Mach, 5) * (-3.7733e-003)
			+ pow(Alpha, 5) * (-9.6648e-019);

		mx_e = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2));

		mx_a = 3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2);

		mx_r = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		mxy = +3.82E-01 - 1.06E-01 * Mach
			+ 1.94E-03 * Alpha - 8.15E-05 * (Alpha * Mach)
			+ 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2)
			- 1.02E-03 * pow(Mach, 3) - 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4)
			+ 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		mxx = +(0) + Mach * (0) + Alpha * (-1.2668e-005)
			+ (Alpha * Mach) * (1.7282e-005) + pow(Mach, 2) * (-1.0966e-001)
			+ pow(Alpha, 2) * (1.0751e-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.0989e-006)
			+ pow(Alpha * Alpha * Mach, 2) * (6.1850e-009)
			+ pow(Alpha * Mach, 2) * (8.6481e-006)
			+ pow(Alpha * Mach, 4) * (-4.3707e-010)
			+ pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-1.1567e-005) + pow(Mach, 4) * (2.6725e-002)
			+ pow(Alpha, 4) * (1.5082e-006) + pow(Mach, 5) * (-5.0800e-003)
			+ pow(Alpha, 5) * (-6.1276e-008);

		myb = +(0) + Mach * (0) + Alpha * (-2.3745e-003)
			+ (Alpha * Mach) * (8.5307e-004)
			+ pow(Mach, 2) * (1.4474e-001)
			+ pow(Alpha, 2) * (5.3105e-004)
			+ pow(Alpha * Mach * Mach, 2) * (-8.3462e-007)
			+ pow(Alpha * Alpha * Mach, 2) * (1.3335e-007)
			+ pow(Alpha * Mach, 2) * (-2.7081e-005)
			+ pow(Alpha * Mach, 4) * (-1.3450e-009)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.1046e-005)
			+ pow(Mach, 4) * (-3.9519e-002) + pow(Alpha, 4) * (-1.5141e-006)
			+ pow(Mach, 5) * (7.7646e-003) + pow(Alpha, 5) * (1.7278e-007);

		my_a = +2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(LE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * LE, 2);

		my_e = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(RE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * RE, 2));

		my_r = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD) - 3.13E-07 * (Alpha * Mach) * RUD;

		myx = +(1.7000e-001) + Alpha * (-6.4056e-018)
			+ Mach * (1.1333e-002) + (Alpha * Mach) * (2.3467e-018)
			+ pow(Alpha, 2) * (2.0917e-019)
			+ pow(Mach, 2) * (-5.3333e-003)
			+ pow(Alpha * Mach, 2) * (-5.0665e-020);

		myy = +(0) + Mach * (0) + Alpha * (-1.3332e-003)
			+ (Alpha * Mach) * (6.6899e-004)
			+ pow(Mach, 2) * (-1.0842e+000)
			+ pow(Alpha, 2) * (1.6434e-003)
			+ pow(Alpha * Mach * Mach, 2) * (-4.4258e-006)
			+ pow(Alpha * Alpha * Mach, 2) * (1.2017e-007)
			+ pow(Alpha * Mach, 2) * (1.0819e-005)
			+ pow(Alpha * Mach, 4) * (-2.8899e-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8118e-004) + pow(Mach, 4) * (2.7379e-001)
			+ pow(Alpha, 4) * (6.7994e-005) + pow(Mach, 5) * (-5.2435e-002)
			+ pow(Alpha, 5) * (-2.5848e-006);

		mz0 = +(-5.7643e-001) + Mach * (1.0553e+000) + CL0 * (-3.7951e-001)
			+ (CL0 * Mach) * (1.0483e-001) + pow(Mach, 2) * (-7.4344e-001)
			+ pow(CL0, 2) * (-1.5412e-001)
			+ pow(CL0 * Mach * Mach, 2) * (-2.1133e-003)
			+ pow(CL0 * CL0 * Mach, 2) * (-1.7858e-001)
			+ pow(CL0 * Mach, 2) * (5.7805e-002)
			+ pow(CL0 * Mach, 4) * (-3.8875e-003)
			+ pow(Mach, 3) * (2.5341e-001)
			+ pow(CL0, 3) * (-4.9731e-001) + pow(Mach, 4) * (-4.1938e-002)
			+ pow(CL0, 4) * (7.1784e+000) + pow(Mach, 5) * (2.7017e-003)
			+ pow(CL0, 5) * (-1.0331e+001);

		mz_e = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE) + 9.72E-08 * (Alpha * Mach) * RE;

		mz_a = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE) + 9.72E-08 * (Alpha * Mach) * LE;

		mz_r = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

		mzz = +(0) + Mach * (0) + Alpha * (-1.0828e-002)
			+ (Alpha * Mach) * (4.2311e-003)
			+ pow(Mach, 2) * (-6.1171e-001)
			+ pow(Alpha, 2) * (4.6974e-003)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1593e-005)
			+ pow(Alpha * Alpha * Mach, 2) * (2.5378e-007)
			+ pow(Alpha * Mach, 2) * (-7.0964e-005)
			+ pow(Alpha * Mach, 4) * (4.1284e-008)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-1.1414e-003)
			+ pow(Mach, 4) * (1.5903e-001)
			+ pow(Alpha, 4) * (1.1176e-004) + pow(Mach, 5) * (-3.0665e-002)
			+ pow(Alpha, 5) * (-3.8123e-006);
	}
	else {
		CL0 = -8.19E-02 + 4.70E-02 * Mach + 1.86E-02 * Alpha
			- 4.73E-04 * (Alpha * Mach) - 9.19E-03 * pow(Mach, 2) - 1.52E-04 * pow(Alpha, 2)
			+ 5.99E-07 * pow(Alpha * Mach, 2) + 7.74E-04 * pow(Mach, 3)
			+ 4.08E-06 * pow(Alpha, 3) - 2.93E-05 * pow(Mach, 4) - 3.91E-07 * pow(Alpha, 4)
			+ 4.12E-07 * pow(Mach, 5) + 1.30E-08 * pow(Alpha, 5);

		CL_e = -1.45E-05 + 1.01E-04 * Alpha + 7.10E-06 * Mach - 4.14E-04 * RE
			- 3.51E-06 * (Alpha * RE) + 4.70E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * RE) - 1.70E-07 * (Alpha * Mach) * RE;

		CL_a = -1.45E-05 + 1.01E-04 * Alpha + 7.10E-06 * Mach - 4.14E-04 * LE
			- 3.51E-06 * (Alpha * LE) + 4.70E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * LE) - 1.70E-07 * (Alpha * Mach) * LE;

		CD0 = +8.717E-02 - 3.307E-02 * Mach + 3.179E-03 * Alpha
			- 1.250E-04 * (Alpha * Mach) + 5.036E-03 * pow(Mach, 2)
			- 1.100E-03 * pow(Alpha, 2) + 1.405E-07 * pow(Alpha * Mach, 2)
			- 3.658E-04 * pow(Mach, 3) + 3.175E-04 * pow(Alpha, 3) + 1.274E-05 * pow(Mach, 4)
			- 2.985E-05 * pow(Alpha, 4) - 1.705E-07 * pow(Mach, 5) + 9.766E-07 * pow(Alpha, 5);

		CD_e = +1 * (4.5548e-004) + Alpha * (2.5411e-005) + Mach * (-1.1436e-004)
			+ RE * (-3.6417e-005) + ((Alpha * Mach) * RE) * (-5.3015e-007)
			+ pow(Alpha, 2) * (3.2187e-006) + pow(Mach, 2) * (3.0140e-006)
			+ pow(RE, 2) * (6.9629e-006)
			+ pow(Alpha * Mach * RE, 2) * (2.1026e-012);

		CD_a = +1 * (4.5548e-004) + Alpha * (2.5411e-005) + Mach * (-1.1436e-004)
			+ LE * (-3.6417e-005) + ((Alpha * Mach) * LE) * (-5.3015e-007)
			+ pow(Alpha, 2) * (3.2187e-006) + pow(Mach, 2) * (3.0140e-006)
			+ pow(LE, 2) * (6.9629e-006)
			+ pow(Alpha * Mach * LE, 2) * (2.1026e-012);

		CD_r = +7.50E-04 - 2.29E-05 * Alpha - 9.69E-05 * Mach - 1.83E-06 * RUD
			+ 9.13E-09 * (Alpha * Mach) * RUD + 8.76E-07 * pow(Alpha, 2)
			+ 2.70E-06 * pow(Mach, 2) + 1.97E-06 * pow(RUD, 2)
			- 1.77E-11 * pow(Alpha * Mach * RUD, 2);

		CNb = +(0) + Mach * (-2.9253e-001) + Alpha * (2.8803e-003)
			+ (Alpha * Mach) * (-2.8943e-004) + pow(Mach, 2) * (5.4822e-002)
			+ pow(Alpha, 2) * (7.3535e-004)
			+ pow(Alpha * Mach * Mach, 2) * (-4.6490e-009)
			+ pow(Alpha * Alpha * Mach, 2) * (-2.0675e-008)
			+ pow(Alpha * Mach, 2) * (4.6205e-006)
			+ pow(Alpha * Mach, 4) * (2.6144e-011)
			+ pow(Mach, 3) * (-4.3203e-003)
			+ pow(Alpha, 3) * (-3.7405e-004) + pow(Mach, 4) * (1.5495e-004)
			+ pow(Alpha, 4) * (2.8183e-005) + pow(Mach, 5) * (-2.0829e-006)
			+ pow(Alpha, 5) * (-5.2083e-007);

		CN_e = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(RE, 2) + 4.50E-14 * pow(Alpha * Mach * RE, 2));

		CN_a = -1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(LE, 2) + 4.50E-14 * pow(Alpha * Mach * LE, 2);

		CN_r = -1.43E-18 + 4.86E-20 * Alpha + 1.86E-19 * Mach + 3.84E-04 * RUD
			- 1.17E-05 * (Alpha * RUD) - 1.07E-05 * (Mach * RUD)
			+ 2.60E-07 * (Alpha * Mach) * RUD;

		mxb = -1.402E-01 + 3.326E-02 * Mach - 7.590E-04 * Alpha
			+ 8.596E-06 * (Alpha * Mach) - 3.794E-03 * pow(Mach, 2)
			+ 2.354E-06 * pow(Alpha, 2) - 1.044E-08 * pow(Alpha * Mach, 2)
			+ 2.219E-04 * pow(Mach, 3) - 8.964E-18 * pow(Alpha, 3) - 6.462E-06 * pow(Mach, 4)
			+ 3.803E-19 * pow(Alpha, 4) + 7.419E-08 * pow(Mach, 5) - 3.353E-21 * pow(Alpha, 5);

		mx_a = +3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2);

		mx_e = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2));

		mx_r = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		mxy = +3.82E-01 - 1.06E-01 * Mach + 1.94E-03 * Alpha
			- 8.15E-05 * (Alpha * Mach) + 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2) - 1.02E-03 * pow(Mach, 3)
			- 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4) + 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		mxx = -2.99E-01 + 7.47E-02 * Mach + 1.38E-03 * Alpha
			- 8.78E-05 * (Alpha * Mach) - 9.13E-03 * pow(Mach, 2) - 2.04E-04 * pow(Alpha, 2)
			- 1.52E-07 * pow(Alpha * Mach, 2) + 5.73E-04 * pow(Mach, 3)
			- 3.86E-05 * pow(Alpha, 3) - 1.79E-05 * pow(Mach, 4) + 4.21E-06 * pow(Alpha, 4)
			+ 2.20E-07 * pow(Mach, 5) - 1.15E-07 * pow(Alpha, 5);

		myb = +(0) + Alpha * (6.9980e-004) + Mach * (5.9115e-002)
			+ (Alpha * Mach) * (-7.5250e-005) + pow(Alpha, 2) * (2.5160e-004)
			+ pow(Mach, 2) * (-1.4824e-002)
			+ pow(Alpha * Mach, 2) * (-2.1924e-007)
			+ pow(Alpha, 3) * (-1.0777e-004) + pow(Mach, 3) * (1.2692e-003)
			+ pow(Alpha * Mach, 3) * (1.0707e-008)
			+ pow(Alpha, 4) * (9.4989e-006) + pow(Mach, 4) * (-4.7098e-005)
			+ pow(Alpha * Mach, 4) * (-5.5472e-011)
			+ pow(Alpha, 5) * (-2.5953e-007) + pow(Mach, 5) * (6.4284e-007)
			+ pow(Alpha * Mach, 5) * (8.5863e-014);

		my_e = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(RE, 2) + 1.41E-11 * pow(Alpha * Mach * RE, 2));

		my_a = 2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(LE, 2) + 1.41E-11 * pow(Alpha * Mach * LE, 2);

		my_r = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD)
			- 3.13E-07 * (Alpha * Mach) * RUD;

		myx = +3.68E-01 - 9.79E-02 * Mach + 7.61E-16 * Alpha + 1.24E-02 * pow(Mach, 2)
			- 4.64E-16 * pow(Alpha, 2) - 8.05E-04 * pow(Mach, 3) + 1.01E-16 * pow(Alpha, 3)
			+ 2.57E-05 * pow(Mach, 4)
			- 9.18E-18 * pow(Alpha, 4) - 3.20E-07 * pow(Mach, 5) + 2.96E-19 * pow(Alpha, 5);

		myy = -2.41E+00 + 5.96E-01 * Mach - 2.74E-03 * Alpha
			+ 2.09E-04 * (Alpha * Mach) - 7.57E-02 * pow(Mach, 2)
			+ 1.15E-03 * pow(Alpha, 2) - 6.53E-08 * pow(Alpha * Mach, 2)
			+ 4.90E-03 * pow(Mach, 3) - 3.87E-04 * pow(Alpha, 3) - 1.57E-04 * pow(Mach, 4)
			+ 3.60E-05 * pow(Alpha, 4) + 1.96E-06 * pow(Mach, 5) - 1.18E-06 * pow(Alpha, 5);

		mz0 = -2.192E-02 + 7.739E-03 * Mach - 2.260E-03 * Alpha
			+ 1.808E-04 * (Alpha * Mach) - 8.849E-04 * pow(Mach, 2)
			+ 2.616E-04 * pow(Alpha, 2) - 2.880E-07 * pow(Alpha * Mach, 2)
			+ 4.617E-05 * pow(Mach, 3) - 7.887E-05 * pow(Alpha, 3) - 1.143E-06 * pow(Mach, 4)
			+ 8.288E-06 * pow(Alpha, 4) + 1.082E-08 * pow(Mach, 5) - 2.789E-07 * pow(Alpha, 5);

		mz_a = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE)
			+ 9.72E-08 * (Alpha * Mach) * LE;

		mz_e = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE)
			+ 9.72E-08 * (Alpha * Mach) * RE;

		mz_r = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

		mzz = -1.36E+00 + 3.86E-01 * Mach + 7.85E-04 * Alpha
			+ 1.40E-04 * (Alpha * Mach) - 5.42E-02 * pow(Mach, 2)
			+ 2.36E-03 * pow(Alpha, 2) - 1.95E-06 * pow(Alpha * Mach, 2)
			+ 3.80E-03 * pow(Mach, 3) - 1.48E-03 * pow(Alpha, 3) - 1.30E-04 * pow(Mach, 4)
			+ 1.69E-04 * pow(Alpha, 4) + 1.71E-06 * pow(Mach, 5) - 5.93E-06 * pow(Alpha, 5);
	}

	double CD = CD0 + CD_e + CD_a + CD_r;
	double CL = CL0 + CL_e + CL_a;
	double CN = CNb * Beta + CN_e + CN_a + CN_r;
	double mx = mxb * Beta + mx_e + mx_a + mx_r + (mxx * omega_x + mxy * omega_y) * B / (2 * Vel);
	double my = myb * Beta + my_e + my_a + my_r + (myx * omega_x + myy * omega_y) * B / (2 * Vel);
	double mz = mz0 + mz_e + mz_a + mz_r + mzz * omega_z * C / (2 * Vel);

	ans[0] = CD;
	ans[1] = CL;
	ans[2] = CN;
	ans[3] = mx;
	ans[4] = my;
	ans[5] = mz;

	return;
}