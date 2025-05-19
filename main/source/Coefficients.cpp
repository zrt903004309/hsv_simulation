#include <Coefficients.h>

void getCoefficients(
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
	std::vector<double>& ans,
	int& flag) {

	// 和《Six DoF Nonlinear Equations of Motion for a Generic Hypersonic Vehicle》中的左右襟翼相反
	// delta_a对应 右升降副翼 delta_e对应 左升降副翼 delta_r 对应 方向舵

	// 这里为什么要进行角度的一个简单映射？
	// 因为控制律计算出的是三通道虚拟控制量(等效舵偏)，针对alpha，beta，mu三通道，
	// 但是计算气动系数需要的是飞行器的可操控舵面(实际舵偏)，所以这两个之间有个转换矩阵

	// double RE = Delta_a - Delta_e;		// Right Elevon 右襟翼
	// double LE = -Delta_a - Delta_e;		// Left Elevon 左襟翼
	// double RUD = Delta_r;
	//double RE = Delta_phi - Delta_alpha;		// Right Elevon 右襟翼
	//double LE = -Delta_phi - Delta_alpha;		// Left Elevon 左襟翼
	//double RUD = -Delta_beta;
	double LE = Delta_e;		// Left Elevon 左襟翼
	double RE = Delta_a;		// Right Elevon 右襟翼
	double RUD = Delta_r;

	double CLbv, CL_RE, CL_LE, CDbv, CD_RE, CD_LE, CD_RUD, CYB, CY_RE, CY_LE, CY_RUD;
	double Cllbv, Cll_RE, Cll_LE, Cll_RUD, Cllr, Cllp, Cnbv, Cn_RE, Cn_LE, Cn_RUD, Cnp, Cnr, Cmbv, Cm_RE, Cm_LE, Cm_RUD, Cm_q;
	if (Mach < 0 || Mach > 30) {
		std::cout << Delta_e << " " << Mach << "坠毁了" << std::endl;
		exit(1);
	}
	else if (Mach <= 1.25) {
		CLbv = -5.2491E-04 + Alpha * 1.5746E-02 + (Alpha * Mach) * 6.0213E-03
			- 3.4437E-04 * pow(Alpha, 2) + pow(Alpha * Mach, 2) * 1.4471E-04
			- 5.1952E-05 * pow(Alpha, 3) + 3.4771E-05 * pow(Alpha, 4)
			+ 2.7717E-03 * pow(Mach, 4) - 2.3034E-06 * pow(Alpha, 5);

		CL_RE = -5.119E-04 + 1.0E-03 * Alpha - 1.406E-04 * (Alpha * RE)
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

		CL_LE = -5.119E-04 + 1.0E-03 * Alpha - 1.406E-04 * (Alpha * LE)
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

		CDbv = +1.1457E-02 + CLbv * (-2.4645E-02) + Mach * (0)
			+ (CLbv * Mach) * (4.9698E-02) + pow(CLbv, 2) * (-1.9112e+0)
			+ pow(Mach, 2) * (0) + pow(CLbv * Mach, 2) * (3.5404e+0)
			+ pow(CLbv, 3) * (4.4334e+01) + pow(Mach, 3) * (0)
			+ pow(CLbv * Mach, 3) * (-7.0367e+01)	
			+ pow(CLbv, 4) * (-2.3841e+02) + pow(Mach, 4) * (0)
			+ pow(CLbv * Mach, 4) * (4.1750e+02) + pow(CLbv, 5) * (4.1734e+02)
			+ pow(Mach, 5) * (5.4910E-02)
			+ pow(CLbv * Mach, 5) * (-7.9055e+02);

		CD_LE = -5.184E-04 + 1.100E-03 * Alpha + 3.38E-07 * (Alpha * LE)
			- 1.36E-03 * (Alpha * Mach) - 2.79E-04 * (Mach * LE)
			- 1.53E-04 * (Alpha * Mach) * LE + 1.29E-03 * pow(Mach, 2)
			- 1.02E-04 * pow(Alpha, 2) + 9.39E-08 * pow(LE, 2)
			- 5.69E-07 * pow(Alpha * Mach * LE, 2) + 4.14E-07 * pow(Alpha * LE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * LE, 2)
			- 1.84E-06 * pow(LE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(LE, 4)
			+ 5.71E-09 * pow(LE, 5) - 8.93E-15 * pow(Alpha * Mach * LE, 5)
			- 7.58E-12 * pow(Alpha * Mach * LE, 4) - 3.94E-10 * pow(Alpha * Mach * LE, 3);

		CD_RE = -5.184E-04 + 1.100E-03 * Alpha + 3.38E-07 * (Alpha * RE)
			- 1.36E-03 * (Alpha * Mach) - 2.79E-04 * (Mach * RE)
			- 1.53E-04 * (Alpha * Mach) * RE + 1.29E-03 * pow(Mach, 2)
			- 1.02E-04 * pow(Alpha, 2) + 9.39E-08 * pow(RE, 2)
			- 5.69E-07 * pow(Alpha * Mach * RE, 2) + 4.14E-07 * pow(Alpha * RE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * RE, 2)
			- 1.84E-06 * pow(RE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(RE, 4)
			+ 5.71E-09 * pow(RE, 5) - 8.93E-15 * pow(Alpha * Mach * RE, 5)
			- 7.58E-12 * pow(Alpha * Mach * RE, 4) - 3.94E-10 * pow(Alpha * Mach * RE, 3);

		CD_RUD = +2.47E-04 - 1.93E-04 * Alpha + 7.27E-05 * (Alpha * Mach)
			+ 4.73E-05 * pow(Mach, 2) + 1.50E-05 * pow(Alpha, 2) + 5.03E-06 * pow(RUD, 2)
			- 1.30E-07 * pow(Alpha * Mach * RUD, 2) - 3.50E-08 * pow(Alpha * RUD, 2)
			- 1.68E-06 * pow(Alpha * Mach, 2) + 4.53E-06 * pow(Mach * RUD, 2)
			- 1.98E-11 * pow(Alpha, 3) - 2.63E-08 * pow(Alpha, 4) + 7.54E-09 * pow(RUD, 4)
			+ 3.12E-12 * pow(Alpha * Mach * RUD, 4);

		CYB = -4.750E-01 - 5.000E-02 * Mach;

		CY_LE = -(-1.845E-04 * Mach - 2.13E-07 * (Alpha * LE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * LE)
			+ 6.17E-08 * (Alpha * Mach) * LE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(LE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(LE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * LE, 2) - 2.66E-09 * pow(Alpha * LE, 2)
			+ 3.92E-08 * pow(Mach * LE, 2) + 5.42E-11 * pow(Alpha * Mach * LE, 3)
			- 4.73E-10 * pow(LE, 4) + 7.35E-14 * pow(Alpha * Mach * LE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(LE, 5)
			- 1.11E-15 * pow(Alpha * Mach * LE, 5));

		CY_RE = -1.845E-04 * Mach - 2.13E-07 * (Alpha * RE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * RE)
			+ 6.17E-08 * (Alpha * Mach) * RE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(RE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(RE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * RE, 2) - 2.66E-09 * pow(Alpha * RE, 2)
			+ 3.92E-08 * pow(Mach * RE, 2) + 5.42E-11 * pow(Alpha * Mach * RE, 3)
			- 4.73E-10 * pow(RE, 4) + 7.35E-14 * pow(Alpha * Mach * RE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(RE, 5)
			- 1.11E-15 * pow(Alpha * Mach * RE, 5);

		CY_RUD = +2.440E-03 * RUD;

		Cllbv = -9.380E-02 - 1.250E-02 * Mach;

		Cll_LE = -(5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * LE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * LE)
			- 8.500E-06 * (Alpha * Mach) * LE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(LE, 2)
			- 5.302E-08 * pow(Alpha * Mach * LE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * LE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(LE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * LE, 3) - 3.576E-08 * pow(Alpha * LE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(LE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * LE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(LE, 5) - 2.984E-14 * pow(Alpha * Mach * LE, 5));

		Cll_RE = 5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * RE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * RE)
			- 8.500E-06 * (Alpha * Mach) * RE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(RE, 2)
			- 5.302E-08 * pow(Alpha * Mach * RE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * RE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(RE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * RE, 3) - 3.576E-08 * pow(Alpha * RE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(RE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * RE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(RE, 5) - 2.984E-14 * pow(Alpha * Mach * RE, 5);

		Cll_RUD = +7.0E-04 * RUD;
		Cllr = +2.6250E-01 + 2.50E-02 * (Mach);
		Cllp = -1.33750E-01 - 1.250E-02 * (Mach);

		Cnbv = +1.062E-01 + 6.250E-02 * Mach;

		Cn_LE = -(-2.7E-7 * (Alpha * LE) - 1.008E-05 * (Mach * LE)
			+ 3.564E-07 * (Alpha * Mach) * LE + 1.1E-7 * pow(LE, 3) + 1.11E-07 * pow(LE, 3)
			- 9.32E-12 * pow(Alpha * Mach * LE, 3) - 1.9910E-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(LE, 4) + 1.82E-28 * pow(Alpha * Mach * LE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046E-010 * pow(LE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * LE, 5));

		Cn_RE = -2.7E-7 * (Alpha * RE) - 1.008E-05 * (Mach * RE)
			+ 3.564E-07 * (Alpha * Mach) * RE + 1.1E-7 * pow(RE, 3) + 1.11E-07 * pow(RE, 3)
			- 9.32E-12 * pow(Alpha * Mach * RE, 3) - 1.9910E-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(RE, 4) + 1.82E-28 * pow(Alpha * Mach * RE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046E-010 * pow(RE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * RE, 5);

		Cn_RUD = -3.000E-03 * RUD;
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

		Cm_RE = +2.88E-04 - 5.351E-04 * Alpha + 4.55E-05 * (Alpha * RE)
			+ 3.379E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * RE)
			- 2.770E-05 * (Alpha * Mach) * RE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(RE, 2)
			- 1.000E-07 * pow(Alpha * Mach * RE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * RE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(RE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * RE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(RE, 4) + 1.29E-12 * pow(Alpha * Mach * RE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(RE, 5)
			- 4.67E-14 * pow(Alpha * Mach * RE, 5);

		Cm_LE = +2.88E-04 - 5.351E-04 * Alpha + 4.55E-05 * (Alpha * LE)
			+ 3.379E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * LE)
			- 2.770E-05 * (Alpha * Mach) * LE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(LE, 2)
			- 1.000E-07 * pow(Alpha * Mach * LE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * LE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(LE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * LE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(LE, 4) + 1.29E-12 * pow(Alpha * Mach * LE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(LE, 5)
			- 4.67E-14 * pow(Alpha * Mach * LE, 5);

		Cm_RUD = -1.841E-04 + 3.5E-06 * Alpha + 2.762E-04 * Mach - 1.0E-07 * RUD
			- 4.0E-07 * pow(Alpha, 2) + 5.8E-06 * pow(RUD, 2)
			+ 6.482E-09 * pow(Alpha * Mach * RUD, 2);

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

		CL_RE = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.3093E-005) + (Alpha * Mach) * (0)
			+ (Mach * RE) * (-1.4287E-004)
			+ ((Alpha * Mach) * RE) * (6.1071E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (2.7242E-004)
			+ pow(Alpha * Mach * RE, 2) * (-9.1890E-008)
			+ pow(Alpha * RE, 2) * (3.4060E-007)
			+ pow(Alpha * Mach, 2) * (-6.5093E-006)
			+ pow(Mach * RE, 2) * (-6.3863E-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092E-004)
			+ pow(RE, 3) * (3.8067E-006)
			+ pow(Alpha * Mach * RE, 3) * (2.3165E-011)
			+ pow(Mach, 4) * (-1.0680E-003)
			+ pow(Alpha, 4) * (-2.1893E-005) + pow(RE, 4) * (-3.7716E-007)
			+ pow(Alpha * Mach * RE, 4) * (7.9006E-014)
			+ pow(Mach, 5) * (2.6056E-004)
			+ pow(Alpha, 5) * (9.2099E-007) + pow(RE, 5) * (-8.5345E-009)
			+ pow(Alpha * Mach * RE, 5) * (-2.5698E-017);

		CL_LE = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.3093E-005) + (Alpha * Mach) * (0)
			+ (Mach * LE) * (-1.4287E-004)
			+ ((Alpha * Mach) * LE) * (6.1071E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (2.7242E-004)
			+ pow(Alpha * Mach * LE, 2) * (-9.1890E-008)
			+ pow(Alpha * LE, 2) * (3.4060E-007)
			+ pow(Alpha * Mach, 2) * (-6.5093E-006)
			+ pow(Mach * LE, 2) * (-6.3863E-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092E-004)
			+ pow(LE, 3) * (3.8067E-006)
			+ pow(Alpha * Mach * LE, 3) * (2.3165E-011)
			+ pow(Mach, 4) * (-1.0680E-003)
			+ pow(Alpha, 4) * (-2.1893E-005) + pow(LE, 4) * (-3.7716E-007)
			+ pow(Alpha * Mach * LE, 4) * (7.9006E-014)
			+ pow(Mach, 5) * (2.6056E-004)
			+ pow(Alpha, 5) * (9.2099E-007) + pow(LE, 5) * (-8.5345E-009)
			+ pow(Alpha * Mach * LE, 5) * (-2.5698E-017);

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

		CD_RE = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.6923E-005) + (Alpha * Mach) * (1.5100E-005)
			+ (Mach * RE) * (1.3641E-007)
			+ ((Alpha * Mach) * RE) * (5.1142E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (1.2125E-005)
			+ pow(Alpha * Mach * RE, 2) * (3.5662E-009)
			+ pow(Alpha * RE, 2) * (-1.3848E-008)
			+ pow(Alpha * Mach, 2) * (-4.7972E-007)
			+ pow(Mach * RE, 2) * (-3.3763E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045E-008)
			+ pow(RE, 3) * (3.9119E-008)
			+ pow(Alpha * Mach * RE, 3) * (-9.7714E-013)
			+ pow(Mach, 4) * (9.6475E-007)
			+ pow(Alpha, 4) * (1.5015E-008) + pow(RE, 4) * (4.5137E-009)
			+ pow(Alpha * Mach * RE, 4) * (-6.6207E-016)
			+ pow(Mach, 5) * (-3.2682E-007)
			+ pow(Alpha, 5) * (-3.5360E-010) + pow(RE, 5) * (-1.1538E-010)
			+ pow(Alpha * Mach * RE, 5) * (4.1917E-019);

		CD_LE = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.6923E-005) + (Alpha * Mach) * (1.5100E-005)
			+ (Mach * LE) * (1.3641E-007)
			+ ((Alpha * Mach) * LE) * (5.1142E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (1.2125E-005)
			+ pow(Alpha * Mach * LE, 2) * (3.5662E-009)
			+ pow(Alpha * LE, 2) * (-1.3848E-008)
			+ pow(Alpha * Mach, 2) * (-4.7972E-007)
			+ pow(Mach * LE, 2) * (-3.3763E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045E-008)
			+ pow(LE, 3) * (3.9119E-008)
			+ pow(Alpha * Mach * LE, 3) * (-9.7714E-013)
			+ pow(Mach, 4) * (9.6475E-007)
			+ pow(Alpha, 4) * (1.5015E-008) + pow(LE, 4) * (4.5137E-009)
			+ pow(Alpha * Mach * LE, 4) * (-6.6207E-016)
			+ pow(Mach, 5) * (-3.2682E-007)
			+ pow(Alpha, 5) * (-3.5360E-010) + pow(LE, 5) * (-1.1538E-010)
			+ pow(Alpha * Mach * LE, 5) * (4.1917E-019);

		CD_RUD = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.6425E-021)
			+ (Alpha * Mach) * (-9.8380E-006)
			+ (Mach * RUD) * (1.8193E-020)
			+ ((Alpha * Mach) * RUD) * (1.0319E-021)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RUD, 2) * (8.7608E-006)
			+ pow(Alpha * Mach * RUD, 2) * (5.4045E-010)
			+ pow(Alpha * RUD, 2) * (-2.8939E-008)
			+ pow(Alpha * Mach, 2) * (2.1842E-007)
			+ pow(Mach * RUD, 2) * (-2.9646E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.0067E-007)
			+ pow(RUD, 3) * (-8.8556E-022)
			+ pow(Alpha * Mach * RUD, 3) * (-5.2022E-027)
			+ pow(Mach, 4) * (1.3388E-006) + pow(Alpha, 4) * (1.6460E-007)
			+ pow(RUD, 4) * (4.6754E-010)
			+ pow(Alpha * Mach * RUD, 4) * (2.6560E-016)
			+ pow(Mach, 5) * (-2.5185E-007)
			+ pow(Alpha, 5) * (-7.2766E-009) + pow(RUD, 5) * (1.5611E-024)
			+ pow(Alpha * Mach * RUD, 5) * (5.4442E-033);

		CYB = +(0) + Mach * (0) + Alpha * (-1.1185E-002)
			+ (Alpha * Mach) * (3.0432E-003) + pow(Mach, 2) * (-3.7586E-001)
			+ pow(Alpha, 2) * (3.4004E-003)
			+ pow(Alpha * Mach * Mach, 2) * (-2.4047E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (3.6104E-007)
			+ pow(Alpha * Mach, 2) * (-8.7176E-005)
			+ pow(Alpha * Mach, 4) * (-5.3622E-010) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8160E-004) + pow(Mach, 4) * (9.4289E-002)
			+ pow(Alpha, 4) * (4.4848E-005) + pow(Mach, 5) * (-1.8384E-002)
			+ pow(Alpha, 5) * (-1.3021E-006);

		CY_LE = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(LE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * LE, 2));

		CY_RE = -1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(RE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * RE, 2);

		CY_RUD = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.0067E-005)
			+ (Alpha * Mach) * (0) + (Mach * RUD) * (-5.7185E-004)
			+ ((Alpha * Mach) * RUD) * (-1.5307E-005) + pow(Mach, 2) * (0)
			+ pow(Alpha, 2) * (0) + pow(RUD, 2) * (1.9243E-019)
			+ pow(Alpha * Mach * RUD, 2) * (2.8011E-022)
			+ pow(Alpha * RUD, 2) * (-2.0404E-021)
			+ pow(Alpha * Mach, 2) * (-1.2673E-020)
			+ pow(Mach * RUD, 2) * (-1.7950E-020)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.9873E-019)
			+ pow(RUD, 3) * (3.2768E-005)
			+ pow(Alpha * Mach * RUD, 3) * (1.2674E-012)
			+ pow(Mach, 4) * (-3.8438E-020)
			+ pow(Alpha, 4) * (1.9239E-019) + pow(RUD, 4) * (7.7275E-023)
			+ pow(Alpha * Mach * RUD, 4) * (-3.2592E-029)
			+ pow(Mach, 5) * (3.1048E-020)
			+ pow(Alpha, 5) * (-9.0794E-021) + pow(RUD, 5) * (-6.5825E-008)
			+ pow(Alpha * Mach * RUD, 5) * (1.2684E-017);

		Cllbv = +(0) + Mach * (0) + Alpha * (5.9211E-004)
			+ (Alpha * Mach) * (-3.1579E-004) + pow(Mach, 2) * (-8.7296E-002)
			+ pow(Alpha, 2) * (-5.7398E-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1037E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (-6.8068E-008)
			+ pow(Alpha * Mach, 2) * (2.0549E-005)
			+ pow(Alpha * Mach, 4) * (3.6561E-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-2.8226E-016) + pow(Mach, 4) * (2.0334E-002)
			+ pow(Alpha, 4) * (1.9013E-007) + pow(Mach, 5) * (-3.7733E-003)
			+ pow(Alpha, 5) * (-9.6648E-019);

		Cll_LE = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2));

		Cll_RE = 3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2);

		Cll_RUD = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		Cllr = +3.82E-01 - 1.06E-01 * Mach
			+ 1.94E-03 * Alpha - 8.15E-05 * (Alpha * Mach)
			+ 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2)
			- 1.02E-03 * pow(Mach, 3) - 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4)
			+ 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		Cllp = +(0) + Mach * (0) + Alpha * (-1.2668E-005)
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

		Cn_RE = +2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(RE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * RE, 2);

		Cn_LE = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(LE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * LE, 2));

		Cn_RUD = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD) - 3.13E-07 * (Alpha * Mach) * RUD;

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

		Cm_RE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE) + 9.72E-08 * (Alpha * Mach) * RE;

		Cm_LE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE) + 9.72E-08 * (Alpha * Mach) * LE;

		Cm_RUD = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

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
		CLbv = -8.19E-02 + 4.70E-02 * Mach + 1.86E-02 * Alpha
			- 4.73E-04 * Alpha * Mach - 9.19E-03 * pow(Mach, 2) - 1.52E-04 * pow(Alpha, 2)
			+ 5.99E-07 * pow(Alpha * Mach, 2) + 7.74E-04 * pow(Mach, 3)
			+ 4.08E-06 * pow(Alpha, 3) - 2.93E-05 * pow(Mach, 4) - 3.91E-07 * pow(Alpha, 4)
			+ 4.12E-07 * pow(Mach, 5) + 1.3E-08 * pow(Alpha, 5);

		CL_RE = -1.45E-05 + 1.01E-04 * Alpha + 7.1E-06 * Mach - 4.14E-04 * RE
			- 3.51E-06 * (Alpha * RE) + 4.7E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * RE) - 1.7E-07 * (Alpha * Mach * RE);

		CL_LE = -1.45E-05 + 1.01E-04 * Alpha + 7.1E-06 * Mach - 4.14E-04 * LE
			- 3.51E-06 * (Alpha * LE) + 4.7E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * LE) - 1.7E-07 * (Alpha * Mach * LE);

		CDbv = 8.717E-02 - 3.307E-02 * Mach + 3.179E-03 * Alpha
			- 1.25E-04 * Alpha * Mach + 5.036E-03 * pow(Mach, 2)
			- 1.1E-03 * pow(Alpha, 2) + 1.405E-07 * pow(Alpha * Mach, 2)
			- 3.658E-04 * pow(Mach, 3) + 3.175E-04 * pow(Alpha, 3) + 1.274E-05 * pow(Mach, 4)
			- 2.985E-05 * pow(Alpha, 4) - 1.705E-07 * pow(Mach, 5) + 9.766E-07 * pow(Alpha, 5);

		CD_RE = 4.5548E-04 + 2.5411E-005 * Alpha - 1.1436E-04 * Mach
			+ RE * (-3.6417E-05) + (Alpha * Mach * RE) * (-5.3015E-07)
			+ pow(Alpha, 2) * (3.2187E-06) + pow(Mach, 2) * (3.014E-06)
			+ pow(RE, 2) * (6.9629E-06)
			+ pow(Alpha * Mach * RE, 2) * (2.1026E-12);

		CD_LE = 4.5548E-004 + Alpha * (2.5411E-005) + Mach * (-1.1436E-004)
			+ LE * (-3.6417E-05) + (Alpha * Mach * LE) * (-5.3015E-07)
			+ pow(Alpha, 2) * (3.2187E-06) + pow(Mach, 2) * (3.0140E-06)
			+ pow(LE, 2) * (6.9629E-06)
			+ pow(Alpha * Mach * LE, 2) * (2.1026E-12);

		CD_RUD = 7.5E-04 - 2.29E-05 * Alpha - 9.69E-05 * Mach - 1.83E-06 * RUD
			+ 9.13E-09 * Alpha * Mach * RUD + 8.76E-07 * pow(Alpha, 2)
			+ 2.7E-06 * pow(Mach, 2) + 1.97E-06 * pow(RUD, 2)
			- 1.77E-11 * pow(Alpha * Mach * RUD, 2);

		CYB = Mach * (-2.9253E-01) + Alpha * (2.8803E-03)
			+ (Alpha * Mach) * (-2.8943E-04) + pow(Mach, 2) * (5.4822E-02)
			+ pow(Alpha, 2) * (7.3535E-04)
			+ pow(Alpha * Mach * Mach, 2) * (-4.6490E-09)
			+ pow(Alpha * Alpha * Mach, 2) * (-2.0675E-08)
			+ pow(Alpha * Mach, 2) * (4.6205E-06)
			+ pow(Alpha * Mach, 4) * (2.6144E-11)
			+ pow(Mach, 3) * (-4.3203E-03)
			+ pow(Alpha, 3) * (-3.7405E-04) + pow(Mach, 4) * (1.5495E-04)
			+ pow(Alpha, 4) * (2.8183E-05) + pow(Mach, 5) * (-2.0829E-06)
			+ pow(Alpha, 5) * (-5.2083E-07);

		CY_LE = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(LE, 2) + 4.50E-14 * pow(Alpha * Mach * LE, 2));

		CY_RE = -1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(RE, 2) + 4.50E-14 * pow(Alpha * Mach * RE, 2);

		CY_RUD = -1.43E-18 + 4.86E-20 * Alpha + 1.86E-19 * Mach + 3.84E-04 * RUD
			- 1.17E-05 * (Alpha * RUD) - 1.07E-05 * (Mach * RUD)
			+ 2.60E-07 * (Alpha * Mach) * RUD;

		Cllbv = -1.402E-01 + 3.326E-02 * Mach - 7.590E-04 * Alpha
			+ 8.596E-06 * (Alpha * Mach) - 3.794E-03 * pow(Mach, 2)
			+ 2.354E-06 * pow(Alpha, 2) - 1.044E-08 * pow(Alpha * Mach, 2)
			+ 2.219E-04 * pow(Mach, 3) - 8.964E-18 * pow(Alpha, 3) - 6.462E-06 * pow(Mach, 4)
			+ 3.803E-19 * pow(Alpha, 4) + 7.419E-08 * pow(Mach, 5) - 3.353E-21 * pow(Alpha, 5);

		Cll_RE = +3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2);

		Cll_LE = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2));

		Cll_RUD = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		Cllr = +3.82E-01 - 1.06E-01 * Mach + 1.94E-03 * Alpha
			- 8.15E-05 * (Alpha * Mach) + 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2) - 1.02E-03 * pow(Mach, 3)
			- 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4) + 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		Cllp = -2.99E-01 + 7.47E-02 * Mach + 1.38E-03 * Alpha
			- 8.78E-05 * (Alpha * Mach) - 9.13E-03 * pow(Mach, 2) - 2.04E-04 * pow(Alpha, 2)
			- 1.52E-07 * pow(Alpha * Mach, 2) + 5.73E-04 * pow(Mach, 3)
			- 3.86E-05 * pow(Alpha, 3) - 1.79E-05 * pow(Mach, 4) + 4.21E-06 * pow(Alpha, 4)
			+ 2.20E-07 * pow(Mach, 5) - 1.15E-07 * pow(Alpha, 5);

		Cnbv = +(0) + Alpha * (6.9980E-004) + Mach * (5.9115E-002)
			+ (Alpha * Mach) * (-7.5250E-005) + pow(Alpha, 2) * (2.5160E-004)
			+ pow(Mach, 2) * (-1.4824E-002)
			+ pow(Alpha * Mach, 2) * (-2.1924E-007)
			+ pow(Alpha, 3) * (-1.0777E-004) + pow(Mach, 3) * (1.2692E-003)
			+ pow(Alpha * Mach, 3) * (1.0707E-008)
			+ pow(Alpha, 4) * (9.4989E-006) + pow(Mach, 4) * (-4.7098E-005)
			+ pow(Alpha * Mach, 4) * (-5.5472E-011)
			+ pow(Alpha, 5) * (-2.5953E-007) + pow(Mach, 5) * (6.4284E-007)
			+ pow(Alpha * Mach, 5) * (8.5863E-014);

		Cn_LE = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(LE, 2) + 1.41E-11 * pow(Alpha * Mach * LE, 2));

		Cn_RE = 2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(RE, 2) + 1.41E-11 * pow(Alpha * Mach * RE, 2);

		Cn_RUD = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD)
			- 3.13E-07 * (Alpha * Mach) * RUD;

		Cnp = +3.68E-01 - 9.79E-02 * Mach + 7.61E-16 * Alpha + 1.24E-02 * pow(Mach, 2)
			- 4.64E-16 * pow(Alpha, 2) - 8.05E-04 * pow(Mach, 3) + 1.01E-16 * pow(Alpha, 3)
			+ 2.57E-05 * pow(Mach, 4)
			- 9.18E-18 * pow(Alpha, 4) - 3.20E-07 * pow(Mach, 5) + 2.96E-19 * pow(Alpha, 5);

		Cnr = -2.41E+00 + 5.96E-01 * Mach - 2.74E-03 * Alpha
			+ 2.09E-04 * (Alpha * Mach) - 7.57E-02 * pow(Mach, 2)
			+ 1.15E-03 * pow(Alpha, 2) - 6.53E-08 * pow(Alpha * Mach, 2)
			+ 4.90E-03 * pow(Mach, 3) - 3.87E-04 * pow(Alpha, 3) - 1.57E-04 * pow(Mach, 4)
			+ 3.60E-05 * pow(Alpha, 4) + 1.96E-06 * pow(Mach, 5) - 1.18E-06 * pow(Alpha, 5);

		Cmbv = -2.192E-02 + 7.739E-03 * Mach - 2.260E-03 * Alpha
			+ 1.808E-04 * (Alpha * Mach) - 8.849E-04 * pow(Mach, 2)
			+ 2.616E-04 * pow(Alpha, 2) - 2.880E-07 * pow(Alpha * Mach, 2)
			+ 4.617E-05 * pow(Mach, 3) - 7.887E-05 * pow(Alpha, 3) - 1.143E-06 * pow(Mach, 4)
			+ 8.288E-06 * pow(Alpha, 4) + 1.082E-08 * pow(Mach, 5) - 2.789E-07 * pow(Alpha, 5);

		Cm_LE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE)
			+ 9.72E-08 * (Alpha * Mach) * LE;

		Cm_RE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE)
			+ 9.72E-08 * (Alpha * Mach) * RE;

		Cm_RUD = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

		Cm_q = -1.36E+00 + 3.86E-01 * Mach + 7.85E-04 * Alpha
			+ 1.40E-04 * (Alpha * Mach) - 5.42E-02 * pow(Mach, 2)
			+ 2.36E-03 * pow(Alpha, 2) - 1.95E-06 * pow(Alpha * Mach, 2)
			+ 3.80E-03 * pow(Mach, 3) - 1.48E-03 * pow(Alpha, 3) - 1.30E-04 * pow(Mach, 4)
			+ 1.69E-04 * pow(Alpha, 4) + 1.71E-06 * pow(Mach, 5) - 5.93E-06 * pow(Alpha, 5);
	}

	// drag force 牵引力 x轴
	double CD = CDbv + CD_RE + CD_LE + CD_RUD;
	// lift force 升力 z轴
	double CL = CLbv + CL_RE + CL_LE;
	// side force 侧向力 y轴
	double CN = CYB * Beta + CY_RE + CY_LE + CY_RUD;
	// rolling moment 滚转力矩 绕x轴 改变滚转角速度p
	double mx = Cllbv * Beta + Cll_RE + Cll_LE + Cll_RUD + (Cllp * omega_x + Cllr * omega_z) * B / (2 * Vel);
	// pitching moment 俯仰力矩 绕y轴 改变俯仰角速度q
	double my = Cmbv + Cm_RE + Cm_LE + Cm_RUD + Cm_q * omega_y * C / (2 * Vel); 
	// yawing moment 偏航力矩 绕z轴 改变偏航加速度r
	double mz = Cnbv * Beta + Cn_RE + Cn_LE + Cn_RUD + (Cnp * omega_x + Cnr * omega_z) * B / (2 * Vel);

	if (flag == 1) {
		srand((unsigned int)time(NULL));
		/*double D = (double)(rand() % 21 - 10) / 100, L = (double)(rand() % 21 - 10) / 100, N = (double)(rand() % 21 - 10) / 100;
		double l = (double)(rand() % 21 - 10) / 100, m = (double)(rand() % 21 - 10) / 100, n = (double)(rand() % 21 - 10) / 100;*/
		double D = (double)(rand() % 41 - 20) / 100, L = (double)(rand() % 41 - 20) / 100, N = (double)(rand() % 41 - 20) / 100;
		double l = (double)(rand() % 41 - 20) / 100, m = (double)(rand() % 41 - 20) / 100, n = (double)(rand() % 41 - 20) / 100;
		ans[0] = CD * (1 + D);
		ans[1] = CL * (1 + L);
		ans[2] = CN * (1 + N);
		ans[3] = mx * (1 + l);
		ans[4] = my * (1 + m);
		ans[5] = mz * (1 + n);
	}
	else {
		ans[0] = CD;
		ans[1] = CL;
		ans[2] = CN;
		ans[3] = mx;	// 引起相对于机体坐标系x轴的转动，产生滚转角速度p
		ans[4] = my;	// 引起相对于机体坐标系y轴的转动，产生俯仰角速度q
		ans[5] = mz;	// 引起相对于机体坐标系z轴的转动，产生偏航角速度r
	}

	return;
}

void getCoefficientsDetails(const double& Mach, const double& Alpha, const double& Beta, const double& Delta_e, const double& Delta_a, const double& Delta_r, const double& C, const double& B, const double& Vel, const double& omega_x, const double& omega_y, const double& omega_z, std::vector<std::vector<double>>& ans)
{
	double LE = Delta_e;		// Left Elevon 左襟翼
	double RE = Delta_a;		// Right Elevon 右襟翼
	double RUD = Delta_r;

	double CLbv, CL_RE, CL_LE, CDbv, CD_RE, CD_LE, CD_RUD, CYB, CY_RE, CY_LE, CY_RUD;
	double Cllbv, Cll_RE, Cll_LE, Cll_RUD, Cllr, Cllp, Cnbv, Cn_RE, Cn_LE, Cn_RUD, Cnp, Cnr, Cmbv, Cm_RE, Cm_LE, Cm_RUD, Cm_q;
	if (Mach < 0 || Mach > 30) {
		std::cout << Delta_e << " " << Mach << "坠毁了" << std::endl;
		exit(1);
	}
	else if (Mach <= 1.25) {
		CLbv = -5.2491E-04 + Alpha * 1.5746E-02 + (Alpha * Mach) * 6.0213E-03
			- 3.4437E-04 * pow(Alpha, 2) + pow(Alpha * Mach, 2) * 1.4471E-04
			- 5.1952E-05 * pow(Alpha, 3) + 3.4771E-05 * pow(Alpha, 4)
			+ 2.7717E-03 * pow(Mach, 4) - 2.3034E-06 * pow(Alpha, 5);

		CL_RE = -5.119E-04 + 1.0E-03 * Alpha - 1.406E-04 * (Alpha * RE)
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

		CL_LE = -5.119E-04 + 1.0E-03 * Alpha - 1.406E-04 * (Alpha * LE)
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

		CDbv = +1.1457E-02 + CLbv * (-2.4645E-02) + Mach * (0)
			+ (CLbv * Mach) * (4.9698E-02) + pow(CLbv, 2) * (-1.9112e+0)
			+ pow(Mach, 2) * (0) + pow(CLbv * Mach, 2) * (3.5404e+0)
			+ pow(CLbv, 3) * (4.4334e+01) + pow(Mach, 3) * (0)
			+ pow(CLbv * Mach, 3) * (-7.0367e+01)
			+ pow(CLbv, 4) * (-2.3841e+02) + pow(Mach, 4) * (0)
			+ pow(CLbv * Mach, 4) * (4.1750e+02) + pow(CLbv, 5) * (4.1734e+02)
			+ pow(Mach, 5) * (5.4910E-02)
			+ pow(CLbv * Mach, 5) * (-7.9055e+02);

		CD_LE = -5.184E-04 + 1.100E-03 * Alpha + 3.38E-07 * (Alpha * LE)
			- 1.36E-03 * (Alpha * Mach) - 2.79E-04 * (Mach * LE)
			- 1.53E-04 * (Alpha * Mach) * LE + 1.29E-03 * pow(Mach, 2)
			- 1.02E-04 * pow(Alpha, 2) + 9.39E-08 * pow(LE, 2)
			- 5.69E-07 * pow(Alpha * Mach * LE, 2) + 4.14E-07 * pow(Alpha * LE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * LE, 2)
			- 1.84E-06 * pow(LE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(LE, 4)
			+ 5.71E-09 * pow(LE, 5) - 8.93E-15 * pow(Alpha * Mach * LE, 5)
			- 7.58E-12 * pow(Alpha * Mach * LE, 4) - 3.94E-10 * pow(Alpha * Mach * LE, 3);

		CD_RE = -5.184E-04 + 1.100E-03 * Alpha + 3.38E-07 * (Alpha * RE)
			- 1.36E-03 * (Alpha * Mach) - 2.79E-04 * (Mach * RE)
			- 1.53E-04 * (Alpha * Mach) * RE + 1.29E-03 * pow(Mach, 2)
			- 1.02E-04 * pow(Alpha, 2) + 9.39E-08 * pow(RE, 2)
			- 5.69E-07 * pow(Alpha * Mach * RE, 2) + 4.14E-07 * pow(Alpha * RE, 2)
			+ 1.81E-04 * pow(Alpha * Mach, 2) - 1.68E-05 * pow(Mach * RE, 2)
			- 1.84E-06 * pow(RE, 3) + 6.40E-08 * pow(Alpha, 4) + 5.76E-08 * pow(RE, 4)
			+ 5.71E-09 * pow(RE, 5) - 8.93E-15 * pow(Alpha * Mach * RE, 5)
			- 7.58E-12 * pow(Alpha * Mach * RE, 4) - 3.94E-10 * pow(Alpha * Mach * RE, 3);

		CD_RUD = +2.47E-04 - 1.93E-04 * Alpha + 7.27E-05 * (Alpha * Mach)
			+ 4.73E-05 * pow(Mach, 2) + 1.50E-05 * pow(Alpha, 2) + 5.03E-06 * pow(RUD, 2)
			- 1.30E-07 * pow(Alpha * Mach * RUD, 2) - 3.50E-08 * pow(Alpha * RUD, 2)
			- 1.68E-06 * pow(Alpha * Mach, 2) + 4.53E-06 * pow(Mach * RUD, 2)
			- 1.98E-11 * pow(Alpha, 3) - 2.63E-08 * pow(Alpha, 4) + 7.54E-09 * pow(RUD, 4)
			+ 3.12E-12 * pow(Alpha * Mach * RUD, 4);

		CYB = -4.750E-01 - 5.000E-02 * Mach;

		CY_LE = -(-1.845E-04 * Mach - 2.13E-07 * (Alpha * LE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * LE)
			+ 6.17E-08 * (Alpha * Mach) * LE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(LE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(LE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * LE, 2) - 2.66E-09 * pow(Alpha * LE, 2)
			+ 3.92E-08 * pow(Mach * LE, 2) + 5.42E-11 * pow(Alpha * Mach * LE, 3)
			- 4.73E-10 * pow(LE, 4) + 7.35E-14 * pow(Alpha * Mach * LE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(LE, 5)
			- 1.11E-15 * pow(Alpha * Mach * LE, 5));

		CY_RE = -1.845E-04 * Mach - 2.13E-07 * (Alpha * RE)
			+ 3.740E-05 * (Alpha * Mach) + 1.990E-05 * (Mach * RE)
			+ 6.17E-08 * (Alpha * Mach) * RE + 3.39E-06 * pow(Alpha, 2)
			+ 1.37E-07 * pow(RE, 2) - 2.14E-06 * pow(Alpha * Mach, 2) - 1.11E-06 * pow(Alpha, 3)
			- 3.40E-07 * pow(RE, 3) + 1.09E-07 * pow(Alpha, 4)
			+ 3.53E-09 * pow(Alpha * Mach * RE, 2) - 2.66E-09 * pow(Alpha * RE, 2)
			+ 3.92E-08 * pow(Mach * RE, 2) + 5.42E-11 * pow(Alpha * Mach * RE, 3)
			- 4.73E-10 * pow(RE, 4) + 7.35E-14 * pow(Alpha * Mach * RE, 4)
			- 3.45E-09 * pow(Alpha, 5) + 6.53E-10 * pow(RE, 5)
			- 1.11E-15 * pow(Alpha * Mach * RE, 5);

		CY_RUD = +2.440E-03 * RUD;

		Cllbv = -9.380E-02 - 1.250E-02 * Mach;

		Cll_LE = -(5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * LE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * LE)
			- 8.500E-06 * (Alpha * Mach) * LE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(LE, 2)
			- 5.302E-08 * pow(Alpha * Mach * LE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * LE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(LE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * LE, 3) - 3.576E-08 * pow(Alpha * LE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(LE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * LE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(LE, 5) - 2.984E-14 * pow(Alpha * Mach * LE, 5));

		Cll_RE = 5.310E-05 - 5.272E-04 * Alpha + 3.690E-05 * (Alpha * RE)
			+ 2.680E-05 * (Alpha * Mach) + 1.926E-04 * (Mach * RE)
			- 8.500E-06 * (Alpha * Mach) * RE - 4.097E-04 * pow(Mach, 2)
			+ 1.258E-04 * pow(Alpha, 2) + 3.762E-06 * pow(RE, 2)
			- 5.302E-08 * pow(Alpha * Mach * RE, 2) + 5.100E-06 * pow(Alpha * Mach, 2)
			+ 2.100E-06 * pow(Mach * RE, 2) - 8.700E-06 * pow(Alpha, 3) + 8.400E-06 * pow(RE, 3)
			+ 1.153E-09 * pow(Alpha * Mach * RE, 3) - 3.576E-08 * pow(Alpha * RE, 2)
			+ 1.384E-08 * pow(Alpha, 4) - 1.137E-08 * pow(RE, 4)
			+ 1.011E-12 * pow(Alpha * Mach * RE, 4) + 1.381E-08 * pow(Alpha, 5)
			- 1.676E-08 * pow(RE, 5) - 2.984E-14 * pow(Alpha * Mach * RE, 5);

		Cll_RUD = +7.0E-04 * RUD;
		Cllr = +2.6250E-01 + 2.50E-02 * (Mach);
		Cllp = -1.33750E-01 - 1.250E-02 * (Mach);

		Cnbv = +1.062E-01 + 6.250E-02 * Mach;

		Cn_LE = -(-2.7E-7 * (Alpha * LE) - 1.008E-05 * (Mach * LE)
			+ 3.564E-07 * (Alpha * Mach) * LE + 1.1E-7 * pow(LE, 3) + 1.11E-07 * pow(LE, 3)
			- 9.32E-12 * pow(Alpha * Mach * LE, 3) - 1.9910E-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(LE, 4) + 1.82E-28 * pow(Alpha * Mach * LE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046E-010 * pow(LE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * LE, 5));

		Cn_RE = -2.7E-7 * (Alpha * RE) - 1.008E-05 * (Mach * RE)
			+ 3.564E-07 * (Alpha * Mach) * RE + 1.1E-7 * pow(RE, 3) + 1.11E-07 * pow(RE, 3)
			- 9.32E-12 * pow(Alpha * Mach * RE, 3) - 1.9910E-021 * pow(Alpha, 4)
			+ 2.89E-25 * pow(RE, 4) + 1.82E-28 * pow(Alpha * Mach * RE, 4)
			+ 6.95E-23 * pow(Alpha, 5) - 2.2046E-010 * pow(RE, 5)
			+ 2.22E-16 * pow(Alpha * Mach * RE, 5);

		Cn_RUD = -3.000E-03 * RUD;
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

		Cm_RE = +2.88E-04 - 5.351E-04 * Alpha + 4.55E-05 * (Alpha * RE)
			+ 3.379E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * RE)
			- 2.770E-05 * (Alpha * Mach) * RE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(RE, 2)
			- 1.000E-07 * pow(Alpha * Mach * RE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * RE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(RE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * RE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(RE, 4) + 1.29E-12 * pow(Alpha * Mach * RE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(RE, 5)
			- 4.67E-14 * pow(Alpha * Mach * RE, 5);

		Cm_LE = +2.88E-04 - 5.351E-04 * Alpha + 4.55E-05 * (Alpha * LE)
			+ 3.379E-04 * (Alpha * Mach) + 6.665E-04 * (Mach * LE)
			- 2.770E-05 * (Alpha * Mach) * LE - 6.027E-04 * pow(Mach, 2)
			+ 2.660E-05 * pow(Alpha, 2) - 1.600E-06 * pow(LE, 2)
			- 1.000E-07 * pow(Alpha * Mach * LE, 2) - 1.910E-05 * pow(Alpha * Mach, 2)
			+ 2.300E-06 * pow(Mach * LE, 2) + 1.300E-05 * pow(Alpha, 3) + 1.920E-05 * pow(LE, 3)
			+ 1.90E-09 * pow(Alpha * Mach * LE, 3) - 1.861200E-06 * pow(Alpha, 4)
			- 4.69E-10 * pow(LE, 4) + 1.29E-12 * pow(Alpha * Mach * LE, 4)
			+ 7.29E-08 * pow(Alpha, 5) - 3.87E-08 * pow(LE, 5)
			- 4.67E-14 * pow(Alpha * Mach * LE, 5);

		Cm_RUD = -1.841E-04 + 3.5E-06 * Alpha + 2.762E-04 * Mach - 1.0E-07 * RUD
			- 4.0E-07 * pow(Alpha, 2) + 5.8E-06 * pow(RUD, 2)
			+ 6.482E-09 * pow(Alpha * Mach * RUD, 2);

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

		CL_RE = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.3093E-005) + (Alpha * Mach) * (0)
			+ (Mach * RE) * (-1.4287E-004)
			+ ((Alpha * Mach) * RE) * (6.1071E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (2.7242E-004)
			+ pow(Alpha * Mach * RE, 2) * (-9.1890E-008)
			+ pow(Alpha * RE, 2) * (3.4060E-007)
			+ pow(Alpha * Mach, 2) * (-6.5093E-006)
			+ pow(Mach * RE, 2) * (-6.3863E-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092E-004)
			+ pow(RE, 3) * (3.8067E-006)
			+ pow(Alpha * Mach * RE, 3) * (2.3165E-011)
			+ pow(Mach, 4) * (-1.0680E-003)
			+ pow(Alpha, 4) * (-2.1893E-005) + pow(RE, 4) * (-3.7716E-007)
			+ pow(Alpha * Mach * RE, 4) * (7.9006E-014)
			+ pow(Mach, 5) * (2.6056E-004)
			+ pow(Alpha, 5) * (9.2099E-007) + pow(RE, 5) * (-8.5345E-009)
			+ pow(Alpha * Mach * RE, 5) * (-2.5698E-017);

		CL_LE = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.3093E-005) + (Alpha * Mach) * (0)
			+ (Mach * LE) * (-1.4287E-004)
			+ ((Alpha * Mach) * LE) * (6.1071E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (2.7242E-004)
			+ pow(Alpha * Mach * LE, 2) * (-9.1890E-008)
			+ pow(Alpha * LE, 2) * (3.4060E-007)
			+ pow(Alpha * Mach, 2) * (-6.5093E-006)
			+ pow(Mach * LE, 2) * (-6.3863E-006)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (1.4092E-004)
			+ pow(LE, 3) * (3.8067E-006)
			+ pow(Alpha * Mach * LE, 3) * (2.3165E-011)
			+ pow(Mach, 4) * (-1.0680E-003)
			+ pow(Alpha, 4) * (-2.1893E-005) + pow(LE, 4) * (-3.7716E-007)
			+ pow(Alpha * Mach * LE, 4) * (7.9006E-014)
			+ pow(Mach, 5) * (2.6056E-004)
			+ pow(Alpha, 5) * (9.2099E-007) + pow(LE, 5) * (-8.5345E-009)
			+ pow(Alpha * Mach * LE, 5) * (-2.5698E-017);

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

		CD_RE = +(0) * 1 + Mach * (0) + Alpha * (0) + RE * (0)
			+ (Alpha * RE) * (-3.6923E-005) + (Alpha * Mach) * (1.5100E-005)
			+ (Mach * RE) * (1.3641E-007)
			+ ((Alpha * Mach) * RE) * (5.1142E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RE, 2) * (1.2125E-005)
			+ pow(Alpha * Mach * RE, 2) * (3.5662E-009)
			+ pow(Alpha * RE, 2) * (-1.3848E-008)
			+ pow(Alpha * Mach, 2) * (-4.7972E-007)
			+ pow(Mach * RE, 2) * (-3.3763E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045E-008)
			+ pow(RE, 3) * (3.9119E-008)
			+ pow(Alpha * Mach * RE, 3) * (-9.7714E-013)
			+ pow(Mach, 4) * (9.6475E-007)
			+ pow(Alpha, 4) * (1.5015E-008) + pow(RE, 4) * (4.5137E-009)
			+ pow(Alpha * Mach * RE, 4) * (-6.6207E-016)
			+ pow(Mach, 5) * (-3.2682E-007)
			+ pow(Alpha, 5) * (-3.5360E-010) + pow(RE, 5) * (-1.1538E-010)
			+ pow(Alpha * Mach * RE, 5) * (4.1917E-019);

		CD_LE = +(0) * 1 + Mach * (0) + Alpha * (0) + LE * (0)
			+ (Alpha * LE) * (-3.6923E-005) + (Alpha * Mach) * (1.5100E-005)
			+ (Mach * LE) * (1.3641E-007)
			+ ((Alpha * Mach) * LE) * (5.1142E-006)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(LE, 2) * (1.2125E-005)
			+ pow(Alpha * Mach * LE, 2) * (3.5662E-009)
			+ pow(Alpha * LE, 2) * (-1.3848E-008)
			+ pow(Alpha * Mach, 2) * (-4.7972E-007)
			+ pow(Mach * LE, 2) * (-3.3763E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-4.6045E-008)
			+ pow(LE, 3) * (3.9119E-008)
			+ pow(Alpha * Mach * LE, 3) * (-9.7714E-013)
			+ pow(Mach, 4) * (9.6475E-007)
			+ pow(Alpha, 4) * (1.5015E-008) + pow(LE, 4) * (4.5137E-009)
			+ pow(Alpha * Mach * LE, 4) * (-6.6207E-016)
			+ pow(Mach, 5) * (-3.2682E-007)
			+ pow(Alpha, 5) * (-3.5360E-010) + pow(LE, 5) * (-1.1538E-010)
			+ pow(Alpha * Mach * LE, 5) * (4.1917E-019);

		CD_RUD = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.6425E-021)
			+ (Alpha * Mach) * (-9.8380E-006)
			+ (Mach * RUD) * (1.8193E-020)
			+ ((Alpha * Mach) * RUD) * (1.0319E-021)
			+ pow(Mach, 2) * (0) + pow(Alpha, 2) * (0) + pow(RUD, 2) * (8.7608E-006)
			+ pow(Alpha * Mach * RUD, 2) * (5.4045E-010)
			+ pow(Alpha * RUD, 2) * (-2.8939E-008)
			+ pow(Alpha * Mach, 2) * (2.1842E-007)
			+ pow(Mach * RUD, 2) * (-2.9646E-007)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.0067E-007)
			+ pow(RUD, 3) * (-8.8556E-022)
			+ pow(Alpha * Mach * RUD, 3) * (-5.2022E-027)
			+ pow(Mach, 4) * (1.3388E-006) + pow(Alpha, 4) * (1.6460E-007)
			+ pow(RUD, 4) * (4.6754E-010)
			+ pow(Alpha * Mach * RUD, 4) * (2.6560E-016)
			+ pow(Mach, 5) * (-2.5185E-007)
			+ pow(Alpha, 5) * (-7.2766E-009) + pow(RUD, 5) * (1.5611E-024)
			+ pow(Alpha * Mach * RUD, 5) * (5.4442E-033);

		CYB = +(0) + Mach * (0) + Alpha * (-1.1185E-002)
			+ (Alpha * Mach) * (3.0432E-003) + pow(Mach, 2) * (-3.7586E-001)
			+ pow(Alpha, 2) * (3.4004E-003)
			+ pow(Alpha * Mach * Mach, 2) * (-2.4047E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (3.6104E-007)
			+ pow(Alpha * Mach, 2) * (-8.7176E-005)
			+ pow(Alpha * Mach, 4) * (-5.3622E-010) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-5.8160E-004) + pow(Mach, 4) * (9.4289E-002)
			+ pow(Alpha, 4) * (4.4848E-005) + pow(Mach, 5) * (-1.8384E-002)
			+ pow(Alpha, 5) * (-1.3021E-006);

		CY_LE = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(LE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * LE, 2));

		CY_RE = -1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2)
			- 2.36E-08 * pow(Mach, 2) - 5.04E-08 * pow(RE, 2)
			+ 4.50E-14 * pow(Alpha * Mach * RE, 2);

		CY_RUD = +(0) * 1 + Mach * (0) + Alpha * (0) + RUD * (0)
			+ (Alpha * RUD) * (2.0067E-005)
			+ (Alpha * Mach) * (0) + (Mach * RUD) * (-5.7185E-004)
			+ ((Alpha * Mach) * RUD) * (-1.5307E-005) + pow(Mach, 2) * (0)
			+ pow(Alpha, 2) * (0) + pow(RUD, 2) * (1.9243E-019)
			+ pow(Alpha * Mach * RUD, 2) * (2.8011E-022)
			+ pow(Alpha * RUD, 2) * (-2.0404E-021)
			+ pow(Alpha * Mach, 2) * (-1.2673E-020)
			+ pow(Mach * RUD, 2) * (-1.7950E-020)
			+ pow(Mach, 3) * (0) + pow(Alpha, 3) * (-9.9873E-019)
			+ pow(RUD, 3) * (3.2768E-005)
			+ pow(Alpha * Mach * RUD, 3) * (1.2674E-012)
			+ pow(Mach, 4) * (-3.8438E-020)
			+ pow(Alpha, 4) * (1.9239E-019) + pow(RUD, 4) * (7.7275E-023)
			+ pow(Alpha * Mach * RUD, 4) * (-3.2592E-029)
			+ pow(Mach, 5) * (3.1048E-020)
			+ pow(Alpha, 5) * (-9.0794E-021) + pow(RUD, 5) * (-6.5825E-008)
			+ pow(Alpha * Mach * RUD, 5) * (1.2684E-017);

		Cllbv = +(0) + Mach * (0) + Alpha * (5.9211E-004)
			+ (Alpha * Mach) * (-3.1579E-004) + pow(Mach, 2) * (-8.7296E-002)
			+ pow(Alpha, 2) * (-5.7398E-005)
			+ pow(Alpha * Mach * Mach, 2) * (-1.1037E-006)
			+ pow(Alpha * Alpha * Mach, 2) * (-6.8068E-008)
			+ pow(Alpha * Mach, 2) * (2.0549E-005)
			+ pow(Alpha * Mach, 4) * (3.6561E-009) + pow(Mach, 3) * (0)
			+ pow(Alpha, 3) * (-2.8226E-016) + pow(Mach, 4) * (2.0334E-002)
			+ pow(Alpha, 4) * (1.9013E-007) + pow(Mach, 5) * (-3.7733E-003)
			+ pow(Alpha, 5) * (-9.6648E-019);

		Cll_LE = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2));

		Cll_RE = 3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2);

		Cll_RUD = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		Cllr = +3.82E-01 - 1.06E-01 * Mach
			+ 1.94E-03 * Alpha - 8.15E-05 * (Alpha * Mach)
			+ 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2)
			- 1.02E-03 * pow(Mach, 3) - 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4)
			+ 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		Cllp = +(0) + Mach * (0) + Alpha * (-1.2668E-005)
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

		Cn_RE = +2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(RE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * RE, 2);

		Cn_LE = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2)
			+ 8.16E-07 * pow(Mach, 2) + 1.97E-06 * pow(LE, 2)
			+ 1.41E-11 * pow(Alpha * Mach * LE, 2));

		Cn_RUD = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD) - 3.13E-07 * (Alpha * Mach) * RUD;

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

		Cm_RE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE) + 9.72E-08 * (Alpha * Mach) * RE;

		Cm_LE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE) + 9.72E-08 * (Alpha * Mach) * LE;

		Cm_RUD = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

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
		CLbv = -8.19E-02 + 4.70E-02 * Mach + 1.86E-02 * Alpha
			- 4.73E-04 * Alpha * Mach - 9.19E-03 * pow(Mach, 2) - 1.52E-04 * pow(Alpha, 2)
			+ 5.99E-07 * pow(Alpha * Mach, 2) + 7.74E-04 * pow(Mach, 3)
			+ 4.08E-06 * pow(Alpha, 3) - 2.93E-05 * pow(Mach, 4) - 3.91E-07 * pow(Alpha, 4)
			+ 4.12E-07 * pow(Mach, 5) + 1.3E-08 * pow(Alpha, 5);

		CL_RE = -1.45E-05 + 1.01E-04 * Alpha + 7.1E-06 * Mach - 4.14E-04 * RE
			- 3.51E-06 * (Alpha * RE) + 4.7E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * RE) - 1.7E-07 * (Alpha * Mach * RE);

		CL_LE = -1.45E-05 + 1.01E-04 * Alpha + 7.1E-06 * Mach - 4.14E-04 * LE
			- 3.51E-06 * (Alpha * LE) + 4.7E-06 * (Alpha * Mach)
			+ 8.72E-06 * (Mach * LE) - 1.7E-07 * (Alpha * Mach * LE);

		CDbv = 8.717E-02 - 3.307E-02 * Mach + 3.179E-03 * Alpha
			- 1.25E-04 * Alpha * Mach + 5.036E-03 * pow(Mach, 2)
			- 1.1E-03 * pow(Alpha, 2) + 1.405E-07 * pow(Alpha * Mach, 2)
			- 3.658E-04 * pow(Mach, 3) + 3.175E-04 * pow(Alpha, 3) + 1.274E-05 * pow(Mach, 4)
			- 2.985E-05 * pow(Alpha, 4) - 1.705E-07 * pow(Mach, 5) + 9.766E-07 * pow(Alpha, 5);

		CD_RE = 4.5548E-04 + 2.5411E-005 * Alpha - 1.1436E-04 * Mach
			+ RE * (-3.6417E-05) + (Alpha * Mach * RE) * (-5.3015E-07)
			+ pow(Alpha, 2) * (3.2187E-06) + pow(Mach, 2) * (3.014E-06)
			+ pow(RE, 2) * (6.9629E-06)
			+ pow(Alpha * Mach * RE, 2) * (2.1026E-12);

		CD_LE = 4.5548E-004 + Alpha * (2.5411E-005) + Mach * (-1.1436E-004)
			+ LE * (-3.6417E-05) + (Alpha * Mach * LE) * (-5.3015E-07)
			+ pow(Alpha, 2) * (3.2187E-06) + pow(Mach, 2) * (3.0140E-06)
			+ pow(LE, 2) * (6.9629E-06)
			+ pow(Alpha * Mach * LE, 2) * (2.1026E-12);

		CD_RUD = 7.5E-04 - 2.29E-05 * Alpha - 9.69E-05 * Mach - 1.83E-06 * RUD
			+ 9.13E-09 * Alpha * Mach * RUD + 8.76E-07 * pow(Alpha, 2)
			+ 2.7E-06 * pow(Mach, 2) + 1.97E-06 * pow(RUD, 2)
			- 1.77E-11 * pow(Alpha * Mach * RUD, 2);

		CYB = Mach * (-2.9253E-01) + Alpha * (2.8803E-03)
			+ (Alpha * Mach) * (-2.8943E-04) + pow(Mach, 2) * (5.4822E-02)
			+ pow(Alpha, 2) * (7.3535E-04)
			+ pow(Alpha * Mach * Mach, 2) * (-4.6490E-09)
			+ pow(Alpha * Alpha * Mach, 2) * (-2.0675E-08)
			+ pow(Alpha * Mach, 2) * (4.6205E-06)
			+ pow(Alpha * Mach, 4) * (2.6144E-11)
			+ pow(Mach, 3) * (-4.3203E-03)
			+ pow(Alpha, 3) * (-3.7405E-04) + pow(Mach, 4) * (1.5495E-04)
			+ pow(Alpha, 4) * (2.8183E-05) + pow(Mach, 5) * (-2.0829E-06)
			+ pow(Alpha, 5) * (-5.2083E-07);

		CY_LE = -(-1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * LE
			+ 4.11E-09 * (Alpha * Mach) * LE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(LE, 2) + 4.50E-14 * pow(Alpha * Mach * LE, 2));

		CY_RE = -1.02E-06 - 1.12E-07 * Alpha + 4.48E-07 * Mach + 2.27E-07 * RE
			+ 4.11E-09 * (Alpha * Mach) * RE + 2.82E-09 * pow(Alpha, 2) - 2.36E-08 * pow(Mach, 2)
			- 5.04E-08 * pow(RE, 2) + 4.50E-14 * pow(Alpha * Mach * RE, 2);

		CY_RUD = -1.43E-18 + 4.86E-20 * Alpha + 1.86E-19 * Mach + 3.84E-04 * RUD
			- 1.17E-05 * (Alpha * RUD) - 1.07E-05 * (Mach * RUD)
			+ 2.60E-07 * (Alpha * Mach) * RUD;

		Cllbv = -1.402E-01 + 3.326E-02 * Mach - 7.590E-04 * Alpha
			+ 8.596E-06 * (Alpha * Mach) - 3.794E-03 * pow(Mach, 2)
			+ 2.354E-06 * pow(Alpha, 2) - 1.044E-08 * pow(Alpha * Mach, 2)
			+ 2.219E-04 * pow(Mach, 3) - 8.964E-18 * pow(Alpha, 3) - 6.462E-06 * pow(Mach, 4)
			+ 3.803E-19 * pow(Alpha, 4) + 7.419E-08 * pow(Mach, 5) - 3.353E-21 * pow(Alpha, 5);

		Cll_RE = +3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * RE
			+ 2.794E-08 * (Alpha * Mach) * RE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(RE, 2)
			- 4.641E-11 * pow(Alpha * Mach * RE, 2);

		Cll_LE = -(3.570E-04 - 9.569E-05 * Alpha - 3.598E-05 * Mach + 1.170E-04 * LE
			+ 2.794E-08 * (Alpha * Mach) * LE + 4.950E-06 * pow(Alpha, 2)
			+ 1.411E-06 * pow(Mach, 2) - 1.160E-06 * pow(LE, 2)
			- 4.641E-11 * pow(Alpha * Mach * LE, 2));

		Cll_RUD = -5.0103E-19 + 6.2723E-20 * Alpha + 2.3418E-20 * Mach
			+ 0.00011441 * RUD - 2.6824E-06 * (Alpha * RUD)
			- 3.4201E-21 * (Alpha * Mach) - 3.5496E-06 * (Mach * RUD)
			+ 5.5547E-08 * (Alpha * Mach) * RUD;

		Cllr = +3.82E-01 - 1.06E-01 * Mach + 1.94E-03 * Alpha
			- 8.15E-05 * (Alpha * Mach) + 1.45E-02 * pow(Mach, 2) - 9.76E-06 * pow(Alpha, 2)
			+ 4.49E-08 * pow(Alpha * Mach, 2) - 1.02E-03 * pow(Mach, 3)
			- 2.70E-07 * pow(Alpha, 3) + 3.56E-05 * pow(Mach, 4) + 3.19E-08 * pow(Alpha, 4)
			- 4.81E-07 * pow(Mach, 5) - 1.06E-09 * pow(Alpha, 5);

		Cllp = -2.99E-01 + 7.47E-02 * Mach + 1.38E-03 * Alpha
			- 8.78E-05 * (Alpha * Mach) - 9.13E-03 * pow(Mach, 2) - 2.04E-04 * pow(Alpha, 2)
			- 1.52E-07 * pow(Alpha * Mach, 2) + 5.73E-04 * pow(Mach, 3)
			- 3.86E-05 * pow(Alpha, 3) - 1.79E-05 * pow(Mach, 4) + 4.21E-06 * pow(Alpha, 4)
			+ 2.20E-07 * pow(Mach, 5) - 1.15E-07 * pow(Alpha, 5);

		Cnbv = +(0) + Alpha * (6.9980E-004) + Mach * (5.9115E-002)
			+ (Alpha * Mach) * (-7.5250E-005) + pow(Alpha, 2) * (2.5160E-004)
			+ pow(Mach, 2) * (-1.4824E-002)
			+ pow(Alpha * Mach, 2) * (-2.1924E-007)
			+ pow(Alpha, 3) * (-1.0777E-004) + pow(Mach, 3) * (1.2692E-003)
			+ pow(Alpha * Mach, 3) * (1.0707E-008)
			+ pow(Alpha, 4) * (9.4989E-006) + pow(Mach, 4) * (-4.7098E-005)
			+ pow(Alpha * Mach, 4) * (-5.5472E-011)
			+ pow(Alpha, 5) * (-2.5953E-007) + pow(Mach, 5) * (6.4284E-007)
			+ pow(Alpha * Mach, 5) * (8.5863E-014);

		Cn_LE = -(2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * LE
			- 8.93E-08 * (Alpha * Mach) * LE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(LE, 2) + 1.41E-11 * pow(Alpha * Mach * LE, 2));

		Cn_RE = 2.10E-04 + 1.83E-05 * Alpha - 3.56E-05 * Mach - 1.30E-05 * RE
			- 8.93E-08 * (Alpha * Mach) * RE - 6.39E-07 * pow(Alpha, 2) + 8.16E-07 * pow(Mach, 2)
			+ 1.97E-06 * pow(RE, 2) + 1.41E-11 * pow(Alpha * Mach * RE, 2);

		Cn_RUD = +2.85E-18 - 3.59E-19 * Alpha - 1.26E-19 * Mach - 5.28E-04 * RUD
			+ 1.39E-05 * (Alpha * RUD) + 1.57E-20 * (Alpha * Mach)
			+ 1.65E-05 * (Mach * RUD)
			- 3.13E-07 * (Alpha * Mach) * RUD;

		Cnp = +3.68E-01 - 9.79E-02 * Mach + 7.61E-16 * Alpha + 1.24E-02 * pow(Mach, 2)
			- 4.64E-16 * pow(Alpha, 2) - 8.05E-04 * pow(Mach, 3) + 1.01E-16 * pow(Alpha, 3)
			+ 2.57E-05 * pow(Mach, 4)
			- 9.18E-18 * pow(Alpha, 4) - 3.20E-07 * pow(Mach, 5) + 2.96E-19 * pow(Alpha, 5);

		Cnr = -2.41E+00 + 5.96E-01 * Mach - 2.74E-03 * Alpha
			+ 2.09E-04 * (Alpha * Mach) - 7.57E-02 * pow(Mach, 2)
			+ 1.15E-03 * pow(Alpha, 2) - 6.53E-08 * pow(Alpha * Mach, 2)
			+ 4.90E-03 * pow(Mach, 3) - 3.87E-04 * pow(Alpha, 3) - 1.57E-04 * pow(Mach, 4)
			+ 3.60E-05 * pow(Alpha, 4) + 1.96E-06 * pow(Mach, 5) - 1.18E-06 * pow(Alpha, 5);

		Cmbv = -2.192E-02 + 7.739E-03 * Mach - 2.260E-03 * Alpha
			+ 1.808E-04 * (Alpha * Mach) - 8.849E-04 * pow(Mach, 2)
			+ 2.616E-04 * pow(Alpha, 2) - 2.880E-07 * pow(Alpha * Mach, 2)
			+ 4.617E-05 * pow(Mach, 3) - 7.887E-05 * pow(Alpha, 3) - 1.143E-06 * pow(Mach, 4)
			+ 8.288E-06 * pow(Alpha, 4) + 1.082E-08 * pow(Mach, 5) - 2.789E-07 * pow(Alpha, 5);

		Cm_LE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * LE
			+ 4.48E-06 * (Alpha * LE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * LE)
			+ 9.72E-08 * (Alpha * Mach) * LE;

		Cm_RE = -5.67E-05 - 6.59E-05 * Alpha - 1.51E-06 * Mach + 2.89E-04 * RE
			+ 4.48E-06 * (Alpha * RE) - 4.46E-06 * (Alpha * Mach)
			- 5.87E-06 * (Mach * RE)
			+ 9.72E-08 * (Alpha * Mach) * RE;

		Cm_RUD = -2.79E-05 * Alpha - 5.89E-08 * pow(Alpha, 2) + 1.58E-03 * pow(Mach, 2)
			+ 6.42E-08 * pow(Alpha, 3) - 6.69E-04 * pow(Mach, 3) - 2.10E-08 * pow(Alpha, 4)
			+ 1.05E-04 * pow(Mach, 4) + 1.43E-07 * pow(RUD, 4) + 3.14E-09 * pow(Alpha, 5)
			- 7.74E-06 * pow(Mach, 5) - 4.77E-22 * pow(RUD, 5) - 2.18E-10 * pow(Alpha, 6)
			+ 2.70E-07 * pow(Mach, 6) - 3.38E-10 * pow(RUD, 6) + 5.74E-12 * pow(Alpha, 7)
			- 3.58E-09 * pow(Mach, 7) + 2.63E-24 * pow(RUD, 7);

		Cm_q = -1.36E+00 + 3.86E-01 * Mach + 7.85E-04 * Alpha
			+ 1.40E-04 * (Alpha * Mach) - 5.42E-02 * pow(Mach, 2)
			+ 2.36E-03 * pow(Alpha, 2) - 1.95E-06 * pow(Alpha * Mach, 2)
			+ 3.80E-03 * pow(Mach, 3) - 1.48E-03 * pow(Alpha, 3) - 1.30E-04 * pow(Mach, 4)
			+ 1.69E-04 * pow(Alpha, 4) + 1.71E-06 * pow(Mach, 5) - 5.93E-06 * pow(Alpha, 5);
	}

	// drag force 牵引力 x轴
	double CD = CDbv + CD_RE + CD_LE + CD_RUD;
	// lift force 升力 z轴
	double CL = CLbv + CL_RE + CL_LE;
	// side force 侧向力 y轴
	double CN = CYB * Beta + CY_RE + CY_LE + CY_RUD;
	// rolling moment 滚转力矩 绕x轴 改变滚转角速度p
	double mx = Cllbv * Beta + Cll_RE + Cll_LE + Cll_RUD + (Cllp * omega_x + Cllr * omega_z) * B / (2 * Vel);
	// pitching moment 俯仰力矩 绕y轴 改变俯仰角速度q
	double my = Cmbv + Cm_RE + Cm_LE + Cm_RUD + Cm_q * omega_y * C / (2 * Vel);
	// yawing moment 偏航力矩 绕z轴 改变偏航加速度r
	double mz = Cnbv * Beta + Cn_RE + Cn_LE + Cn_RUD + (Cnp * omega_x + Cnr * omega_z) * B / (2 * Vel);

	ans[0][0] = CDbv;		ans[0][1] = CD_LE;		ans[0][2] = CD_RE;		ans[0][3] = CD_RUD;		ans[0][4] = CD;
	ans[1][0] = CLbv;		ans[1][1] = CL_LE;		ans[1][2] = CL_RE;		ans[1][3] = CL;
	ans[2][0] = CYB;		ans[2][1] = CY_LE;		ans[2][2] = CY_RE;		ans[2][3] = CY_RUD;
	ans[3][0] = Cllbv;		ans[3][1] = Cll_LE;		ans[3][2] = Cll_RE;		ans[3][3] = Cll_RUD;	// 引起相对于机体坐标系x轴的转动，产生滚转角速度p
	ans[4][0] = Cmbv;		ans[4][1] = Cm_LE;		ans[4][2] = Cm_RE;		ans[4][3] = Cm_RUD;	// 引起相对于机体坐标系y轴的转动，产生俯仰角速度q
	ans[5][0] = Cnbv;		ans[5][1] = Cn_LE;		ans[5][2] = Cn_RE;		ans[5][3] = Cn_RUD;	// 引起相对于机体坐标系z轴的转动，产生偏航角速度r

	return;
}

void loadAeroCoefficients(const double& Mach_step, const double& Alpha_step)
{
	// 清空用于记录的txt
	FILE* File = NULL;

	// 在默认路径下(也可修改路径)，以写入方式打开txt用于记录数据，若不存在该txt则会自动生成
	File = fopen("../data/aero_coefficients.txt", "w");

	double Ma = 5, Alpha = -5;
	std::vector<std::vector<double>> aero(6, std::vector<double>(5, 0));

	// 打1500个点
	for (int i = 0; i < 301; i++) {
		Alpha = -5;
		for (int j = 0; j < 301; j++) {
			getCoefficientsDetails(Ma, Alpha, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, aero);
			// 攻角引起的阻力系数 攻角引起的升力系数 攻角引起的俯仰力矩系数 升阻比 
			fprintf(File, "%.7f %.7f %.7f %.7f %.7f %.7f ", Ma, Alpha, aero[0][0], aero[1][0], aero[4][0], aero[1][0] / aero[0][0]);
			getCoefficientsDetails(Ma, Alpha, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, aero);
			// 左舵10度引起的三轴力矩系数
			fprintf(File, "%.7f %.7f %.7f ", aero[3][1], aero[4][1], aero[5][1]);
			getCoefficientsDetails(Ma, Alpha, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, aero);
			// 方向舵10度引起的三轴力矩系数
			fprintf(File, "%.7f %.7f %.7f\n", aero[3][3], aero[4][3], aero[5][3]);
			Alpha += Alpha_step;
		}
		Ma += Mach_step;
	}
}
