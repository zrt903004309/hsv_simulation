#include "Controller_State.h"
#include <Vehicle_State.h>

int sgn_s(double s);
double sat_s(double s, double k_sat);
double limit(double x, double Min, double Max);

void ControllerState::Controller_State_Update(const VehicleState& Vehicle_State, const VehiclePara& Vehicle_Para, const ModelConfig& Model_Config)
{
	
	double sat_limit_a = 0.01,sat_limit_b = 0.01,sat_limit_m = 0.01;
	double X,Y,Z,V,Gamma,Chi,Alpha,Beta,Mu,p,q,r,Time;

	double Mass, Jx, Jy, Jz, B, C, S;
	
	double Alpha_Ref;
	double Beta_Ref;
	double Mu_Ref;
	double dd_Alpha_Ref;
	double dd_Beta_Ref;
	double dd_Mu_Ref;
	
	double e_Alpha,de_Alpha;
	double e_Beta,de_Beta;
	double e_Mu,de_Mu;

	double f1,f2,f3,f4,f5,f6;
	double rho, g, Ma;

	Eigen::Vector3d F_Back = Eigen::Vector3d::Zero();
	Eigen::Vector3d v_c = Eigen::Vector3d::Zero();
	Eigen::Vector3d M_control = Eigen::Vector3d::Zero();
	Eigen::Matrix3d E_Back = Eigen::Matrix3d::Zero();


	// 读飞行器状态
	Time = Vehicle_State.Time;
	V = Vehicle_State.V;
	Alpha = Vehicle_State.Alpha;
	Beta = Vehicle_State.Beta;
	Gamma = Vehicle_State.Gamma;
	Mu = Vehicle_State.Mu;
	Chi = Vehicle_State.Chi;
	p = Vehicle_State.p;
	q = Vehicle_State.q;
	r = Vehicle_State.r;
	X = Vehicle_State.X;
	Y = Vehicle_State.Y;
	Z = Vehicle_State.Z; 

	// 读飞行器参数
	Mass = Vehicle_Para.Mass;
	Jx = Vehicle_Para.Jx;
	Jy = Vehicle_Para.Jy;
	Jz = Vehicle_Para.Jz;
	B = Vehicle_Para.B;
	C = Vehicle_Para.C;
	S = Vehicle_Para.S;

	// 期望的姿态角以及他们的二阶导参考值
	Alpha_Ref = Model_Config.alpha_ref;
	Beta_Ref = Model_Config.beta_ref;
	Mu_Ref = Model_Config.mu_ref;
	dd_Alpha_Ref = Model_Config.dd_alpha_ref;
	dd_Beta_Ref = Model_Config.dd_beta_ref;
	dd_Mu_Ref = Model_Config.dd_mu_ref;

	// 当前的姿态角的偏差（姿态角的*参考值* - 姿态角的*当前值*）
	e_Alpha = Alpha_Ref-Alpha;
	e_Beta = Beta_Ref-Beta;
	e_Mu = Mu_Ref-Mu;

	// 输出姿态角误差的一阶导（一阶微分迭代公式） 增加了一个 大的一阶滤波
	de_Alpha = (e_Alpha - e_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Alpha_Pre * 0.9851;
	de_Beta = (e_Beta - e_Beta_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Beta_Pre * 0.9851;
	de_Mu = (e_Mu - e_Mu_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Mu_Pre * 0.9851;
	
	// 输入高度，输出声速a，马赫数Ma，重力加速度g和大气密度rho
	AtmoPara atmo;
	atmo.GetAtmoPara(Z);
	rho = atmo.rho;
	g = atmo.g;
	Ma = V / atmo.a;

	// 计算状态空间方程中的f(x)
	f1 = (Jy - Jz) * q * r / Jx;
	f2 = (Jz - Jx) * p * r / Jy;
	f3 = (Jx - Jy) * p * q / Jz;
	f4 = -p * cos(Alpha) * tan(Beta) + q - r * sin(Alpha) * tan(Beta);
	f5 = p * sin(Alpha) - r * cos(Alpha);
	f6 = (p * cos(Alpha) + r * sin(Alpha)) / cos(Beta);

	// 对应二阶导数中的F
	F_Back(0) = -f1 * cos(Alpha) * tan(Beta) + f2 - f3 * sin(Alpha) * tan(Beta) + f4 * (p * sin(Alpha) * tan(Beta) - r * cos(Alpha) * tan(Beta)) - f5 * (p * cos(Alpha) + r * sin(Alpha)) / (cos(Beta) * cos(Beta));
	F_Back(1) = f1 * sin(Alpha) - f3 * cos(Alpha) + f4 * (p * cos(Alpha) + r * sin(Alpha));
	F_Back(2) = f1 * cos(Alpha) / cos(Beta) + f3 * sin(Alpha) / cos(Beta) + f4 * (-p * sin(Alpha) / cos(Beta) + r * cos(Alpha) / cos(Beta)) + f5 * (p * cos(Alpha) * tan(Beta) / cos(Beta) + r * sin(Alpha) * tan(Beta) / cos(Beta));

	// 对应二阶导数中的E
	E_Back(0, 0) = -cos(Alpha) * tan(Beta) / Jx;
	E_Back(0, 1) = 1 / Jy;
	E_Back(0, 2) = -sin(Alpha) * tan(Beta) / Jz;
	E_Back(1, 0) = sin(Alpha) / Jx;
	E_Back(1, 1) = 0;
	E_Back(1, 2) = -cos(Alpha) / Jz;
	E_Back(2, 0) = (cos(Alpha) / cos(Beta)) / Jx;
	E_Back(2, 1) = 0;
	E_Back(2, 2) = (sin(Alpha) / cos(Beta)) / Jz;

	/*
	s_Alpha_0 = Controller_State.s_Alpha_0;
	s_Beta_0 = Controller_State.s_Beta_0;
	s_Mu_0 = Controller_State.s_Mu_0;
	if(Time <= 0.001)
    {  //定义s0=0
	   s_Alpha_0 = -(de_Alpha+A_k*e_Alpha);
	   s_Beta_0 = -(de_Beta+A_k*e_Beta);
	   s_Mu_0 = -(de_Mu+A_k*e_Mu);
	}

	//构造时变滑模函数
	s_Alpha = A_k*e_Alpha + de_Alpha + s_Alpha_0*exp(-A_lambda*Time);
	s_Beta = B_k*e_Beta + de_Beta + s_Beta_0*exp(-B_lambda*Time);
	s_Mu = M_k*e_Mu + de_Mu + s_Mu_0*exp(-M_lambda*Time);

	//计算辅助控制率
	v_c[0] = A_k*de_Alpha + A_eta*sat_s(s_Alpha,sat_limit_a) + dd_Alpha_Ref - s_Alpha_0*A_lambda*exp(-A_lambda*Time);
	v_c[1] = B_k*de_Beta + B_eta*sat_s(s_Beta,sat_limit_b) + dd_Beta_Ref - s_Beta_0*B_lambda*exp(-B_lambda*Time);
	v_c[2] = M_k*de_Mu + M_eta*sat_s(s_Mu,sat_limit_m) + dd_Mu_Ref - s_Mu_0*M_lambda*exp(-M_lambda*Time);
	*/
	
	
	if (controlMode == "快速滑模") {
		// 1. 快速光滑二阶滑模
		double A_k1 = 0.5, A_k2 = 0.2, A_k3 = 0.1, A_k4 = 0.01, A_k5 = 0.01, A_k6 = 0.01, A_eta = 0.001;
		double B_k1 = 0.5, B_k2 = 0.2, B_k3 = 0.1, B_k4 = 0.01, B_k5 = 0.01, B_k6 = 0.01, B_eta = 0.001;
		double M_k1 = 0.5, M_k2 = 0.2, M_k3 = 0.1, M_k4 = 0.01, M_k5 = 0.01, M_k6 = 0.01, M_eta = 0.001;
		double p1 = 5;

		i_e_Alpha = i_e_Alpha * 0.995 + e_Alpha * 0.005;
		i_e_Beta = i_e_Beta * 0.995 + e_Beta * 0.005;
		i_e_Mu = i_e_Mu * 0.995 + e_Mu * 0.005;

		s_Alpha = de_Alpha + A_k1 * e_Alpha + A_k2 * i_e_Alpha;
		s_Beta = de_Beta + B_k1 * e_Beta + B_k2 * i_e_Beta;
		s_Mu = de_Mu + M_k1 * e_Mu + M_k2 * i_e_Mu;

		i_s5_Alpha = i_s5_Alpha * 0.995 + A_k5 * pow(fabs(s_Alpha), (p1 - 2) / p1) * sat_s(s_Alpha, sat_limit_a) * 0.005;
		i_s5_Beta = i_s5_Beta * 0.995 + B_k5 * pow(fabs(s_Beta), (p1 - 2) / p1) * sat_s(s_Beta, sat_limit_b) * 0.005;
		i_s5_Mu = i_s5_Mu * 0.995 + M_k5 * pow(fabs(s_Mu), (p1 - 2) / p1) * sat_s(s_Mu, sat_limit_m) * 0.005;
		i_s6_Alpha = i_s6_Alpha * 0.995 + A_k6 * s_Alpha * 0.005;
		i_s6_Beta = i_s6_Beta * 0.995 + B_k6 * s_Beta * 0.005;
		i_s6_Mu = i_s6_Mu * 0.995 + M_k6 * s_Mu * 0.005;

		v_c(0) = dd_Alpha_Ref - A_k1 * de_Alpha - A_k2 * e_Alpha - A_k3 * pow(fabs(s_Alpha), (p1 - 1) / p1) * sat_s(s_Alpha, sat_limit_a) - A_k4 * s_Alpha - i_s5_Alpha - i_s6_Alpha - A_eta * sat_s(s_Alpha, sat_limit_a);
		v_c(1) = dd_Beta_Ref - B_k1 * de_Beta - B_k2 * e_Beta - B_k3 * pow(fabs(s_Beta), (p1 - 1) / p1) * sat_s(s_Beta, sat_limit_b) - B_k4 * s_Beta - i_s5_Beta - i_s6_Beta - B_eta * sat_s(s_Beta, sat_limit_b);
		v_c(2) = dd_Mu_Ref - M_k1 * de_Mu - M_k2 * e_Mu - M_k3 * pow(fabs(s_Mu), (p1 - 1) / p1) * sat_s(s_Mu, sat_limit_m) - M_k4 * s_Mu - i_s5_Mu - i_s6_Mu - M_eta * sat_s(s_Mu, sat_limit_m);
	}	
	else if (controlMode == "自适应滑模") {
		// 2. 自适应终端滑模
		double A_k1 = 0.01, A_k2 = 0.01, A_k3 = 0.005, A_k4 = 0.005, A_k5 = 0.1;
		double B_k1 = 0.01, B_k2 = 0.01, B_k3 = 0.005, B_k4 = 0.005, B_k5 = 0.1;
		double M_k1 = 0.01, M_k2 = 0.01, M_k3 = 0.005, M_k4 = 0.005, M_k5 = 0.1;
		double a0 = 5, p1 = 0.9, epson = 0.9;
		double a1 = 1 + 1 / a0, a2 = 1 - 1 / a0;

		i_s2_Alpha = i_s2_Alpha * 0.995 + A_k2 * pow(fabs(de_Alpha), a1) * sat_s(de_Alpha, sat_limit_a) * 0.005;
		i_s2_Beta = i_s2_Beta * 0.995 + B_k2 * pow(fabs(de_Beta), a1) * sat_s(de_Beta, sat_limit_b) * 0.005;
		i_s2_Mu = i_s2_Mu * 0.995 + M_k2 * pow(fabs(de_Mu), a1) * sat_s(de_Mu, sat_limit_m) * 0.005;
		i_s3_Alpha = i_s3_Alpha * 0.995 + A_k3 * pow(fabs(de_Alpha), a2) * sat_s(de_Alpha, sat_limit_a) * 0.005;
		i_s3_Beta = i_s3_Beta * 0.995 + B_k3 * pow(fabs(de_Beta), a2) * sat_s(de_Beta, sat_limit_b) * 0.005;
		i_s3_Mu = i_s3_Mu * 0.995 + M_k3 * pow(fabs(de_Mu), a2) * sat_s(de_Mu, sat_limit_m) * 0.005;

		s_Alpha = de_Alpha + A_k1 * e_Alpha + i_s2_Alpha + i_s3_Alpha;
		s_Beta = de_Beta + B_k1 * e_Beta + i_s2_Beta + i_s3_Beta;
		s_Mu = de_Mu + M_k1 * e_Mu + i_s2_Mu + i_s3_Mu;

		v_c(0) = dd_Alpha_Ref - A_k1 * de_Alpha - A_k2 * pow(fabs(de_Alpha), a1) * sat_s(de_Alpha, sat_limit_a) - A_k3 * pow(fabs(de_Alpha), a2) * sat_s(de_Alpha, sat_limit_a) - A_k4 * s_Alpha - A_k5 * pow(fabs(s_Alpha), p1) * sat_s(s_Alpha, sat_limit_a) - epson * s_Alpha;
		v_c(1) = dd_Beta_Ref - B_k1 * de_Beta - B_k2 * pow(fabs(de_Beta), a1) * sat_s(de_Beta, sat_limit_b) - B_k3 * pow(fabs(de_Beta), a2) * sat_s(de_Beta, sat_limit_b) - B_k4 * s_Beta - B_k5 * pow(fabs(s_Beta), p1) * sat_s(s_Beta, sat_limit_b) - epson * s_Beta;
		v_c(2) = dd_Mu_Ref - M_k1 * de_Mu - M_k2 * pow(fabs(de_Mu), a1) * sat_s(de_Mu, sat_limit_m) - M_k3 * pow(fabs(de_Mu), a2) * sat_s(de_Mu, sat_limit_m) - M_k4 * s_Mu - M_k5 * pow(fabs(s_Mu), p1) * sat_s(s_Mu, sat_limit_m) - epson * s_Mu;
	}	
	else {
		v_c(0) = 0;
		v_c(1) = 0;
		v_c(2) = 0;
	}

	v_c *= -1;
	
	// 计算行列式
	double determinant = E_Back.determinant();
	//std::cout << "Determinant: " << determinant << std::endl;

	// 检查矩阵是否可逆
	if (determinant == 0.0) {
		std::cerr << "Error: The matrix is singular and cannot be inverted." << std::endl;
	}

	Eigen::Matrix3d E_Inv = E_Back.inverse();	// 计算E逆
	M_control = E_Inv * (v_c - F_Back);				// 计算E逆*（v-F）得到控制力矩

	
	/************************************************** 控制分配 **************************************************************/
	// 舵偏计算    执行机构指令限幅
	double q_dyn = 0.5 * V * V * rho;			// 计算动压
	Eigen::Vector3d expected_dd_x_y_z = -v_c;
	Eigen::Vector3d inertia(Jx, Jy, Jz);
	Eigen::Array3d expected_Moment = expected_dd_x_y_z.array() * inertia.array();
	
	
	double delta = 0.01;
	vector<double> C_mxyz0(6);
	CoefficientsSixDoF_FandM(Ma, Alpha, Beta, 0, 0, 0, C, B, V, 0, 0, 0, C_mxyz0);
	vector<double> B_C_myz_da(6);
	CoefficientsSixDoF_FandM(Ma, Alpha, Beta, 0, delta, 0, C, B, V, 0, 0, 0, B_C_myz_da);
	vector<double> B_C_myz_dr(6);
	CoefficientsSixDoF_FandM(Ma, Alpha, Beta, 0, 0, delta, C, B, V, 0, 0, 0, B_C_myz_dr);
	vector<double> B_C_myz_de(6);
	CoefficientsSixDoF_FandM(Ma, Alpha, Beta, delta, 0, 0, C, B, V, 0, 0, 0, B_C_myz_de);

	Eigen::Map<Eigen::VectorXd> C_mxyz0_vec(C_mxyz0.data(), C_mxyz0.size());
	Eigen::Map<Eigen::VectorXd> B_C_myz_da_vec(B_C_myz_da.data(), B_C_myz_da.size());
	Eigen::Map<Eigen::VectorXd> B_C_myz_dr_vec(B_C_myz_dr.data(), B_C_myz_dr.size());
	Eigen::Map<Eigen::VectorXd> B_C_myz_de_vec(B_C_myz_de.data(), B_C_myz_de.size());

	B_C_myz_da_vec -= C_mxyz0_vec;
	B_C_myz_dr_vec -= C_mxyz0_vec;
	B_C_myz_de_vec -= C_mxyz0_vec;

	// 求在当前情况下气动系数关于delta_a, delta_r, delta_e的导数
	// 第一行是C_mx^(delta_a) C_my^(delta_a) C_mz^(delta_a)
	// 第二行是C_mx^(delta_r) C_my^(delta_r) C_mz^(delta_r)
	// 第三行是C_mx^(delta_e) C_my^(delta_e) C_mz^(delta_e)
	vector<vector<double>> B_C_myz(3, vector<double>(3));
	for (int i = 0; i < 3; ++i) {
		B_C_myz[i][0] = B_C_myz_da_vec(i + 3) / delta;
		B_C_myz[i][1] = B_C_myz_dr_vec(i + 3) / delta;
		B_C_myz[i][2] = B_C_myz_de_vec(i + 3) / delta;
	}
	vector<vector<double>> D_B_C_myz(3, vector<double>(3));
	for (int i = 0; i < 3; ++i) {
		D_B_C_myz[i][i] = B_C_myz[i][i];
	}
	vector<vector<double>> Inv_D_B_C_myz(D_B_C_myz);
	for (int i = 0; i < 3; ++i) {
		Inv_D_B_C_myz[i][i] = 1 / D_B_C_myz[i][i];
	}
	for (int i = 0; i < 3; ++i) {
		expected_Moment[i] /= q_dyn * S * B;
		expected_Moment[i] -= C_mxyz0[i + 3];
	}
	vector<double> expected_Delta(3);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			expected_Delta[i] += Inv_D_B_C_myz[i][j] * expected_Moment[j];
		}
	}

	double delta_aileron = expected_Delta[0];
	double delta_rudder = expected_Delta[1];
	double delta_elevator = expected_Delta[2];

	Delta[0] = limit(delta_elevator, -Model_Config.delta_e_limit, Model_Config.delta_e_limit);
	Delta[1] = limit(delta_aileron, -Model_Config.delta_a_limit, Model_Config.delta_a_limit);
	Delta[2] = limit(delta_rudder, -Model_Config.delta_r_limit, Model_Config.delta_r_limit);

	// 这是舵面角的幅值限制，还需要添加舵面角的角速度限制不能超过 200deg/s

	if (situation == " ") {
		Delta[0] = 20 / rad + Delta[0];
		Delta[1] = 10 * Delta[1];
		Delta[2] = -10 /rad + 3 * Delta[2];
	}
	else {
		Delta[0] = -20 / rad + Delta[0];
		Delta[1] = -15 / rad + 10 * Delta[1];
		Delta[2] = 3 * Delta[2];
	}
	
	
	if (situation == "舵面损伤") {
		if (Time > 10 && Time < 10.2) {
			Delta[0] += 1.0 * (rand() % 100) / 15 / rad;
			Delta[1] += 1.0 * (rand() % 100) / 15 / rad;
			Delta[2] += 1.0 * (rand() % 100) / 15 / rad;
		}
	}
	else if (situation == "舵面卡死") {
		if (Time > 10 && Time < 10.2) {
			Delta[1] += 1.0 * (rand() % 100) / 15 / rad;
			Delta[2] += 1.0 * (rand() % 100) / 15 / rad;
		}
		if (Time > 10) {
			Delta[1] += 6 / rad;
			Delta[2] += -16 / rad;
		}
		if (Time > 10.2) Delta[2] += 2 / rad;
	}
	
	if (situation == " ") {
		vector<double> C_FM(6);
		CoefficientsSixDoF_FandM(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, C_FM);
		M_c[0] = C_FM[3] * 0.5 * rho * V * V * S * B;
		M_c[1] = C_FM[4] * 0.5 * rho * V * V * S * B;
		M_c[2] = C_FM[5] * 0.5 * rho * V * V * S * B;
	}
	else {
		M_c[0] = M_control(0);
		M_c[1] = M_control(1);
		M_c[2] = M_control(2);
	}
	
	// 保存控制器参数中的过去值
	e_Alpha_Pre=e_Alpha;
	de_Alpha_Pre=de_Alpha;
	e_Beta_Pre=e_Beta;
	de_Beta_Pre=de_Beta;
	e_Mu_Pre=e_Mu;
	de_Mu_Pre=de_Mu;


	if (Time > 10) {
		if (situation == "舵面损伤") {
			Delta[0] = 0.6 * Delta[0];
			Delta[1] = 0.5 * Delta[1];
			Delta[2] = 0.4 * Delta[2];
		}
		else if (situation == "舵面卡死") {
			Delta[0] = 10 / rad;
			Delta[1] = Delta[1];
			Delta[2] = Delta[2];
		}
	}
}

ControllerState::ControllerState()  //控制器状态初始化
{

	e_Alpha_Pre = 0.0;						// 攻角的偏差（参考值-实际值）
	de_Alpha_Pre = 0.0;						// 攻角的偏差的导数

	e_Beta_Pre = 0.0;						// 侧滑角的偏差（参考值-实际值）
	de_Beta_Pre = 0.0;						// 侧滑角的偏差的导数

	e_Mu_Pre = 0.0;							// 倾侧角的偏差（参考值-实际值）
	de_Mu_Pre = 0.0;						// 倾侧角的偏差的导数

	s_Alpha = 0.0;							// 攻角的滑膜函数
	s_Beta = 0.0;							// 侧滑角的滑膜函数
	s_Mu = 0.0;								// 倾侧角的滑膜函数

	s_Alpha_0 = 0.0;						// 攻角的滑膜函数初值
	s_Beta_0 = 0.0;							// 侧滑角的滑膜函数初值
	s_Mu_0 = 0.0;							// 倾侧角的滑膜函数初值

	M_c = std::vector<double>(3, 0.0);		// 控制力矩
	Delta = std::vector<double>(3, 0.0);	// 舵面角度

	i_e_Alpha = 0.0;						// 攻角误差的积分（积分时间为控制周期）
	i_e_Beta = 0.0;							// 侧滑角误差的积分（积分时间为控制周期）
	i_e_Mu = 0.0;							// 倾侧角误差的积分（积分时间为控制周期）

	i_s5_Alpha = 0.0;
	i_s5_Beta = 0.0;
	i_s5_Mu = 0.0;
	i_s6_Alpha = 0.0;
	i_s6_Beta = 0.0;
	i_s6_Mu = 0.0;

	i_s2_Alpha = 0.0;
	i_s2_Beta = 0.0;
	i_s2_Mu = 0.0;
	i_s3_Alpha = 0.0;
	i_s3_Beta = 0.0;
	i_s3_Mu = 0.0;

}

int sgn_s(double s) //符号函数
{
	if(s > 0)
		return 1;
	else if(s < 0)
		return -1;
	else
		return 0;
}

double sat_s(double s,double k_sat)  //饱和函数
{
	double d;
	d = s / k_sat;
	if(fabs(d) <= 1.0)
		return d;	
	else		
		return	sgn_s(s);
	
}

double limit(double x, double Min, double Max) {
	if (Min <= x && x <= Max) return x;
	else if (x < Min) return Min;
	else return Max;
}
