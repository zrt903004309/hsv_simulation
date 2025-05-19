#include <Controller_State.h>

int sgn_s(double s);
double sat_s(double s, double k_sat);
double limit(const double& x, const double& x_pre, const double& amplitude_limit, const double& velocity_limit);

double norm(const double& value, const double& mean, const double& std);
double rnorm(const double& value, const double& mean, const double& std);

void ControllerState::Controller_State_Update(const VehicleState& Vehicle_State, const GuidanceState Guidance_State, const VehiclePara& Vehicle_Para, ModelConfig& Model_Config, torch::jit::script::Module Module, CSerialPort& mySerialPort)
{
	// 飞行器状态(直接读取)
	double X, Y, Z, V, Gamma, Chi, Alpha, Beta, Mu, p, q, r, Time;
	// 飞行器参数(直接读取)
	double Mass, Jx, Jy, Jz, B, C, S, Xcg;
	// 跟踪姿态(直接读取)
	double Alpha_Ref, Beta_Ref, Mu_Ref, dd_Alpha_Ref, dd_Beta_Ref, dd_Mu_Ref;
	// 常量获取	
	double rho, g, Ma;
	
	// 控制器参数，需要计算
	double e_Alpha, de_Alpha, dde_Alpha, e_Beta, de_Beta, dde_Beta ,e_Mu, de_Mu, dde_Mu;
	double f1, f2, f3, f4, f5, f6;

	std::vector<double> C_FM(6, 0);

	// 干扰观测器归一化

	Eigen::Vector3d F_Back = Eigen::Vector3d::Zero();
	Eigen::Vector3d v_c = Eigen::Vector3d::Zero();
	Eigen::Vector3d M_e = Eigen::Vector3d::Zero();
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
	Xcg = Vehicle_Para.Xcg;

	// 输入高度，输出声速a，马赫数Ma，重力加速度g和大气密度rho
	AtmoPara atmo;
	atmo.GetAtmoPara(-Z);
	rho = atmo.rho;
	g = atmo.g;
	Ma = V / atmo.a;

	// 期望的姿态角以及他们的二阶导参考值
	Alpha_Ref = Guidance_State.alpha_ref;
	Beta_Ref = Guidance_State.beta_ref;
	Mu_Ref = Guidance_State.mu_ref;
	dd_Alpha_Ref = Guidance_State.dd_alpha_ref;
	dd_Beta_Ref = Guidance_State.dd_beta_ref;
	dd_Mu_Ref = Guidance_State.dd_mu_ref;

	if (Time <= 0.001)
	{  // 定义e0
		e_Alpha_0 = Alpha - Alpha_Ref;
		e_Beta_0 = Beta - Beta_Ref;
		e_Mu_0 = Mu - Mu_Ref;
	}

	// 当前的姿态角的偏差（姿态角的*参考值* - 姿态角的*当前值*）
	e_Alpha = Alpha - Alpha_Ref;
	e_Beta = Beta - Beta_Ref;
	e_Mu = Mu - Mu_Ref;

	// 输出姿态角误差的一阶导（一阶微分迭代公式） 增加了一个 大的一阶滤波
	/*de_Alpha = (e_Alpha - e_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Alpha_Pre * 0.9851;
	de_Beta = (e_Beta - e_Beta_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Beta_Pre * 0.9851;
	de_Mu = (e_Mu - e_Mu_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Mu_Pre * 0.9851;*/
	de_Alpha = (e_Alpha - e_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl);
	de_Beta = (e_Beta - e_Beta_Pre) / (Model_Config.step * Model_Config.tctrl);
	de_Mu = (e_Mu - e_Mu_Pre) / (Model_Config.step * Model_Config.tctrl);

	dde_Alpha = (de_Alpha - de_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl);
	dde_Beta = (de_Beta - de_Beta_Pre) / (Model_Config.step * Model_Config.tctrl);
	dde_Mu = (de_Mu - de_Mu_Pre) / (Model_Config.step * Model_Config.tctrl);

	// 输出姿态角误差的一阶导（一阶微分迭代公式） 
	/*i_e_Alpha = e_Alpha * (Model_Config.step * Model_Config.tctrl) * 0.01 + i_e_Alpha * 0.99;
	i_e_Beta = e_Beta * (Model_Config.step * Model_Config.tctrl) * 0.01 + i_e_Beta * 0.99;
	i_e_Mu = e_Mu * (Model_Config.step * Model_Config.tctrl) * 0.01 + i_e_Mu * 0.99;*/
	i_e_Alpha += e_Alpha * (Model_Config.step * Model_Config.tctrl);
	i_e_Beta += e_Beta * (Model_Config.step * Model_Config.tctrl);
	i_e_Mu += e_Mu * (Model_Config.step * Model_Config.tctrl);
	
	//printf("三轴误差:%f %f %f 三轴误差的微分:%f %f %f 三轴误差的积分:%f %f %f\n", e_Alpha, e_Beta, e_Mu, de_Alpha, de_Beta, de_Mu, i_e_Alpha, i_e_Beta, i_e_Mu);

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
	E_Back(0, 1) = sin(Alpha) / Jx;
	E_Back(0, 2) = (cos(Alpha) / cos(Beta)) / Jx;
	E_Back(1, 0) = 1 / Jy;
	E_Back(1, 1) = 0;
	E_Back(1, 2) = 0;
	E_Back(2, 0) = -sin(Alpha) * tan(Beta) / Jz;
	E_Back(2, 1) = -cos(Alpha) / Jz;
	E_Back(2, 2) = (sin(Alpha) / cos(Beta)) / Jz;

	if (controlMode == "非线性pid全局自适应滑模容错") {
		double A_k0 = 6.25, B_k0 = 3.25, M_k0 = 5.8;
		double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		double A_c = 100, B_c = 100, Mu_c = 100;

		double A_kL = 5.0625, B_kL = 7.0625, M_kL = 7.0625;
		double A_ksi = 0.025, B_ksi = 0.01, M_ksi = 0.05;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k_pre = A_k, B_k_pre = B_k, M_k_pre = M_k;

		A_k = A_k0 + A_b * exp(-A_c * pow(e_Alpha / e_Alpha_0, 2));
		B_k = B_k0 + B_b * exp(-B_c * pow(e_Beta / e_Beta_0, 2));
		M_k = M_k0 + M_b * exp(-Mu_c * pow(e_Mu / e_Mu_0, 2));

		// 容错部分
		double gamma = 0.1;
		double d_p1 = gamma * s_Alpha * v_c[0], d_p2 = gamma * s_Beta * v_c[1], d_p3 = gamma * s_Mu * v_c[2];

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}
		else {
			s_Alpha_0 -= (A_k - A_k_pre) * e_Alpha * (Model_Config.step * Model_Config.tctrl);
			s_Beta_0 -= (B_k - B_k_pre) * e_Beta * (Model_Config.step * Model_Config.tctrl);
			s_Mu_0 -= (M_k - M_k_pre) * e_Mu * (Model_Config.step * Model_Config.tctrl);
		}

		double s_Alpha_pre = s_Alpha, s_Beta_pre = s_Beta, s_Mu_pre = s_Mu;
		// 构造时变滑模函数
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		//double dA_eta = abs(s_Alpha) / A_ksi, dB_eta = abs(s_Beta) / B_ksi, dM_eta = abs(s_Mu) / M_ksi;
		double dA_eta = (-0.005 * A_eta + abs(s_Alpha)) / A_ksi, dB_eta = (-0.005 * B_eta + abs(s_Beta)) / B_ksi, dM_eta = (-0.005 * M_eta + abs(s_Mu)) / M_ksi;

		printf("系数:%.10f %.10f %.10f \n", A_eta, B_eta, M_eta);
		printf("滑模函数值:%.10f %.10f %.10f \n", s_Alpha, s_Beta, s_Mu);

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha) + s_Alpha - s_Alpha_pre) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta) + s_Beta - s_Beta_pre) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu) + s_Mu - s_Mu_pre) + dd_Mu_Ref;

		A_eta += dA_eta * (Model_Config.step * Model_Config.tctrl), B_eta += dB_eta * (Model_Config.step * Model_Config.tctrl), M_eta += dM_eta * (Model_Config.step * Model_Config.tctrl);

		// 容错功能
		double A_gamma = 0.001, B_gamma = 0.001, M_gamma = 0.001;
		double dpA = -A_gamma * s_Alpha * v_c[0];
		double dpB = -B_gamma * s_Beta * v_c[1];
		double dpM = -M_gamma * s_Mu * v_c[2];
		pA += dpA * (Model_Config.step * Model_Config.tctrl);
		pB += dpB * (Model_Config.step * Model_Config.tctrl);
		pM += dpM * (Model_Config.step * Model_Config.tctrl);
		v_c[0] = F_Back(0) + pA * v_c[0];
		v_c[1] = F_Back(1) + pB * v_c[1];
		v_c[2] = F_Back(2) + pM * v_c[2];
	}
	else if (controlMode == "非线性pid全局自适应滑模") {
		// 1. 非线性pid全局自适应滑模
		//double A_k0 = 2.25, B_k0 = 2.25, M_k0 = 2.25;
		//double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		//double A_c = 100, B_c = 100, Mu_c = 100;

		//double A_kL = 9.0625, B_kL = 7.0625, M_kL = 7.0625;
		//// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		//double A_ksi = 0.1, B_ksi = 0.005, M_ksi = 0.01;
		//double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k0 = 6.25, B_k0 = 3.25, M_k0 = 5.8;
		double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		double A_c = 100, B_c = 100, Mu_c = 100;

		double A_kL = 5.0625, B_kL = 7.0625, M_kL = 7.0625;
		double A_ksi = 0.025, B_ksi = 0.01, M_ksi = 0.05;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k_pre = A_k, B_k_pre = B_k, M_k_pre = M_k;

		A_k = A_k0 + A_b * exp(-A_c * pow(e_Alpha / e_Alpha_0, 2));
		B_k = B_k0 + B_b * exp(-B_c * pow(e_Beta / e_Beta_0, 2));
		M_k = M_k0 + M_b * exp(-Mu_c * pow(e_Mu / e_Mu_0, 2));

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}
		else {
			s_Alpha_0 -= (A_k - A_k_pre) * e_Alpha * (Model_Config.step * Model_Config.tctrl);
			s_Beta_0 -= (B_k - B_k_pre) * e_Beta * (Model_Config.step * Model_Config.tctrl);
			s_Mu_0 -= (M_k - M_k_pre) * e_Mu * (Model_Config.step * Model_Config.tctrl);
		}

		double s_Alpha_pre = s_Alpha, s_Beta_pre = s_Beta, s_Mu_pre = s_Mu;
		// 构造时变滑模函数
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		//double dA_eta = abs(s_Alpha) / A_ksi, dB_eta = abs(s_Beta) / B_ksi, dM_eta = abs(s_Mu) / M_ksi;
		double dA_eta = (-0.005 * A_eta + abs(s_Alpha)) / A_ksi, dB_eta = (-0.005 * B_eta + abs(s_Beta)) / B_ksi, dM_eta = (-0.005 * M_eta + abs(s_Mu)) / M_ksi;


		printf("系数:%.10f %.10f %.10f \n", A_eta, B_eta, M_eta);
		printf("滑模函数值:%.10f %.10f %.10f \n", s_Alpha, s_Beta, s_Mu);

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha) + s_Alpha - s_Alpha_pre) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta) + s_Beta - s_Beta_pre) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu) + s_Mu - s_Mu_pre) + dd_Mu_Ref;

		A_eta += dA_eta * (Model_Config.step * Model_Config.tctrl), B_eta += dB_eta * (Model_Config.step * Model_Config.tctrl), M_eta += dM_eta * (Model_Config.step * Model_Config.tctrl);
	}	
	else if (controlMode == "非线性全局pid滑模") {
		// 2. 非线性全局pid滑模
		//double A_k0 = 2.25, B_k0 = 2.25, M_k0 = 2.25;
		//double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		//double A_c = 100, B_c = 100, Mu_c = 100;

		//double A_kL = 9.0625, B_kL = 7.0625, M_kL = 7.0625;
		//// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		//double A_eta = 2, B_eta = 2.5, M_eta = 1.5;

		double A_k0 = 6.25, B_k0 = 3.25, M_k0 = 4.25;
		double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		double A_c = 100, B_c = 100, Mu_c = 100;

		double A_kL = 5.0625, B_kL = 7.0625, M_kL = 7.0625;
		// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		A_eta = 7, B_eta = 2.5, M_eta = 1.5;

		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k_pre = A_k, B_k_pre = B_k, M_k_pre = M_k;

		A_k = A_k0 + A_b * exp(-A_c * pow(e_Alpha / e_Alpha_0, 2));
		B_k = B_k0 + B_b * exp(-B_c * pow(e_Beta / e_Beta_0, 2));
		M_k = M_k0 + M_b * exp(-Mu_c * pow(e_Mu / e_Mu_0, 2));

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}
		else {
			s_Alpha_0 -= (A_k - A_k_pre) * e_Alpha * (Model_Config.step * Model_Config.tctrl);
			s_Beta_0  -= (B_k - B_k_pre) * e_Beta * (Model_Config.step * Model_Config.tctrl);
			s_Mu_0    -= (M_k - M_k_pre) * e_Mu * (Model_Config.step * Model_Config.tctrl);
		}

		////构造时变滑模函数
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha)) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta)) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu)) + dd_Mu_Ref;
	}	
	else if (controlMode == "全局pid滑模"){
		//double A_k = 3.15, B_k = 3.15, M_k = 3.15;是这样的

		//A_k = 4.95, B_k = 4.95, M_k = 4.95;
		//double A_kL = 9.0625, B_kL = 7.0625, M_kL = 7.0625;
		//// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		//double A_eta = 2, B_eta = 2.5, M_eta = 1.5;
		//double A_lambda = 100, B_lambda = 100, M_lambda = 100;
		A_k = 15.5, B_k = 6, M_k = 3.95;
		double A_kL = 15.0625, B_kL = 10.0625, M_kL = 10.0625;
		// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		A_eta = 3, B_eta = 3, M_eta = 2;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}

		// 构造时变滑模函数hhh
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		//计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha)) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta)) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu)) + dd_Mu_Ref;
	}
	else if (controlMode == "基于DANN干扰观测器的全局pid滑模") {
		//double A_k = 3.15, B_k = 3.15, M_k = 3.15;是这样的

		//A_k = 4.95, B_k = 4.95, M_k = 4.95;
		//double A_kL = 9.0625, B_kL = 7.0625, M_kL = 7.0625;
		//// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		//double A_eta = 2, B_eta = 2.5, M_eta = 1.5;
		//double A_lambda = 100, B_lambda = 100, M_lambda = 100;
		A_k = 15.5, B_k = 6, M_k = 3.95;
		double A_kL = 15.0625, B_kL = 10.0625, M_kL = 10.0625;
		// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		A_eta = 3, B_eta = 3, M_eta = 2;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		/*s_Alpha_0 = s_Alpha_0;
		s_Beta_0 = s_Beta_0;
		s_Mu_0 = s_Mu_0;*/
		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}

		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);
		//s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0;
		//s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0;
		//s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0;

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha)) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta)) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu)) + dd_Mu_Ref;

		// 创建模型输入
		// Delta1, Delta2, Delta3, Alpha, Beta, Mu, V, Mcx, Mcy, Mcz
		torch::Tensor input = torch::tensor(
			{ norm(Delta[0], -1.2986101244, 1.1425814102),
			norm(Delta[1], -1.6035383472, 0.9798966840),
			norm(Delta[2], 0.0698074510, 0.1292082638),
			norm(Alpha, 6.4139936042, 3.5818639269),
			norm(Beta, -0.0042144665, 0.0262834744),
			norm(Mu, -1.0425787721, 0.1641854773),
			norm(V, 2990.1751482785, 4.3831317984),
			norm(Mc[0], -9352.7491174128, 10212.3247441523),
			norm(Mc[0], 313457.1866609525, 184172.0793187855),
			norm(Mc[0], -5805.9316447228, 13459.4385586405) },
			torch::kFloat64
		).view({ 1, 10 });  // 调整形状为[1, 10]，batch_size=1，可加.to(device)
		//std::cout << "输入张量形状: " << input.sizes() << std::endl;
		std::vector<torch::jit::IValue> inputs;
		inputs.push_back(input);

		// 前向推理
		at::Tensor output = Module.forward(inputs).toTensor();
		//std::cout << "输出张量形状: " << output.sizes() << std::endl;
		//std::cout << "输出值: " << output << std::endl;
		// 读取输出
		double* data = output.data_ptr<double>();
		phi(0) = data[0]; phi(1) = data[1]; phi(2) = data[2]; phi(3) = data[3];
		//phi(0) = 1; phi(1) = 1; phi(2) = 1; phi(3) = 1;

		///*dP = -2 * P + Q - 1/kalman_r *(P * phi * (phi.transpose() * P));
		//std::cout << "矩阵 P:\n" << P << std::endl;
		//std::cout << "矩阵 dP:\n" << dP << std::endl;
		//std::cout << "向量 phi:\n" << phi << std::endl;*/
		////std::cout << "矩阵 Q:\n" << Q << std::endl;
		//
		//mf_d(0) = dde_Alpha - v_c(0); mf_d(1) = dde_Beta - v_c(1); mf_d(2) = dde_Mu - v_c(2);
		////Eigen::Vector3d e_v;
		////e_v << de_Alpha, de_Beta, de_Mu;
		////std::cout << "向量 e_v:\n" << e_v << std::endl;
		///*pf_d = phi.transpose() * theta;
		//pf_d(0) = rnorm(pf_d(0), -0.2256511907, 1.2129823780);
		//pf_d(1) = rnorm(pf_d(1), -0.0090987499, 0.0416966731);
		//pf_d(2) = rnorm(pf_d(2), -0.2256511907, 1.2129823780);*/
		////dtheta = -lambda * theta - 1 / kalman_r * P * phi * (phi.transpose() * theta - mf_d.transpose()) + P * phi * e_v.transpose();
		////std::cout << "矩阵 dtheta:\n" << dtheta << std::endl;
		////std::cout << "矩阵 theta:\n" << theta << std::endl;
		////std::cout << "向量 pf_d:\n" << pf_d << std::endl;

		////v_c[0] -= pf_d(0); v_c[1] -= pf_d(1); v_c[2] -= pf_d(2);

		////P += dP;
		////theta += dtheta;*/

		theta << -3.89878442, 3.59229589, -2.45466025,
			0.84073118, 5.53631673, 0.89045204,
			-4.70226172, -2.51898888, 2.54010773,
			0.13321052, -0.09528556, -0.10822285;
		pf_d = phi.transpose() * theta;
		pf_d(0) = rnorm(pf_d(0), -0.2256511907, 1.2129823780);
		pf_d(1) = rnorm(pf_d(1), -0.0090987499, 0.0416966731);
		pf_d(2) = rnorm(pf_d(2), -0.2256511907, 1.2129823780);
		//dtheta = -lambda * theta - 1 / kalman_r * P * phi * (phi.transpose() * theta - mf_d.transpose()) + P * phi * e_v.transpose();
		//std::cout << "矩阵 dtheta:\n" << dtheta << std::endl;
		//std::cout << "矩阵 theta:\n" << theta << std::endl;
		//std::cout << "向量 pf_d:\n" << pf_d << std::endl;

		v_c[0] = v_c[0] - 0.02 * pf_d(0); v_c[1] = v_c[1] - 0.02 * pf_d(1); v_c[2] = v_c[2] - 0.02 * pf_d(2);
	}
	else if (controlMode == "基于DANN干扰观测器的非线性全局pid滑模") {
		double A_k0 = 6.25, B_k0 = 3.25, M_k0 = 4.25;
		double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		double A_c = 100, B_c = 100, Mu_c = 100;

		double A_kL = 5.0625, B_kL = 7.0625, M_kL = 7.0625;
		// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		A_eta = 7, B_eta = 2.5, M_eta = 1.5;

		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k_pre = A_k, B_k_pre = B_k, M_k_pre = M_k;

		A_k = A_k0 + A_b * exp(-A_c * pow(e_Alpha / e_Alpha_0, 2));
		B_k = B_k0 + B_b * exp(-B_c * pow(e_Beta / e_Beta_0, 2));
		M_k = M_k0 + M_b * exp(-Mu_c * pow(e_Mu / e_Mu_0, 2));

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}
		else {
			s_Alpha_0 -= (A_k - A_k_pre) * e_Alpha * (Model_Config.step * Model_Config.tctrl);
			s_Beta_0 -= (B_k - B_k_pre) * e_Beta * (Model_Config.step * Model_Config.tctrl);
			s_Mu_0 -= (M_k - M_k_pre) * e_Mu * (Model_Config.step * Model_Config.tctrl);
		}

		////构造时变滑模函数
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha)) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta)) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu)) + dd_Mu_Ref;

		// 创建模型输入
		// Delta1, Delta2, Delta3, Alpha, Beta, Mu, V, Mcx, Mcy, Mcz
		torch::Tensor input = torch::tensor(
			{ norm(Delta[0], -1.2986101244, 1.1425814102),
			norm(Delta[1], -1.6035383472, 0.9798966840),
			norm(Delta[2], 0.0698074510, 0.1292082638),
			norm(Alpha, 6.4139936042, 3.5818639269),
			norm(Beta, -0.0042144665, 0.0262834744),
			norm(Mu, -1.0425787721, 0.1641854773),
			norm(V, 2990.1751482785, 4.3831317984),
			norm(Mc[0], -9352.7491174128, 10212.3247441523),
			norm(Mc[0], 313457.1866609525, 184172.0793187855),
			norm(Mc[0], -5805.9316447228, 13459.4385586405) },
			torch::kFloat64
		).view({ 1, 10 });  // 调整形状为[1, 10]，batch_size=1，可加.to(device)
		//std::cout << "输入张量形状: " << input.sizes() << std::endl;
		std::vector<torch::jit::IValue> inputs;
		inputs.push_back(input);

		//// 前向推理
		at::Tensor output = Module.forward(inputs).toTensor();
		//std::cout << "输出张量形状: " << output.sizes() << std::endl;
		//std::cout << "输出值: " << output << std::endl;
		// 读取输出
		double* data = output.data_ptr<double>();
		phi(0) = data[0]; phi(1) = data[1]; phi(2) = data[2]; phi(3) = data[3];
		//phi(0) = 1; phi(1) = 1; phi(2) = 1; phi(3) = 1;

		//dP = -2 * P + Q - 1 / kalman_r * (P * phi * (phi.transpose() * P));
		//std::cout << "矩阵 P:\n" << P << std::endl;
		//std::cout << "矩阵 dP:\n" << dP << std::endl;
		//std::cout << "向量 phi:\n" << phi << std::endl;
		//std::cout << "矩阵 Q:\n" << Q << std::endl;

		mf_d(0) = dde_Alpha - v_c(0); mf_d(1) = dde_Beta - v_c(1); mf_d(2) = dde_Mu - v_c(2);
		//Eigen::Vector3d e_v;
		//e_v << de_Alpha, de_Beta, de_Mu;
		//std::cout << "向量 e_v:\n" << e_v << std::endl;

		// theta矩阵的最小二乘结果
		theta << -3.89878442, 3.59229589, -2.45466025,
			0.84073118, 5.53631673, 0.89045204,
			-4.70226172, -2.51898888, 2.54010773,
			0.13321052, -0.09528556, -0.10822285;
		pf_d = phi.transpose() * theta;
		pf_d(0) = rnorm(pf_d(0), -0.2256511907, 1.2129823780);
		pf_d(1) = rnorm(pf_d(1), -0.0090987499, 0.0416966731);
		pf_d(2) = rnorm(pf_d(2), -0.2256511907, 1.2129823780);
		//dtheta = -lambda * theta - 1 / kalman_r * P * phi * (phi.transpose() * theta - mf_d.transpose()) + P * phi * e_v.transpose();
		//std::cout << "矩阵 dtheta:\n" << dtheta << std::endl;
		//std::cout << "矩阵 theta:\n" << theta << std::endl;
		//std::cout << "向量 pf_d:\n" << pf_d << std::endl;

		v_c[0] = v_c[0] - 0.02 * pf_d(0); v_c[1] = v_c[1] - 0.02 * pf_d(1); v_c[2] = v_c[2] - 0.02 * pf_d(2);

		//P += dP;
		//theta += dtheta;
	}
	else if (controlMode == "基于DANN干扰观测器的非线性pid全局自适应滑模") {
		// 1. 非线性pid全局自适应滑模
		//double A_k0 = 2.25, B_k0 = 2.25, M_k0 = 2.25;
		//double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		//double A_c = 100, B_c = 100, Mu_c = 100;

		//double A_kL = 9.0625, B_kL = 7.0625, M_kL = 7.0625;
		//// 因为滑模函数中有积分项，为了消除稳态误差，可以增大对应通道的eta项
		//double A_ksi = 0.1, B_ksi = 0.005, M_ksi = 0.01;
		//double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k0 = 6.25, B_k0 = 3.25, M_k0 = 5.8;
		double A_b = 2.7, B_b = 2.7, M_b = 2.7;
		double A_c = 100, B_c = 100, Mu_c = 100;

		double A_kL = 5.0625, B_kL = 7.0625, M_kL = 7.0625;
		double A_ksi = 0.025, B_ksi = 0.01, M_ksi = 0.05;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		double A_k_pre = A_k, B_k_pre = B_k, M_k_pre = M_k;

		A_k = A_k0 + A_b * exp(-A_c * pow(e_Alpha / e_Alpha_0, 2));
		B_k = B_k0 + B_b * exp(-B_c * pow(e_Beta / e_Beta_0, 2));
		M_k = M_k0 + M_b * exp(-Mu_c * pow(e_Mu / e_Mu_0, 2));

		if (Time <= 0.001)
		{  // 定义s0
			s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
			s_Beta_0 = -(de_Beta + B_k * e_Beta);
			s_Mu_0 = -(de_Mu + M_k * e_Mu);
		}
		else {
			s_Alpha_0 -= (A_k - A_k_pre) * e_Alpha * (Model_Config.step * Model_Config.tctrl);
			s_Beta_0 -= (B_k - B_k_pre) * e_Beta * (Model_Config.step * Model_Config.tctrl);
			s_Mu_0 -= (M_k - M_k_pre) * e_Mu * (Model_Config.step * Model_Config.tctrl);
		}

		double s_Alpha_pre = s_Alpha, s_Beta_pre = s_Beta, s_Mu_pre = s_Mu;
		// 构造时变滑模函数
		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha + s_Alpha_0 * exp(-A_lambda * Time);
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		//double dA_eta = (-0.005 * A_eta + abs(s_Alpha)) / A_ksi, dB_eta = (-0.005 * B_eta + abs(s_Beta)) / B_ksi, dM_eta = (-0.005 * M_eta + abs(s_Mu)) / M_ksi;
		double dA_eta = abs(s_Alpha) / A_ksi, dB_eta = abs(s_Beta) / B_ksi, dM_eta = abs(s_Mu) / M_ksi;

		// 计算辅助控制率
		v_c[0] = -(A_kL * e_Alpha + A_k * de_Alpha + A_eta * tanh(s_Alpha) + s_Alpha - s_Alpha_pre) + dd_Alpha_Ref;
		v_c[1] = -(B_kL * e_Beta + B_k * de_Beta + B_eta * tanh(s_Beta) + s_Beta - s_Beta_pre) + dd_Beta_Ref;
		v_c[2] = -(M_kL * e_Mu + M_k * de_Mu + M_eta * tanh(s_Mu) + s_Mu - s_Mu_pre) + dd_Mu_Ref;

		// 创建模型输入
		// Delta1, Delta2, Delta3, Alpha, Beta, Mu, V, Mcx, Mcy, Mcz
		torch::Tensor input = torch::tensor(
			{ norm(Delta[0], -1.2986101244, 1.1425814102),
			norm(Delta[1], -1.6035383472, 0.9798966840),
			norm(Delta[2], 0.0698074510, 0.1292082638),
			norm(Alpha, 6.4139936042, 3.5818639269),
			norm(Beta, -0.0042144665, 0.0262834744),
			norm(Mu, -1.0425787721, 0.1641854773),
			norm(V, 2990.1751482785, 4.3831317984),
			norm(Mc[0], -9352.7491174128, 10212.3247441523),
			norm(Mc[0], 313457.1866609525, 184172.0793187855),
			norm(Mc[0], -5805.9316447228, 13459.4385586405) },
			torch::kFloat64
		).view({ 1, 10 });  // 调整形状为[1, 10]，batch_size=1，可加.to(device)
		//std::cout << "输入张量形状: " << input.sizes() << std::endl;
		std::vector<torch::jit::IValue> inputs;
		inputs.push_back(input);

		// 前向推理
		at::Tensor output = Module.forward(inputs).toTensor();
		//std::cout << "输出张量形状: " << output.sizes() << std::endl;
		//std::cout << "输出值: " << output << std::endl;
		// 读取输出
		double* data = output.data_ptr<double>();
		phi(0) = data[0]; phi(1) = data[1]; phi(2) = data[2]; phi(3) = data[3];
		//phi(0) = 1; phi(1) = 1; phi(2) = 1; phi(3) = 1;

		//dP = -2 * P + Q - 1 / kalman_r * (P * phi * (phi.transpose() * P));
		//std::cout << "矩阵 P:\n" << P << std::endl;
		//std::cout << "矩阵 dP:\n" << dP << std::endl;
		//std::cout << "向量 phi:\n" << phi << std::endl;
		////std::cout << "矩阵 Q:\n" << Q << std::endl;

		mf_d(0) = dde_Alpha - v_c(0); mf_d(1) = dde_Beta - v_c(1); mf_d(2) = dde_Mu - v_c(2);
		//Eigen::Vector3d e_v;
		//e_v << de_Alpha, de_Beta, de_Mu;
		//std::cout << "向量 e_v:\n" << e_v << std::endl;

		theta << -3.89878442, 3.59229589, -2.45466025,
			0.84073118, 5.53631673, 0.89045204,
			-4.70226172, -2.51898888, 2.54010773,
			0.13321052, -0.09528556, -0.10822285;
		pf_d = phi.transpose() * theta;
		pf_d(0) = rnorm(pf_d(0), -0.2256511907, 1.2129823780);
		pf_d(1) = rnorm(pf_d(1), -0.0090987499, 0.0416966731);
		pf_d(2) = rnorm(pf_d(2), -0.2256511907, 1.2129823780);
		//dtheta = -lambda * theta - 1 / kalman_r * P * phi * (phi.transpose() * theta - mf_d.transpose()) + P * phi * e_v.transpose();

		v_c[0] = v_c[0] - 0.02 * pf_d(0); v_c[1] = v_c[1] - 0.02 * pf_d(1); v_c[2] = v_c[2] - 0.02 * pf_d(2);

		A_eta += dA_eta * (Model_Config.step * Model_Config.tctrl), B_eta += dB_eta * (Model_Config.step * Model_Config.tctrl), M_eta += dM_eta * (Model_Config.step * Model_Config.tctrl);
	}

	u[0] = v_c[0], u[1] = v_c[1], u[2] = v_c[2];
	Eigen::Matrix3d E_Inv = E_Back.inverse();					// 计算E逆
	M_e = E_Inv.transpose() * (v_c - F_Back);					// 计算E逆 *（v-F）得到期望力矩
	
	/************************************************** 控制分配 **************************************************************/
	// 舵偏计算    执行机构指令限幅
	double q_dyn = 0.5 * V * V * rho;			// 计算动压
	std::vector<std::vector<double>> CoefficientDerivative(6, std::vector<double>(7, 0));
	getCoefficientsDerivative(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, CoefficientDerivative);

	// 期望的由舵产生的力矩
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	temp(0) = M_e(0) - q_dyn * B * S * CoefficientDerivative[3][6];
	temp(1) = M_e(1) - q_dyn * S * (C * CoefficientDerivative[4][5] + Xcg * (CoefficientDerivative[0][0] * sin(Alpha) + CoefficientDerivative[2][0] * cos(Alpha)));
	temp(2) = M_e(2) - q_dyn * S * (B * CoefficientDerivative[5][6] - Xcg * CoefficientDerivative[1][0] * Beta);

	//printf("期望虚拟控制:%f %f %f \n", v_c(0), v_c(1), v_c(2));
	//printf("期望力矩:%f %f %f 期望舵产生的力矩:%f %f %f \n", M_e(0), M_e(1), M_e(2), temp(0), temp(1), temp(2));

	// 控制效率矩阵(对每个舵的导数)  需要修改r那一项
	Eigen::Matrix3d gDelta = Eigen::Matrix3d::Zero();
	gDelta(0, 0) = q_dyn * B * S * CoefficientDerivative[3][1];
	gDelta(0, 1) = q_dyn * S * (C * CoefficientDerivative[4][1] + Xcg * (CoefficientDerivative[0][1] * sin(Alpha) + CoefficientDerivative[2][1] * cos(Alpha)));
	gDelta(0, 2) = q_dyn * S * (B * CoefficientDerivative[5][1] + Xcg * CoefficientDerivative[1][1]);
	gDelta(1, 0) = q_dyn * B * S * CoefficientDerivative[3][2];
	gDelta(1, 1) = q_dyn * S * (C * CoefficientDerivative[4][2] + Xcg * (CoefficientDerivative[0][2] * sin(Alpha) + CoefficientDerivative[2][2] * cos(Alpha)));
	gDelta(1, 2) = q_dyn * S * (B * CoefficientDerivative[5][2] + Xcg * CoefficientDerivative[1][2]);
	gDelta(2, 0) = q_dyn * B * S * CoefficientDerivative[3][3];
	gDelta(2, 1) = q_dyn * S * (C * CoefficientDerivative[4][3] + Xcg * CoefficientDerivative[0][3]);
	gDelta(2, 2) = q_dyn * S * (B * CoefficientDerivative[5][3] + Xcg * CoefficientDerivative[1][3]);

	Eigen::Matrix3d g_Inv = gDelta.inverse();

	Eigen::Vector3d expected_delta = g_Inv.transpose() * temp;
	//printf("期望舵角:%f %f %f \n", expected_delta(0), expected_delta(1), expected_delta(2));


	// 幅值和速度限制
	if (Time == 0) {
		Delta[0] = expected_delta(0);
		Delta[1] = expected_delta(1);
		Delta[2] = expected_delta(2);
	}
	Delta[0] = limit(expected_delta(0), Delta[0], Model_Config.delta_limit, Model_Config.d_delta_limit);
	Delta[1] = limit(expected_delta(1), Delta[1], Model_Config.delta_limit, Model_Config.d_delta_limit);
	Delta[2] = limit(expected_delta(2), Delta[2], Model_Config.delta_limit, Model_Config.d_delta_limit);
	
	if (situation == "舵面损伤" && Time > 15) {
		//Delta[0] *= 0.7;
		Delta[0] *= Model_Config.DecreaseFactor;
		Delta[1] *= Model_Config.DecreaseFactor;
		Delta[2] *= Model_Config.DecreaseFactor;
	}

	// 硬件
	if (Model_Config.hardware_en == 1) {
		if (Time <= Model_Config.step * Model_Config.tctrl) {
			mySerialPort.WriteData(0x04, (float)(Delta[0] + 180), 8);
			Sleep(5000);
			double res = (double)(mySerialPort.readFloatData() - 180);
			if (abs(res - Delta[0]) < 1) Delta[0] = res;
			std::cout << "真实舵偏角:" << Delta[0] << std::endl;
		}
		else {
			mySerialPort.WriteData(0x04, (float)(Delta[0] + 180), 8);
			Sleep(300);					// 暂停时间30ms再开始读取数据
			double res = (double)(mySerialPort.readFloatData() - 180);
			if (abs(res - Delta[0]) < 1) Delta[0] = res;
			std::cout << "真实舵偏角:" << Delta[0] << std::endl;
		}
	}
	
	//if (situation == "舵面损伤") {
	if (situation == "") {
		getCoefficients(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, C_FM, Model_Config.disturb_flag);
		M[0] = C_FM[3] * 0.5 * rho * V * V * S * B;
		M[1] = C_FM[4] * 0.5 * rho * V * V * S * B;
		M[2] = C_FM[5] * 0.5 * rho * V * V * S * B;
	}
	else {
		M[0] = M_e(0);
		M[1] = M_e(1);
		M[2] = M_e(2);
		Mc_e[0] = temp(0);
		Mc_e[1] = temp(1);
		Mc_e[2] = temp(2);
		getCoefficients(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, C_FM, Model_Config.disturb_flag);
		Mc[0] = C_FM[3] * 0.5 * rho * V * V * S * B;
		Mc[1] = C_FM[4] * 0.5 * rho * V * V * S * B;
		Mc[2] = C_FM[5] * 0.5 * rho * V * V * S * B;
	}
	
	// 保存控制器参数中的过去值
	e_Alpha_Pre=e_Alpha;
	de_Alpha_Pre=de_Alpha;
	e_Beta_Pre=e_Beta;
	de_Beta_Pre=de_Beta;
	e_Mu_Pre=e_Mu;
	de_Mu_Pre=de_Mu;
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

	e_Alpha_0 = 0.0;						
	e_Beta_0 = 0.0;							
	e_Mu_0 = 0.0;							

	M = std::vector<double>(3, 0.0);		// 期望力矩
	Mc_e = std::vector<double>(3, 0.0);		// 期望控制力矩
	Mc = std::vector<double>(3, 0.0);		// 实际控制力矩
	Delta = std::vector<double>(3, 0.0);	// 舵面角度
	u = std::vector<double>(3, 0.0);		// 虚拟控制量

	i_e_Alpha = 0.0;						// 攻角误差的积分（积分时间为控制周期）
	i_e_Beta = 0.0;							// 侧滑角误差的积分（积分时间为控制周期）
	i_e_Mu = 0.0;							// 倾侧角误差的积分（积分时间为控制周期）

	A_k = 0;
	B_k = 0;
	M_k = 0;
	A_eta = 0;
	B_eta = 0;
	M_eta = 0;
	pA = 1;
	pB = 1;
	pM = 1;

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

	// 干扰观测器相关参数设置
	mf_d = Eigen::Vector3d::Zero();
	pf_d = Eigen::Vector3d::Zero();
	phi = Eigen::Vector4d::Zero();
	theta = Eigen::MatrixXd::Zero(4, 3);
	/*theta << 1, 2, 3,
			 4, 5, 6,
			 7, 8, 9,
			 10, 11, 12;*/
	dtheta = Eigen::MatrixXd::Zero(4, 3);
	//P = Eigen::Matrix4d::Zero();
	P = 0.1 * Eigen::Matrix4d::Ones();
	dP = Eigen::Matrix4d::Ones();

	lambda = 0.01;
	kalman_q = 0.5;	// 0.1的数量级
	kalman_r = 10000;
	//Eigen::Vector4d diag_elements(kalman_q, kalman_q, kalman_q, kalman_q);
	//Q = diag_elements.asDiagonal();
	Q << 0.5, 0, 0, 0,
		0, 0.5, 0, 0,
		0, 0, 0.5, 0,
		0, 0, 0, 0.5;
	//Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();  // 对角元素为q, q, q, q
	//Q.diagonal() << kalman_q, kalman_q, kalman_q, kalman_q;
	//Q.diagonal() << 0.5, 0.5, 0.5, 0.5;
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

double sat_s(double s, double k_sat)  //饱和函数
{
	double d;
	d = s / k_sat;
	if(fabs(d) <= 1.0)
		return d;	
	else		
		return	sgn_s(s);
		//return	tanh(s);
	
}

double limit(const double& x, const double& x_pre, const double& amplitude_limit, const double& velocity_limit) {
	double velocity = (x - x_pre);		// 一个时间步的变化量
	if (velocity >= -velocity_limit && velocity <= velocity_limit) {
		// 检查速度限制，符合
		if (-amplitude_limit <= x && x <= amplitude_limit) {
			return x;
		}
		else if (x < -amplitude_limit) {
			return -amplitude_limit;
		}
		else return amplitude_limit;
	}
	else if (velocity < -velocity_limit) {
		// 调整 x 以满足速度下限
		double new_x = x_pre - velocity_limit;
		// 再检查幅值限制
		if (-amplitude_limit <= new_x && new_x <= amplitude_limit) {
			return new_x;
		}
		else if (new_x < -amplitude_limit) {
			return -amplitude_limit;
		}
		else return amplitude_limit;
	}
	else {
		// 调整 x 以满足速度上限
		double new_x = x_pre + velocity_limit;
		// 再检查幅值限制
		if (-amplitude_limit <= new_x && new_x <= amplitude_limit) {
			return new_x;
		}
		else if (new_x < -amplitude_limit) {
			return -amplitude_limit;
		}
		else return amplitude_limit;
	}
}

double norm(const double& value, const double& mean, const double& std) {
	double res = (value - mean) / std;
	return res;
}

double rnorm(const double& value, const double& mean, const double& std) {
	double res = value * std + mean;
	return res;
}
