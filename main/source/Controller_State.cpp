#include <Controller_State.h>

int sgn_s(double s);
double sat_s(double s, double k_sat);
double limit(const double& x, const double& x_pre, const double& amplitude_limit, const double& velocity_limit);

void ControllerState::Controller_State_Update(const VehicleState& Vehicle_State, const GuidanceState Guidance_State, const VehiclePara& Vehicle_Para, const ModelConfig& Model_Config)
{
	
	double sat_limit_a = 0.02,sat_limit_b = 0.02,sat_limit_m = 0.02;
	// 飞行器状态(直接读取)
	double X, Y, Z, V, Gamma, Chi, Alpha, Beta, Mu, p, q, r, Time;
	// 飞行器参数(直接读取)
	double Mass, Jx, Jy, Jz, B, C, S, Xcg;
	// 跟踪姿态(直接读取)
	double Alpha_Ref, Beta_Ref, Mu_Ref, dd_Alpha_Ref, dd_Beta_Ref, dd_Mu_Ref;
	// 常量获取	
	double rho, g, Ma;
	
	// 控制器参数，需要计算
	double e_Alpha, de_Alpha, e_Beta, de_Beta, e_Mu, de_Mu;
	double f1, f2, f3, f4, f5, f6;

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

	// 当前的姿态角的偏差（姿态角的*参考值* - 姿态角的*当前值*）
	e_Alpha = Alpha_Ref - Alpha;
	e_Beta = Beta_Ref - Beta;
	e_Mu = Mu_Ref - Mu;

	// 输出姿态角误差的一阶导（一阶微分迭代公式） 增加了一个 大的一阶滤波
	/*de_Alpha = (e_Alpha - e_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Alpha_Pre * 0.9851;
	de_Beta = (e_Beta - e_Beta_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Beta_Pre * 0.9851;
	de_Mu = (e_Mu - e_Mu_Pre) / (Model_Config.step * Model_Config.tctrl) + de_Mu_Pre * 0.9851;*/
	de_Alpha = (e_Alpha - e_Alpha_Pre) / (Model_Config.step * Model_Config.tctrl);
	de_Beta = (e_Beta - e_Beta_Pre) / (Model_Config.step * Model_Config.tctrl);
	de_Mu = (e_Mu - e_Mu_Pre) / (Model_Config.step * Model_Config.tctrl);

	// 输出姿态角误差的一阶导（一阶微分迭代公式） 增加了一个 大的一阶滤波
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

		v_c[0] = dd_Alpha_Ref - A_k1 * de_Alpha - A_k2 * e_Alpha - A_k3 * pow(fabs(s_Alpha), (p1 - 1) / p1) * sat_s(s_Alpha, sat_limit_a) - A_k4 * s_Alpha - i_s5_Alpha - i_s6_Alpha - A_eta * sat_s(s_Alpha, sat_limit_a);
		v_c[1] = dd_Beta_Ref - B_k1 * de_Beta - B_k2 * e_Beta - B_k3 * pow(fabs(s_Beta), (p1 - 1) / p1) * sat_s(s_Beta, sat_limit_b) - B_k4 * s_Beta - i_s5_Beta - i_s6_Beta - B_eta * sat_s(s_Beta, sat_limit_b);
		v_c[2] = dd_Mu_Ref - M_k1 * de_Mu - M_k2 * e_Mu - M_k3 * pow(fabs(s_Mu), (p1 - 1) / p1) * sat_s(s_Mu, sat_limit_m) - M_k4 * s_Mu - i_s5_Mu - i_s6_Mu - M_eta * sat_s(s_Mu, sat_limit_m);
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
	else if (controlMode == "全局pid滑模"){
		double A_k = 3.15, B_k = 3.15, M_k = 3.15;
		double A_kL = 5.0625, B_kL = 5.0625, M_kL = 5.0625;
		double A_eta = 0.2, B_eta = 0.05, M_eta = 1;
		double A_lambda = 100, B_lambda = 100, M_lambda = 100;

		//s_Alpha_0 = s_Alpha_0;
		//s_Beta_0 = s_Beta_0;
		//s_Mu_0 = s_Mu_0;
		//if (Time <= 0.001)
		//{  //定义s0=0
		//	s_Alpha_0 = -(de_Alpha + A_k * e_Alpha);
		//	s_Beta_0 = -(de_Beta + B_k * e_Beta);
		//	s_Mu_0 = -(de_Mu + M_k * e_Mu);
		//}

		////构造时变滑模函数
		//s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha  + s_Alpha_0 * exp(-A_lambda * Time);
		//s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta + s_Beta_0 * exp(-B_lambda * Time);
		//s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu + s_Mu_0 * exp(-M_lambda * Time);

		s_Alpha = A_k * e_Alpha + de_Alpha + A_kL * i_e_Alpha;
		s_Beta = B_k * e_Beta + de_Beta + B_kL * i_e_Beta;
		s_Mu = M_k * e_Mu + de_Mu + M_kL * i_e_Mu;

		//计算辅助控制率
		v_c[0] = A_kL * e_Alpha + A_k * de_Alpha + A_eta * sat_s(s_Alpha, sat_limit_a) + dd_Alpha_Ref;
		v_c[1] = B_kL * e_Beta + B_k * de_Beta + B_eta * sat_s(s_Beta, sat_limit_b) + dd_Beta_Ref;
		v_c[2] = M_kL * e_Mu + M_k * de_Mu + M_eta * sat_s(s_Mu, sat_limit_m) + dd_Mu_Ref;
	}
	
	// 计算行列式
	double determinant = E_Back.determinant();

	Eigen::Matrix3d E_Inv = E_Back.inverse();					// 计算E逆
	M_control = E_Inv.transpose() * (v_c - F_Back);				// 计算E逆 *（v-F）得到控制力矩
	
	/************************************************** 控制分配 **************************************************************/
	// 舵偏计算    执行机构指令限幅
	double q_dyn = 0.5 * V * V * rho;			// 计算动压
	vector<vector<double>> CoefficientDerivative(6, vector<double>(7, 0));
	getCoefficientsDerivative(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, CoefficientDerivative);

	// 期望的由舵产生的力矩
	Eigen::Vector3d temp = Eigen::Vector3d::Zero();
	temp(0) = M_control(0) - q_dyn * B * S * CoefficientDerivative[3][6];
	temp(1) = M_control(1) - q_dyn * S * (C * CoefficientDerivative[4][5] + Xcg * (CoefficientDerivative[0][0] * sin(Alpha) + CoefficientDerivative[2][0] * cos(Alpha)));
	temp(2) = M_control(2) - q_dyn * S * (B * CoefficientDerivative[5][6] - Xcg * CoefficientDerivative[1][0] * Beta);

	printf("期望虚拟控制:%f %f %f \n", v_c(0), v_c(1), v_c(2));
	printf("期望控制力矩:%f %f %f 期望舵产生的力矩:%f %f %f \n", M_control(0), M_control(1), M_control(2), temp(0), temp(1), temp(2));

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
	printf("期望舵角:%f %f %f \n", expected_delta(0), expected_delta(1), expected_delta(2));

	/*if (Vehicle_State.Time > 20) {
		expected_delta[0] *= 0.9;
	}*/

	// 幅值和速度限制
	Delta[0] = limit(expected_delta(0), Delta[0], Model_Config.delta_limit, Model_Config.d_delta_limit);
	Delta[1] = limit(expected_delta(1), Delta[1], Model_Config.delta_limit, Model_Config.d_delta_limit);
	Delta[2] = limit(expected_delta(2), Delta[2], Model_Config.delta_limit, Model_Config.d_delta_limit);
	
	/*if (Time > 10) {
		Delta[0] *= 0.7;
			Delta[2] = 0.4 * Delta[2];
		else if (situation == "舵面卡死") {
			Delta[0] = Delta[0] + Model_Config.BiasFactor / rad;
			Delta[1] = Delta[1];
			Delta[2] = Delta[2];
		}
		std::cout << v_c.transpose() << endl;
	}*/
	//if (situation == "舵面损伤") {
	if (situation == "") {
		vector<double> C_FM(6, 0);
		getCoefficients(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, C_FM);
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
