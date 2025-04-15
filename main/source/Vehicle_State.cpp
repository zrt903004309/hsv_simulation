#include <Vehicle_State.h>

void VehicleState::Vehicle_State_Update(ControllerState Controller_State, VehiclePara Vehicle_Para, const double& step)
{
	double F_D, F_Y, F_L;	// 阻力D，侧力Y，升力L
	double M_X, M_Y, M_Z;
	double d_X, d_Y, d_Z, d_V, d_Gamma, d_Chi, d_Alpha,d_Beta, d_Mu, d_p, d_q, d_r, d_Var_theta, d_Psi, d_Phi_b;
	vector<double> Delta(3, 0.0), M_c(3, 0.0);
	double Mass, Jx, Jy, Jz, C, B, S;
	
	// 读飞行器参数
	Mass = Vehicle_Para.Mass;
	Jx = Vehicle_Para.Jx;
	Jy = Vehicle_Para.Jy;
	Jz = Vehicle_Para.Jz;
	B = Vehicle_Para.B;
	C = Vehicle_Para.C;
	S = Vehicle_Para.S;
	
	// 读控制器状态
	Delta[0] = Controller_State.Delta[0];
	Delta[1] = Controller_State.Delta[1];
	Delta[2] = Controller_State.Delta[2];

	// 读力矩
	M_c[0] = Controller_State.M_c[0];
	M_c[1] = Controller_State.M_c[1];
	M_c[2] = Controller_State.M_c[2];

	// 所有单位原始情况都是弧度制，公式中需要的是角度制
	vector<double> C_FM(6);
	getCoefficients(Ma, Alpha, Beta, Delta[0], Delta[1], Delta[2], C, B, V, p, q, r, C_FM);
	F_D = C_FM[0] * 0.5 * rho * V * V * S;
	F_L = C_FM[1] * 0.5 * rho * V * V * S;
	F_Y = C_FM[2] * 0.5 * rho * V * V * S;
	// 计算力矩
	M_X = C_FM[3] * 0.5 * rho * V * V * S * B;
	M_Y = C_FM[4] * 0.5 * rho * V * V * S * C;
	M_Z = C_FM[5] * 0.5 * rho * V * V * S * B;

	printf("马赫:%f 攻角:%f 侧滑角:%f 倾侧角:%f\n", Ma, Alpha*rad, Beta*rad, Mu*rad);
	//printf("X:%f Y:%f H:%f\n",X, Y, -Z);
	//printf("p:%f q:%f r:%f \n", p, q, r);
	//printf("左升降副翼:%f 右升降副翼:%f 方向舵:%f \n", Delta[0]*rad, Delta[1]*rad, Delta[2]*rad);
	
	
	// 飞行器运动非线性方程组
	d_X     = V * cos(Chi) * cos(Gamma);
	d_Y     = V * sin(Chi) * cos(Gamma);
	d_Z     = -V * sin(Gamma);
	d_V     = (-F_D - Mass * g * sin(Gamma)) / Mass;
	d_Gamma = (F_L * cos(Mu) - F_Y * sin(Mu) - Mass * g * cos(Gamma)) / (Mass * V);
	d_Chi   = (F_L * sin(Mu) + F_Y * cos(Mu)) / (Mass * V * cos(Gamma));
	d_Alpha = q - tan(Beta) * (p * cos(Alpha) + r * sin(Alpha)) + (-F_L  + Mass * g * cos(Gamma) * cos(Mu)) / (Mass * V * cos(Beta));
	d_Beta  = p * sin(Alpha) - r * cos(Alpha) + (F_Y + Mass * g * cos(Gamma) * sin(Mu)) / (Mass * V);
	d_Mu    = (p * cos(Alpha) + r * sin(Alpha)) / cos(Beta) + (F_L * tan(Beta) + F_L * tan(Gamma) * sin(Mu) + F_Y * tan(Gamma) * cos(Mu) - Mass * g * cos(Gamma) * cos(Mu) * tan(Beta)) / (Mass * V);
	d_p     = M_X / Jx + q * r * (Jy - Jz) / Jx;
	d_q     = M_Y / Jy + p * r * (Jz - Jx) / Jy;
	d_r     = M_Z / Jz + p * q * (Jx - Jy) / Jz;

	//printf("%f %f %f %f %f %f %f %f %f %f %f %f\n", d_X, d_Y, d_Z, d_V, d_Gamma, d_Chi, d_Alpha, d_Beta, d_Mu, d_p, d_q, d_r);

	// 地面坐标系和机体坐标系下的俯仰角，偏航角，滚转角对控制没有用
	d_Var_theta = q * sin(Phi_b) / cos(Var_theta) + r * cos(Phi_b) / cos(Var_theta);
	d_Psi = q * cos(Phi_b) - r * sin(Phi_b);
	d_Phi_b = p + q * tan(Var_theta) * sin(Phi_b) + r * tan(Var_theta) * cos(Phi_b);
	
	// 欧拉法得状态方程数值解
	Time  = Time + step;
	V     = V + d_V * step;
	Alpha = Alpha + d_Alpha * step;
	Beta  = Beta + d_Beta * step;
	Mu    = Mu + d_Mu * step;
	Gamma = Gamma + d_Gamma * step;
	Chi   = Chi + d_Chi * step;
	p     = p + d_p * step;
	q     = q + d_q * step;
	r     = r + d_r * step;
	X     = X + d_X * step;
	Y     = Y + d_Y * step;
	Z     = Z + d_Z * step;
	Var_theta = Var_theta + d_Var_theta * step;
	Psi = Psi + d_Psi * step;
	Phi_b = Phi_b + d_Phi_b * step;

	// 输入高度，利用新高度输出声速a，马赫数Ma，重力加速度g和大气密度rho
	AtmoPara atmo;
	atmo.GetAtmoPara(-Z);

	Ma = V / atmo.a;
	rho = atmo.rho;
	g = atmo.g;

}

VehicleState::VehicleState(double H, double v, double alpha, double beta, double mu)  //飞行器状态初始化
{
	AtmoPara Atmo;

	Atmo.GetAtmoPara(H);

	h = H;					// 初始高度
	V = v;					// 飞行速度
	Time = 0.0;				// 时间

	Alpha = alpha;			// 初始攻角为 2 度
	Beta = beta;			// 初始侧滑角为 1 度
	Mu = mu;				// 初始倾侧角为 1 度

	Gamma = 0.0;			// 弹道倾角
	Chi = 0.0;				// 弹道偏角

	p = 0.0 / rad;			// 滚转角速度
	q = 0.0 / rad;			// 俯仰角速度
	r = 0.0 / rad;			// 偏航角速度

	X = 0.0;
	Y = 0.0;
	Z = -h;					// 地面坐标系下起始位置高度为-h

	Var_theta = 0.0;		// 俯仰角
	Psi = 0.0;				// 偏航角
	Phi_b = 0.0;			// 滚转角

	Ma = V / Atmo.a;		// 马赫数
	rho = Atmo.rho;			// 大气密度
	g = Atmo.g;				// 重力加速度

}