#include <Record.h>
//记录数据
void Record(const VehicleState& Vehicle_State, const ControllerState& Controller_State, const GuidanceState Guidance_State, FILE *File_Vehicle, FILE *File_Control)
{
	double X, Y, Z, V, Ma, Gamma, Chi, Alpha, Beta, Mu, p, q, r, Time, Var_theta, Psi, Phi_b;
	// 控制变量
	std::vector<double> Delta(3, 0), M(3, 0), Mc_e(3, 0), Mc(3, 0), u(3, 0), measure_fd(3,0), predict_fd(3,0);
	double s_Alpha, s_Beta, s_Mu;
	// 跟踪姿态
	double alpha_ref, beta_ref, mu_ref;
	double dd_alpha_ref, dd_beta_ref, dd_mu_ref;
	double A_eta, B_eta, M_eta;

	Time      = Vehicle_State.Time;
	V         = Vehicle_State.V;
	Alpha     = Vehicle_State.Alpha;
	Beta      = Vehicle_State.Beta;
	Gamma     = Vehicle_State.Gamma;
	Mu        = Vehicle_State.Mu;
	Chi       = Vehicle_State.Chi;
	p         = Vehicle_State.p;
	q         = Vehicle_State.q;
	r         = Vehicle_State.r;
	X         = Vehicle_State.X;
	Y         = Vehicle_State.Y;
	Z         = Vehicle_State.Z; 
	Var_theta = Vehicle_State.Var_theta;
	Psi       = Vehicle_State.Psi;
	Phi_b     = Vehicle_State.Phi_b;
	Ma        = Vehicle_State.Ma;
	/******************/
	Delta[0]  = Controller_State.Delta[0];
	Delta[1]  = Controller_State.Delta[1];
	Delta[2]  = Controller_State.Delta[2];
	u[0] = Controller_State.u[0];
	u[1] = Controller_State.u[1];
	u[2] = Controller_State.u[2];

	M[0]      = Controller_State.M[0];
	M[1]      = Controller_State.M[1];
	M[2]      = Controller_State.M[2];
	Mc_e[0]   = Controller_State.Mc_e[0];
	Mc_e[1]   = Controller_State.Mc_e[1];
	Mc_e[2]   = Controller_State.Mc_e[2];
	Mc[0]     = Controller_State.Mc[0];
	Mc[1]     = Controller_State.Mc[1];
	Mc[2]     = Controller_State.Mc[2];

	s_Alpha   = Controller_State.s_Alpha;
	s_Beta    = Controller_State.s_Beta;
	s_Mu      = Controller_State.s_Mu;

	double A_K = Controller_State.A_k;
	double B_K = Controller_State.B_k;
	double M_K = Controller_State.M_k;

	/******************/
	alpha_ref = Guidance_State.alpha_ref;
	beta_ref  = Guidance_State.beta_ref;
	mu_ref    = Guidance_State.mu_ref;

	dd_alpha_ref = Guidance_State.dd_alpha_ref;
	dd_beta_ref = Guidance_State.dd_beta_ref;
	dd_mu_ref = Guidance_State.dd_mu_ref;

	A_eta = Controller_State.A_eta;
	B_eta = Controller_State.B_eta;
	M_eta = Controller_State.M_eta;

	if (controlMode == "基于DANN干扰观测器的非线性pid全局滑模") {
		measure_fd[0] = Controller_State.mf_d(0);
		measure_fd[1] = Controller_State.mf_d(1);
		measure_fd[2] = Controller_State.mf_d(2);

		predict_fd[0] = Controller_State.pf_d(0);
		predict_fd[1] = Controller_State.pf_d(1);
		predict_fd[2] = Controller_State.pf_d(2);
		fprintf(File_Vehicle, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Time, V, Ma, Alpha * rad, Beta * rad, Mu * rad, X, Y, Z, alpha_ref * rad, beta_ref * rad, mu_ref * rad, Gamma * rad, Chi * rad, p, q, r, dd_alpha_ref, dd_beta_ref, dd_mu_ref);

		fprintf(File_Control, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Delta[0], Delta[1], Delta[2], s_Alpha, s_Beta, s_Mu, M[0], M[1], M[2], Mc_e[0], Mc_e[1], Mc_e[2], Mc[0], Mc[1], Mc[2], u[0], u[1], u[2], measure_fd[0], measure_fd[1], measure_fd[2], predict_fd[0], predict_fd[1], predict_fd[2]);
	}
	else if (controlMode == "基于DANN干扰观测器的非线性pid全局自适应滑模")
	{
		measure_fd[0] = Controller_State.mf_d(0);
		measure_fd[1] = Controller_State.mf_d(1);
		measure_fd[2] = Controller_State.mf_d(2);

		predict_fd[0] = Controller_State.pf_d(0);
		predict_fd[1] = Controller_State.pf_d(1);
		predict_fd[2] = Controller_State.pf_d(2);

		fprintf(File_Vehicle, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Time, V, Ma, Alpha * rad, Beta * rad, Mu * rad, X, Y, Z, alpha_ref * rad, beta_ref * rad, mu_ref * rad, Gamma * rad, Chi * rad, p, q, r, dd_alpha_ref, dd_beta_ref, dd_mu_ref);

		fprintf(File_Control, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Delta[0], Delta[1], Delta[2], s_Alpha, s_Beta, s_Mu, M[0], M[1], M[2], Mc_e[0], Mc_e[1], Mc_e[2], Mc[0], Mc[1], Mc[2], u[0], u[1], u[2], measure_fd[0], measure_fd[1], measure_fd[2], predict_fd[0], predict_fd[1], predict_fd[2], A_eta, B_eta, M_eta);
	}
	else {
		fprintf(File_Vehicle, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Time, V, Ma, Alpha * rad, Beta * rad, Mu * rad, X, Y, Z, alpha_ref * rad, beta_ref * rad, mu_ref * rad, Gamma * rad, Chi * rad, p, q, r, dd_alpha_ref, dd_beta_ref, dd_mu_ref);

		fprintf(File_Control, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
			Delta[0], Delta[1], Delta[2], s_Alpha, s_Beta, s_Mu, M[0], M[1], M[2], Mc_e[0], Mc_e[1], Mc_e[2], Mc[0], Mc[1], Mc[2], u[0], u[1], u[2], A_K, B_K, M_K, A_eta, B_eta, M_eta);
	}
}