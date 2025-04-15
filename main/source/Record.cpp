#include <Record.h>
//记录数据
void Record(const VehicleState& Vehicle_State, const ControllerState& Controller_State, const GuidanceState Guidance_State, FILE *File_Vehicle, FILE *File_Control)
{
	double X, Y, Z, V, Ma, Gamma, Chi, Alpha, Beta, Mu, p, q, r, Time, Var_theta, Psi, Phi_b;
	// 控制变量
	vector<double> Delta(3, 0), M_c(3, 0);
	double s_Alpha, s_Beta, s_Mu;
	// 跟踪姿态
	double alpha_ref, beta_ref, mu_ref;

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

	M_c[0]    = Controller_State.M_c[0];
	M_c[1]    = Controller_State.M_c[1];
	M_c[2]    = Controller_State.M_c[2];

	s_Alpha   = Controller_State.s_Alpha;
	s_Beta    = Controller_State.s_Beta;
	s_Mu      = Controller_State.s_Mu;
	/******************/
	alpha_ref = Guidance_State.alpha_ref;
	beta_ref  = Guidance_State.beta_ref;
	mu_ref    = Guidance_State.mu_ref;

	fprintf(File_Vehicle,"%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
		Time, V, Ma, Alpha*rad, Beta*rad, Mu*rad, X, Y, Z, alpha_ref * rad, beta_ref * rad, mu_ref * rad, Gamma * rad, Chi * rad, p, q, r);

	fprintf(File_Control,"%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f \n",
		Delta[0], Delta[1], Delta[2], s_Alpha, s_Beta, s_Mu, M_c[0], M_c[1], M_c[2]);
}
