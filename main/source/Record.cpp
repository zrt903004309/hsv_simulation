#include<Record.h>
//记录数据
void Record(VehicleState Vehicle_State, ControllerState Controller_State, FILE *File_Vehicle, FILE *File_Control)
{
	double X,Y,Z,V,Gamma,Chi,Alpha,Beta,Mu,p,q,r,Time,Var_theta,Psi,Phi_b;
	double Ma;
	//控制变量定义
	vector<double> Delta(3, 0), M_c(3, 0);
	double s_Alpha, s_Beta, s_Mu;

	Time=Vehicle_State.Time;
	V=Vehicle_State.V;
	Alpha=Vehicle_State.Alpha;
	Beta=Vehicle_State.Beta;
	Gamma=Vehicle_State.Gamma;
	Mu=Vehicle_State.Mu;
	Chi=Vehicle_State.Chi;
	p=Vehicle_State.p;
	q=Vehicle_State.q;
	r=Vehicle_State.r;
	X=Vehicle_State.X;
	Y=Vehicle_State.Y;
	Z=Vehicle_State.Z; 
	Var_theta=Vehicle_State.Var_theta;
	Psi=Vehicle_State.Psi;
	Phi_b=Vehicle_State.Phi_b;
	Ma=Vehicle_State.Ma;
	/******************/
	Delta[0] = Controller_State.Delta[0];
	Delta[1] = Controller_State.Delta[1];
	Delta[2] = Controller_State.Delta[2];

	M_c[0] = Controller_State.M_c[0];
	M_c[1] = Controller_State.M_c[1];
	M_c[2] = Controller_State.M_c[2];


	s_Alpha = Controller_State.s_Alpha;
	s_Beta = Controller_State.s_Beta;
	s_Mu = Controller_State.s_Mu;

	fprintf(File_Vehicle,"%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",Time,V,Ma,Alpha*rad,Beta*rad,Mu*rad,X,Y,Z);
	fprintf(File_Control,"%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",Delta[0]*rad,Delta[1]*rad,Delta[2]*rad, s_Alpha, s_Beta, s_Mu, M_c[0], M_c[1], M_c[2]);
}