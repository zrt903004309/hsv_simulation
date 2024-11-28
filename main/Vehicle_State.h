#pragma once
#ifndef _VEHICLE_STATE_H
#define _VEHICLE_STATE_H

#include <Global.h>

class ControllerState;
struct VehiclePara;

class VehicleState {
public:

	VehicleState(double H = 30000.0, double v = 3000.0, double alpha = 2.0 / rad, double beta = 1.0 / rad, double mu = 1.0 / rad);

	void Vehicle_State_Update(ControllerState Controller_State, VehiclePara Vehicle_Para, const double& step);

	//void Initial_Vehicle_State();

	double h;			// 初始高度
	double V;			// 飞行速度
	double Time;		// 时间

	double X;
	double Y;
	double Z;			// 地面坐标系下起始位置高度为-h

	double Gamma;		// 弹道倾角
	double Chi;			// 弹道偏角

	double Alpha;		// 初始攻角为 2 度
	double Beta;		// 初始侧滑角为 1 度
	double Mu;			// 初始倾侧角为 1 度

	double p;			// 滚转角速度
	double q;			// 俯仰角速度
	double r;			// 偏航角速度

	double Var_theta;	// 俯仰角
	double Psi;			// 偏航角
	double Phi_b;		// 滚转角

	double Ma, rho, g;	// 马赫数和大气密度(未初始化)

};
#endif // _VEHICLE_STATE_H