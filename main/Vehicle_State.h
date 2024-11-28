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

	double h;			// ��ʼ�߶�
	double V;			// �����ٶ�
	double Time;		// ʱ��

	double X;
	double Y;
	double Z;			// ��������ϵ����ʼλ�ø߶�Ϊ-h

	double Gamma;		// �������
	double Chi;			// ����ƫ��

	double Alpha;		// ��ʼ����Ϊ 2 ��
	double Beta;		// ��ʼ�໬��Ϊ 1 ��
	double Mu;			// ��ʼ����Ϊ 1 ��

	double p;			// ��ת���ٶ�
	double q;			// �������ٶ�
	double r;			// ƫ�����ٶ�

	double Var_theta;	// ������
	double Psi;			// ƫ����
	double Phi_b;		// ��ת��

	double Ma, rho, g;	// ������ʹ����ܶ�(δ��ʼ��)

};
#endif // _VEHICLE_STATE_H