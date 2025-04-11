#pragma once
#ifndef _CONTROLLER_STATE_H
#define _CONTROLLER_STATE_H

#include <Global.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

class VehicleState;
class GuidanceState;
struct VehiclePara;
struct ModelConfig;

class ControllerState {
public:

	ControllerState();

	void Controller_State_Update(const VehicleState& Vehicle_State, const GuidanceState Guidance_State, const VehiclePara& Vehicle_Para, const ModelConfig& Model_Config);

	// Pre��ʾ��һʱ�� 
	double e_Alpha_Pre;					// ���ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Alpha_Pre;				// ���ǵ�ƫ��ĵ���
	double e_Beta_Pre;					// �໬�ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Beta_Pre;					// �໬�ǵ�ƫ��ĵ���
	double e_Mu_Pre;					// ���ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Mu_Pre;					// ���ǵ�ƫ��ĵ���

	double s_Alpha;						// ���ǵĻ�Ĥ����
	double s_Beta;						// �໬�ǵĻ�Ĥ����
	double s_Mu;						// ���ǵĻ�Ĥ����

	double s_Alpha_0;					// ���ǵĻ�Ĥ������ֵ
	double s_Beta_0;					// �໬�ǵĻ�Ĥ������ֵ
	double s_Mu_0;						// ���ǵĻ�Ĥ������ֵ

	vector<double> M_c;					// ��������
	vector<double> Delta;				// ����Ƕ�

	double i_e_Alpha;					// �������Ļ��֣�����ʱ��Ϊ�������ڣ�
	double i_e_Beta;					// �໬�����Ļ��֣�����ʱ��Ϊ�������ڣ�
	double i_e_Mu;						// �������Ļ��֣�����ʱ��Ϊ�������ڣ�

	double i_s5_Alpha, i_s5_Beta, i_s5_Mu;
	double i_s6_Alpha, i_s6_Beta, i_s6_Mu;

	double i_s2_Alpha, i_s2_Beta, i_s2_Mu;
	double i_s3_Alpha, i_s3_Beta, i_s3_Mu;
};


#endif