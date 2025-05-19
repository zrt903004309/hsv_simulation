#pragma once
#ifndef _GUIDANCE_STATE_H
#define _GUIDANCE_STATE_H

#include <Global.h>

struct ModelConfig;

class GuidanceState {
public:

	GuidanceState();

	void Guidance_State_Update(const double& time, const ModelConfig& Model_Config);

	//void Initial_Vehicle_State();

	double alpha_ref;			// ��������
	double beta_ref;			// �����໬��
	double mu_ref;				// ��������

	double pre_alpha_ref;		// ��һʱ�̵���������
	double pre_beta_ref;		// ��һʱ�̵������໬��
	double pre_mu_ref;			// ��һʱ�̵���������

	double cur_alpha_ref;		// ��һʱ�̵���������
	double cur_beta_ref;		// ��һʱ�̵������໬��
	double cur_mu_ref;			// ��һʱ�̵���������

	double d_alpha_ref;			// �������ǽǼ��ٶ�
	double d_beta_ref;			// �����໬�ǽǼ��ٶ�
	double d_mu_ref;			// �������ǽǼ��ٶ�

	double pre_d_alpha_ref;		// �������ǽǼ��ٶ�
	double pre_d_beta_ref;		// �����໬�ǽǼ��ٶ�
	double pre_d_mu_ref;		// �������ǽǼ��ٶ�

	double dd_alpha_ref;		// �������ǽǼ��ٶ�
	double dd_beta_ref;			// �����໬�ǽǼ��ٶ�
	double dd_mu_ref;			// �������ǽǼ��ٶ�
};
#endif // _GUIDANCE_STATE_H