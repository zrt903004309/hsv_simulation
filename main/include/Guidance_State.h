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

	double alpha_ref;			// 期望攻角
	double beta_ref;			// 期望侧滑角
	double mu_ref;				// 期望倾侧角

	double pre_alpha_ref;		// 上一时刻的期望攻角
	double pre_beta_ref;		// 上一时刻的期望侧滑角
	double pre_mu_ref;			// 上一时刻的期望倾侧角

	double cur_alpha_ref;		// 上一时刻的期望攻角
	double cur_beta_ref;		// 上一时刻的期望侧滑角
	double cur_mu_ref;			// 上一时刻的期望倾侧角

	double d_alpha_ref;			// 期望攻角角加速度
	double d_beta_ref;			// 期望侧滑角角加速度
	double d_mu_ref;			// 期望倾侧角角加速度

	double pre_d_alpha_ref;		// 期望攻角角加速度
	double pre_d_beta_ref;		// 期望侧滑角角加速度
	double pre_d_mu_ref;		// 期望倾侧角角加速度

	double dd_alpha_ref;		// 期望攻角角加速度
	double dd_beta_ref;			// 期望侧滑角角加速度
	double dd_mu_ref;			// 期望倾侧角角加速度
};
#endif // _GUIDANCE_STATE_H