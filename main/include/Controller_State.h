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

	// Pre表示上一时刻 
	double e_Alpha_Pre;					// 攻角的偏差（参考值-实际值）
	double de_Alpha_Pre;				// 攻角的偏差的导数
	double e_Beta_Pre;					// 侧滑角的偏差（参考值-实际值）
	double de_Beta_Pre;					// 侧滑角的偏差的导数
	double e_Mu_Pre;					// 倾侧角的偏差（参考值-实际值）
	double de_Mu_Pre;					// 倾侧角的偏差的导数

	double s_Alpha;						// 攻角的滑膜函数
	double s_Beta;						// 侧滑角的滑膜函数
	double s_Mu;						// 倾侧角的滑膜函数

	double s_Alpha_0;					// 攻角的滑膜函数初值
	double s_Beta_0;					// 侧滑角的滑膜函数初值
	double s_Mu_0;						// 倾侧角的滑膜函数初值

	vector<double> M_c;					// 控制力矩
	vector<double> Delta;				// 舵面角度

	double i_e_Alpha;					// 攻角误差的积分（积分时间为控制周期）
	double i_e_Beta;					// 侧滑角误差的积分（积分时间为控制周期）
	double i_e_Mu;						// 倾侧角误差的积分（积分时间为控制周期）

	double i_s5_Alpha, i_s5_Beta, i_s5_Mu;
	double i_s6_Alpha, i_s6_Beta, i_s6_Mu;

	double i_s2_Alpha, i_s2_Beta, i_s2_Mu;
	double i_s3_Alpha, i_s3_Beta, i_s3_Mu;
};


#endif