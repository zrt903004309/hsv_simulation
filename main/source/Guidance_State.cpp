#include <Guidance_State.h>

void GuidanceState::Guidance_State_Update(const double& time, const ModelConfig& Model_Config)
{
	if (Model_Config.trajectoryMode == 0) {
		alpha_ref = Model_Config.alpha_ref;
		beta_ref = Model_Config.beta_ref;
		mu_ref = Model_Config.mu_ref;

		dd_alpha_ref = 0.0 / rad;
		dd_beta_ref = 0.0 / rad;
		dd_mu_ref = 0.0 / rad;
	}
	else if (Model_Config.trajectoryMode == 1) {
		if (time >= 0 && time <= 1) {
			cur_alpha_ref = 0.0 / rad;
			cur_mu_ref = 0.0 / rad;
		}
		else if (time > 1 && time <= 10) {
			cur_alpha_ref = (2 * (time - 1)) / rad;
			cur_mu_ref = 0.0 / rad;
		}
		else if (time > 10 && time <= 25) {
			cur_alpha_ref = 18.0 / rad;
			cur_mu_ref = 0.0 / rad;
		}
		else if (time > 25 && time <= 35) {
			cur_alpha_ref = (43 - time) / rad;
			cur_mu_ref = 0.0 / rad;
		}
		else if (time > 35 && time <= 50) {
			cur_alpha_ref = 8 / rad;
			cur_mu_ref = 0.0 / rad;
		}
		else if (time > 50 && time <= 64) {
			cur_alpha_ref = 8 / rad;
			cur_mu_ref = (50 - time) / rad;
		}
		else if (time > 64 && time <= 80) {
			cur_alpha_ref = 8 / rad;
			cur_mu_ref = -14 / rad;
		}
		else if (time > 80 && time <= 90) {
			cur_alpha_ref = 8 / rad;
			cur_mu_ref = (2 * time - 174) / rad;
		}
		else if (time > 90 && time <= 100) {
			cur_alpha_ref = 8 / rad;
			cur_mu_ref = 6 / rad;
		}
		else {
			cur_alpha_ref = 0.0 / rad;
			cur_mu_ref = 0.0 / rad;
		}

		// 增加滤波环节，时间常数0.4s
		alpha_ref = 0.9876 * pre_alpha_ref + 0.0124 * cur_alpha_ref;
		mu_ref = 0.9876 * pre_mu_ref + 0.0124 * cur_mu_ref;

		dd_alpha_ref = (alpha_ref - pre_alpha_ref) / (Model_Config.step * Model_Config.tctrl);
		dd_mu_ref = (mu_ref - pre_mu_ref) / (Model_Config.step * Model_Config.tctrl);

		pre_alpha_ref = alpha_ref;
		pre_mu_ref = mu_ref;
	}
}

GuidanceState::GuidanceState()		// 跟踪姿态初始化
{
	alpha_ref = 0.0 / rad;
	beta_ref = 0.0 / rad;
	mu_ref = 0.0 / rad;

	pre_alpha_ref = 0.0 / rad;
	pre_beta_ref = 0.0 / rad;
	pre_mu_ref = 0.0 / rad;

	cur_alpha_ref = 0.0 / rad;
	cur_beta_ref = 0.0 / rad;
	cur_mu_ref = 0.0 / rad;

	dd_alpha_ref = 0.0 / rad;
	dd_beta_ref = 0.0 / rad;
	dd_mu_ref = 0.0 / rad;

}