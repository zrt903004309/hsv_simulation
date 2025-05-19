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
		if (time >= 0 && time <= 5) {
			cur_alpha_ref = 0.0 / rad;
			cur_mu_ref = -1.0 / rad;
		}
		else if (time > 5 && time <= 10) {
			cur_alpha_ref = (2 * (time - 5)) / rad;
			cur_mu_ref = -1.0 / rad;
		}
		else if (time > 10 && time <= 20) {
			cur_alpha_ref = 10.0 / rad;
			cur_mu_ref = -1.0 / rad;
		}
		else if (time > 20 && time <= 25) {
			cur_alpha_ref = (30 - time) / rad;
			cur_mu_ref = -1.0 / rad;
		}
		else if (time > 25 && time <= 50) {
			cur_alpha_ref = 5 / rad;
			cur_mu_ref = -1.0 / rad;
		}
		else {
			cur_alpha_ref = 0.0 / rad;
			cur_mu_ref = 0.0 / rad;
		}

		// 增加滤波环节，时间常数0.4s
		alpha_ref = 0.9876 * pre_alpha_ref + 0.0124 * cur_alpha_ref;
		mu_ref = 0.9876 * pre_mu_ref + 0.0124 * cur_mu_ref;
		/*alpha_ref = cur_alpha_ref;
		mu_ref = cur_mu_ref;*/

		d_alpha_ref = (alpha_ref - pre_alpha_ref) / (Model_Config.step * Model_Config.tctrl);
		d_mu_ref = (mu_ref - pre_mu_ref) / (Model_Config.step * Model_Config.tctrl);

		dd_alpha_ref = (d_alpha_ref - pre_d_alpha_ref) / (Model_Config.step * Model_Config.tctrl);
		dd_mu_ref = (d_mu_ref - pre_d_mu_ref) / (Model_Config.step * Model_Config.tctrl);

		pre_alpha_ref = alpha_ref;
		pre_mu_ref = mu_ref;

		pre_d_alpha_ref = d_alpha_ref;
		pre_d_mu_ref = d_mu_ref;
	}
}

GuidanceState::GuidanceState()		// 跟踪姿态初始化
{
	alpha_ref = 0.0 / rad;
	beta_ref = 0.0 / rad;
	mu_ref = 0.0 / rad;

	pre_alpha_ref = 0.0 / rad;
	pre_beta_ref = 0.0 / rad;
	pre_mu_ref = -1.0 / rad;

	cur_alpha_ref = 0.0 / rad;
	cur_beta_ref = 0.0 / rad;
	cur_mu_ref = 0.0 / rad;

	d_alpha_ref = 0.0 / rad;
	d_beta_ref = 0.0 / rad;
	d_mu_ref = 0.0 / rad;

	pre_d_alpha_ref = 0.0 / rad;
	pre_d_beta_ref = 0.0 / rad;
	pre_d_mu_ref = 0.0 / rad;

	dd_alpha_ref = 0.0 / rad;
	dd_beta_ref = 0.0 / rad;
	dd_mu_ref = 0.0 / rad;

}