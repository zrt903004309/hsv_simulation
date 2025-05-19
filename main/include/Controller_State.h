#pragma once
#ifndef _CONTROLLER_STATE_H
#define _CONTROLLER_STATE_H

#include <Global.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

#include <SerialPort.h>

class VehicleState;
class GuidanceState;
struct VehiclePara;
struct ModelConfig;

class ControllerState {
public:

	ControllerState();

	void Controller_State_Update(const VehicleState& Vehicle_State, const GuidanceState Guidance_State, const VehiclePara& Vehicle_Para, ModelConfig& Model_Config, torch::jit::script::Module Module, CSerialPort& mySerialPort);
	

	// Pre��ʾ��һʱ�� 
	double e_Alpha_Pre;					// ���ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Alpha_Pre;				// ���ǵ�ƫ��ĵ���
	double e_Beta_Pre;					// �໬�ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Beta_Pre;					// �໬�ǵ�ƫ��ĵ���
	double e_Mu_Pre;					// ���ǵ�ƫ��ο�ֵ-ʵ��ֵ��
	double de_Mu_Pre;					// ���ǵ�ƫ��ĵ���

	double dde_Alpha_Pre;				// ���ǵ�ƫ��ĵ���
	double dde_Beta_Pre;					// �໬�ǵ�ƫ��ĵ���
	double dde_Mu_Pre;					// ���ǵ�ƫ��ĵ���

	double e_Alpha_0;					// ���ǵ�ƫ���ֵ
	double e_Beta_0;					// �໬�ǵ�ƫ���ֵ
	double e_Mu_0;						// ���ǵ�ƫ���ֵ

	double s_Alpha;						// ���ǵĻ�ģ����
	double s_Beta;						// �໬�ǵĻ�ģ����
	double s_Mu;						// ���ǵĻ�ģ����

	double s_Alpha_0;					// ���ǵĻ�ģ������ֵ
	double s_Beta_0;					// �໬�ǵĻ�ģ������ֵ
	double s_Mu_0;						// ���ǵĻ�ģ������ֵ

	std::vector<double> M;					// ��������
	std::vector<double> Mc_e;				// ������������
	std::vector<double> Mc;					// ʵ�ʿ�������
	std::vector<double> Delta;				// ����Ƕ�
	std::vector<double> u;				    // ������

	double i_e_Alpha;					// �������Ļ��֣�����ʱ��Ϊ�������ڣ�
	double i_e_Beta;					// �໬�����Ļ��֣�����ʱ��Ϊ�������ڣ�
	double i_e_Mu;						// �������Ļ��֣�����ʱ��Ϊ�������ڣ�

	double A_k, B_k, M_k;
	double A_eta, B_eta, M_eta;

	double pA, pB, pM;

	double i_s5_Alpha, i_s5_Beta, i_s5_Mu;
	double i_s6_Alpha, i_s6_Beta, i_s6_Mu;

	double i_s2_Alpha, i_s2_Beta, i_s2_Mu;
	double i_s3_Alpha, i_s3_Beta, i_s3_Mu;

	// ���Ź۲�����ز�������
	Eigen::Vector3d mf_d;
	Eigen::Vector3d pf_d;
	Eigen::Vector4d phi;
	Eigen::MatrixXd theta;
	Eigen::MatrixXd dtheta;
	Eigen::Matrix4d P;
	Eigen::Matrix4d Q;
	Eigen::Matrix4d dP;
	double lambda, kalman_r, kalman_q;
};


#endif