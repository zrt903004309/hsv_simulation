#ifndef _INCLUDE_H
#define _INCLUDE_H

//�궨�壬������״̬��������������״̬�ṹ��������ȫ�ֱ���rad��PI������
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<cmath>
#include<time.h>
#include<vector>
#include<Eigen/dense>
#include<iostream>
#include <Coefficients.h>
#include <CoefficientsDerivative.h>

#include<torch/torch.h>
#include<torch/script.h>

//#include<SerialPort.h>

#define rad 57.2958
#define PI	3.14159265

#define controlMode "ȫ��pid��ģ"
//#define controlMode "������ȫ��pid��ģ"
//#define controlMode "������pidȫ������Ӧ��ģ"
//#define controlMode "������pidȫ������Ӧ��ģ�ݴ�"
//#define controlMode "����DANN���Ź۲�����ȫ��pid��ģ"
//#define controlMode "����DANN���Ź۲����ķ�����ȫ��pid��ģ"
//#define controlMode "����DANN���Ź۲����ķ�����pidȫ������Ӧ��ģ"

#define situation "��������"

struct ModelConfig
{
	const double T = 30;										// ����ʱ��
	const double step = 0.001;									// ���沽��
	const unsigned int iters = (int)(T / step);					// ��������
	const int tctrl = 5;										// �������������� * h

	int disturb_flag = 1;										// �Ƿ�������

	int hardware_en = 0;										// �Ƿ����Ӳ������
	int network_en = 1;											// �Ƿ����������

	const double h = 35000.0;									// �߶� - �����һ����26km��40km
	const double v = 3000.0;									// �ٶ�

	const double alpha = 2.0 / rad;								// ��ʼ����
	const double beta = 1.0 / rad;								// ��ʼ�໬��
	const double mu = 2.0 / rad;								// ��ʼ����


	const double dd_alpha_ref = 0.0 / rad;						// �������ǽǼ��ٶ�
	const double dd_beta_ref = 0.0 / rad;						// �����໬�ǽǼ��ٶ�
	const double dd_mu_ref = 0.0 / rad;							// �������ǽǼ��ٶ�

	const double delta_limit = 30.0;							// ������޷�,��deg

	const double d_delta_limit = 200.0 * step * tctrl ;			// ����ǽ��ٶ��޷� 200deg/s

	const double DecreaseFactor = 1.0;							// Ч����ʧϵ��
	const double BiasFactor = 5;								// �����ƫ��ֵ

	// ���ٵ���̬�ǹ켣ģʽ�븳ֵ
	int trajectoryMode = 1;										// 0:ֱ�Ӹ��ٹ̶�ֵ ��Ծ��Ӧ 1:����б���ź�
	const double alpha_ref = 5.0 / rad;							// ��������
	const double beta_ref = 0.0 / rad;							// �����໬��
	const double mu_ref = -1.0 / rad;							// ��������

	std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_vehicle.txt";
	std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_control.txt";
	//std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_vehicle.txt";
	//std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_control.txt";

	//const char* vehicle_filename = "../data/result_vehicle.txt";
	//const char* control_filename = "../data/result_control.txt";
};
//struct ModelConfig
//{
//	const double T = 20;								// ����ʱ��
//	const double step = 0.001;							// ���沽��
//	const unsigned int iters = (int)(T / step);			// ��������
//	const int tctrl = 5;								// �������������� * h
//
//	const double h = 40000.0;							// �߶�
//	const double v = 3000.0;							// �ٶ�
//	const double alpha = 2.0 / rad;						// ��ʼ����
//	const double beta = 1.0 / rad;						// ��ʼ�໬��
//	const double mu = 1.0 / rad;						// ��ʼ����
//
//	const double alpha_ref = 3.0 / rad;					// ��������
//	const double beta_ref = 0.0 / rad;					// �����໬��
//	const double mu_ref = 2.0 / rad;					// ��������
//
//	const double dd_alpha_ref = 0.0 / rad;				// �������ǽǼ��ٶ�
//	const double dd_beta_ref = 0.0 / rad;				// �����໬�ǽǼ��ٶ�
//	const double dd_mu_ref = 0.0 / rad;					// �������ǽǼ��ٶ�
//
//	const double delta_a_limit = 30.0 / rad;			// ������޷�
//	const double delta_r_limit = 30.0 / rad;			// ������޷�
//	const double delta_e_limit = 30.0 / rad;			// ������޷�
//
//	const double d_delta_a_limit = 200.0 / rad;			// ����ǽ��ٶ��޷� 200deg/s
//	const double d_delta_r_limit = 200.0 / rad;			// ����ǽ��ٶ��޷�
//	const double d_delta_e_limit = 200.0 / rad;			// ����ǽ��ٶ��޷�
//
//
//	const double DecreaseFactor = 1;					// Ч����ʧϵ��
//	const double BiasFactor = 5;						// �����ƫ��ֵ
//
//	std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_vehicle.txt";
//	std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_control.txt";
//	//std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_vehicle.txt";
//	//std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_control.txt";
//
//	//const char* vehicle_filename = "../data/result_vehicle.txt";
//	//const char* control_filename = "../data/result_control.txt";
//};


struct VehiclePara
{
	//VehiclePara() {
	//	// ��Ϊ�����о������޶�������η����������������ѿգ����������ĺ�ת�������Ȳ���������Ϊ�ǲ����
	//	Mass = 400;		// ���� - ���ݡ�Hypersonic Vehicle Simulation Model:Winged-Cone Configuration���� Mass Model ������Ϊ 300,000 ��/ 136,078 kg(������) �ݼ��� 140,000 ��/ 63,500 kg(������)

	//	Jx = 150;		// �������������configuration�е�ͼ
	//	Jy = 1800;
	//	Jz = 1800;

	//	B = 4;			// ��չ - 60 ft
	//	C = 4;			// ƽ�������ҳ� - 80 ft
	//	S = 1;			// ������� - 3603.00 ƽ��Ӣ��
	//	Xcg = 4.467;
	//};
	VehiclePara() {
		// ��Ϊ�����о������޶�������η����������������ѿգ����������ĺ�ת�������Ȳ���������Ϊ�ǲ����
		Mass = 63500;		// ���� - ���ݡ�Hypersonic Vehicle Simulation Model:Winged-Cone Configuration���� Mass Model ������Ϊ 300,000 ��/ 136,078 kg(������) �ݼ��� 140,000 ��/ 63,500 kg(������)

		Jx = 569443.39;		// �������������configuration�е�ͼ
		Jy = 5694433.92;
		Jz = 5694433.92;

		B = 18.288;			// ��չ - 60 ft
		C = 24.384;			// ƽ�������ҳ� - 80 ft
		S = 334.73;			// ������� - 3603.00 ƽ��Ӣ��
		Xcg = 4.4668;		// �������ĵ����ĵľ���
	};

	double Mass;		// ����

	double Jx;			// �������
	double Jy;
	double Jz;

	double B;			// ��������չ����(span) - ����ƫ�����غ͹�������ʱ����������
	double C;			// ƽ�������ҳ�(Mean aerodynamic chord) - ���㸩������ʱ���������ȣ��������������չ���(���������һ��)
	double S;			// �������

	double Xcg;			// �������������ľ����������ĵ����ĵľ���
};

struct AtmoPara
{
	void GetAtmoPara(double h);

	double g;			// ��ǰ�߶ȵ��������ٶ�
	double rho;			// ��ǰ�߶ȵĴ����ܶ�
	double a;			// ��ǰ�߶ȵ�����
};

#endif