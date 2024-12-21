#ifndef _INCLUDE_H
#define _INCLUDE_H

//�궨�壬������״̬��������������״̬�ṹ��������ȫ�ֱ���rad��PI������
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<vector>
#include<Eigen/dense>
#include<iostream>
#include "CoefficientsSixDoF_FandM.h"

#define rad 57.2958
#define PI	3.14159265
#define controlMode "���ٻ�ģ"
#define situation "��������"

#define dsp_type float				// ���dspƽ̨����ͨѶ��Ҫʹ��float��������


struct ModelConfig
{
	const double T = 20;								// ����ʱ��
	const double step = 0.001;							// ���沽��
	const unsigned int iters = (int)(T / step);			// ��������
	const int tctrl = 5;								// ��������������*h

	const double h = 30000.0;							// �߶�
	const double v = 3000.0;							// �ٶ�
	const double alpha = 2.0 / rad;						// ��ʼ����
	const double beta = 1.0 / rad;						// ��ʼ�໬��
	const double mu = 1.0 / rad;						// ��ʼ����

	const double alpha_ref = 3.0 / rad;					// ��������
	const double beta_ref = 0.0 / rad;					// �����໬��
	const double mu_ref = 2.0 / rad;					// ��������

	const double dd_alpha_ref = 0.0 / rad;				// �������ǽǼ��ٶ�
	const double dd_beta_ref = 0.0 / rad;				// �����໬�ǽǼ��ٶ�
	const double dd_mu_ref = 0.0 / rad;					// �������ǽǼ��ٶ�

	const double delta_a_limit = 30.0 / rad;			// ������޷�
	const double delta_r_limit = 30.0 / rad;			// ������޷�
	const double delta_e_limit = 30.0 / rad;			// ������޷�


	const char* vehicle_filename = "../data/result_vehicle.txt";
	const char* control_filename = "../data/result_control.txt";
};


struct VehiclePara
{
	VehiclePara() {
		Mass = 63500;		// ����

		Jx = 569443.39;		// �������
		Jy = 5694433.92;
		Jz = 5694433.92;

		B = 18.288;
		C = 24.384;
		S = 334.73;			// �������
	};

	double Mass;		// ����

	double Jx;			// �������
	double Jy;
	double Jz;

	double B;			// ��������չ���� - ����ƫ�����غ͹�������ʱ����������
	double C;			// ƽ�������ҳ� - ���㸩������ʱ����������
	double S;			// �������
};

struct AtmoPara
{
	void GetAtmoPara(double h);

	double g;			// ��ǰ�߶ȵ��������ٶ�
	double rho;			// ��ǰ�߶ȵĴ����ܶ�
	double a;			// ��ǰ�߶ȵ�����
};

#endif