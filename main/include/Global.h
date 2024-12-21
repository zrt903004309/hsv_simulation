#ifndef _INCLUDE_H
#define _INCLUDE_H

//宏定义，飞行器状态、参数，控制器状态结构的声明，全局变量rad和PI的声明
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
#define controlMode "快速滑模"
#define situation "舵面损伤"

#define dsp_type float				// 针对dsp平台串口通讯需要使用float数据类型


struct ModelConfig
{
	const double T = 20;								// 仿真时间
	const double step = 0.001;							// 仿真步长
	const unsigned int iters = (int)(T / step);			// 迭代次数
	const int tctrl = 5;								// 控制器更新周期*h

	const double h = 30000.0;							// 高度
	const double v = 3000.0;							// 速度
	const double alpha = 2.0 / rad;						// 初始攻角
	const double beta = 1.0 / rad;						// 初始侧滑角
	const double mu = 1.0 / rad;						// 初始倾侧角

	const double alpha_ref = 3.0 / rad;					// 期望攻角
	const double beta_ref = 0.0 / rad;					// 期望侧滑角
	const double mu_ref = 2.0 / rad;					// 期望倾侧角

	const double dd_alpha_ref = 0.0 / rad;				// 期望攻角角加速度
	const double dd_beta_ref = 0.0 / rad;				// 期望侧滑角角加速度
	const double dd_mu_ref = 0.0 / rad;					// 期望倾侧角角加速度

	const double delta_a_limit = 30.0 / rad;			// 舵面角限幅
	const double delta_r_limit = 30.0 / rad;			// 舵面角限幅
	const double delta_e_limit = 30.0 / rad;			// 舵面角限幅


	const char* vehicle_filename = "../data/result_vehicle.txt";
	const char* control_filename = "../data/result_control.txt";
};


struct VehiclePara
{
	VehiclePara() {
		Mass = 63500;		// 质量

		Jx = 569443.39;		// 三轴惯量
		Jy = 5694433.92;
		Jz = 5694433.92;

		B = 18.288;
		C = 24.384;
		S = 334.73;			// 特征面积
	};

	double Mass;		// 质量

	double Jx;			// 三轴惯量
	double Jy;
	double Jz;

	double B;			// 飞行器翼展长度 - 计算偏航力矩和滚动力矩时的特征长度
	double C;			// 平均气动弦长 - 计算俯仰力矩时的特征长度
	double S;			// 特征面积
};

struct AtmoPara
{
	void GetAtmoPara(double h);

	double g;			// 当前高度的重力加速度
	double rho;			// 当前高度的大气密度
	double a;			// 当前高度的音速
};

#endif