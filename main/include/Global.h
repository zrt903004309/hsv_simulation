#ifndef _INCLUDE_H
#define _INCLUDE_H

//宏定义，飞行器状态、参数，控制器状态结构的声明，全局变量rad和PI的声明
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

#define controlMode "全局pid滑模"
//#define controlMode "非线性全局pid滑模"
//#define controlMode "非线性pid全局自适应滑模"
//#define controlMode "非线性pid全局自适应滑模容错"
//#define controlMode "基于DANN干扰观测器的全局pid滑模"
//#define controlMode "基于DANN干扰观测器的非线性全局pid滑模"
//#define controlMode "基于DANN干扰观测器的非线性pid全局自适应滑模"

#define situation "舵面损伤"

struct ModelConfig
{
	const double T = 30;										// 仿真时间
	const double step = 0.001;									// 仿真步长
	const unsigned int iters = (int)(T / step);					// 迭代次数
	const int tctrl = 5;										// 控制器更新周期 * h

	int disturb_flag = 1;										// 是否加入干扰

	int hardware_en = 0;										// 是否加入硬件仿真
	int network_en = 1;											// 是否加入神经网络

	const double h = 35000.0;									// 高度 - 再入段一般在26km到40km
	const double v = 3000.0;									// 速度

	const double alpha = 2.0 / rad;								// 初始攻角
	const double beta = 1.0 / rad;								// 初始侧滑角
	const double mu = 2.0 / rad;								// 初始倾侧角


	const double dd_alpha_ref = 0.0 / rad;						// 期望攻角角加速度
	const double dd_beta_ref = 0.0 / rad;						// 期望侧滑角角加速度
	const double dd_mu_ref = 0.0 / rad;							// 期望倾侧角角加速度

	const double delta_limit = 30.0;							// 舵面角限幅,用deg

	const double d_delta_limit = 200.0 * step * tctrl ;			// 舵面角角速度限幅 200deg/s

	const double DecreaseFactor = 1.0;							// 效益损失系数
	const double BiasFactor = 5;								// 舵面恒偏差值

	// 跟踪的姿态角轨迹模式与赋值
	int trajectoryMode = 1;										// 0:直接跟踪固定值 阶跃响应 1:跟踪斜坡信号
	const double alpha_ref = 5.0 / rad;							// 期望攻角
	const double beta_ref = 0.0 / rad;							// 期望侧滑角
	const double mu_ref = -1.0 / rad;							// 期望倾侧角

	std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_vehicle.txt";
	std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(DecreaseFactor * 100)) + "_control.txt";
	//std::string vehicle_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_vehicle.txt";
	//std::string control_filename = "../data/" + std::to_string(int(h / 1000)) + "_" + std::to_string(int(BiasFactor * 10)) + "_control.txt";

	//const char* vehicle_filename = "../data/result_vehicle.txt";
	//const char* control_filename = "../data/result_control.txt";
};
//struct ModelConfig
//{
//	const double T = 20;								// 仿真时间
//	const double step = 0.001;							// 仿真步长
//	const unsigned int iters = (int)(T / step);			// 迭代次数
//	const int tctrl = 5;								// 控制器更新周期 * h
//
//	const double h = 40000.0;							// 高度
//	const double v = 3000.0;							// 速度
//	const double alpha = 2.0 / rad;						// 初始攻角
//	const double beta = 1.0 / rad;						// 初始侧滑角
//	const double mu = 1.0 / rad;						// 初始倾侧角
//
//	const double alpha_ref = 3.0 / rad;					// 期望攻角
//	const double beta_ref = 0.0 / rad;					// 期望侧滑角
//	const double mu_ref = 2.0 / rad;					// 期望倾侧角
//
//	const double dd_alpha_ref = 0.0 / rad;				// 期望攻角角加速度
//	const double dd_beta_ref = 0.0 / rad;				// 期望侧滑角角加速度
//	const double dd_mu_ref = 0.0 / rad;					// 期望倾侧角角加速度
//
//	const double delta_a_limit = 30.0 / rad;			// 舵面角限幅
//	const double delta_r_limit = 30.0 / rad;			// 舵面角限幅
//	const double delta_e_limit = 30.0 / rad;			// 舵面角限幅
//
//	const double d_delta_a_limit = 200.0 / rad;			// 舵面角角速度限幅 200deg/s
//	const double d_delta_r_limit = 200.0 / rad;			// 舵面角角速度限幅
//	const double d_delta_e_limit = 200.0 / rad;			// 舵面角角速度限幅
//
//
//	const double DecreaseFactor = 1;					// 效益损失系数
//	const double BiasFactor = 5;						// 舵面恒偏差值
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
	//	// 因为本文研究的是无动力再入段飞行器，所以油箱已空，质量、质心和转动惯量等参数可以认为是不变的
	//	Mass = 400;		// 质量 - 根据《Hypersonic Vehicle Simulation Model:Winged-Cone Configuration》的 Mass Model ，质量为 300,000 磅/ 136,078 kg(满油箱) 递减到 140,000 磅/ 63,500 kg(空油箱)

	//	Jx = 150;		// 三轴惯量，对照configuration中的图
	//	Jy = 1800;
	//	Jz = 1800;

	//	B = 4;			// 翼展 - 60 ft
	//	C = 4;			// 平均气动弦长 - 80 ft
	//	S = 1;			// 特征面积 - 3603.00 平方英尺
	//	Xcg = 4.467;
	//};
	VehiclePara() {
		// 因为本文研究的是无动力再入段飞行器，所以油箱已空，质量、质心和转动惯量等参数可以认为是不变的
		Mass = 63500;		// 质量 - 根据《Hypersonic Vehicle Simulation Model:Winged-Cone Configuration》的 Mass Model ，质量为 300,000 磅/ 136,078 kg(满油箱) 递减到 140,000 磅/ 63,500 kg(空油箱)

		Jx = 569443.39;		// 三轴惯量，对照configuration中的图
		Jy = 5694433.92;
		Jz = 5694433.92;

		B = 18.288;			// 翼展 - 60 ft
		C = 24.384;			// 平均气动弦长 - 80 ft
		S = 334.73;			// 特征面积 - 3603.00 平方英尺
		Xcg = 4.4668;		// 力矩中心到质心的距离
	};

	double Mass;		// 质量

	double Jx;			// 三轴惯量
	double Jy;
	double Jz;

	double B;			// 飞行器翼展长度(span) - 计算偏航力矩和滚动力矩时的特征长度
	double C;			// 平均气动弦长(Mean aerodynamic chord) - 计算俯仰力矩时的特征长度，和特征面积与翼展相关(矩形翼情况一致)
	double S;			// 特征面积

	double Xcg;			// 由于再入段油箱耗尽，力矩中心到质心的距离
};

struct AtmoPara
{
	void GetAtmoPara(double h);

	double g;			// 当前高度的重力加速度
	double rho;			// 当前高度的大气密度
	double a;			// 当前高度的音速
};

#endif