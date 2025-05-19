#ifndef _CoefficientsDerivative_H
#define _CoefficientsDerivative_H

#include <Global.h>
#include <vector>

// 输入变量
// Mach : 马赫
// Alpha : 攻角
// Beta : 倾侧角
// Delta_e : 左升降副翼角度（弧度制）
// Delta_a : 右升降副翼角度（弧度制）
// Delta_r : 方向舵角度（弧度制）
// C : 平均气动弦长
// B : 翼展
// V : 速度
// p : 滚转角速度
// q : 俯仰角速度
// r : 偏航角速度

// 输出变量
// y : 二维向量，六行分别代表 CD, CY, CL, Cl, Cm, Cn

// 获取气动系数对方向舵的偏导等相关系数
void getCoefficientsDerivative(const double& Mach, 
							   const double& Alpha, 
							   const double& Beta, 
							   const double& Delta_e, 
							   const double& Delta_a, 
							   const double& Delta_r, 
							   const double& C, 
							   const double& B, 
							   const double& V, 
							   const double& p, 
							   const double& q, 
							   const double& r, 
							   std::vector<std::vector<double>>& y);

#endif
