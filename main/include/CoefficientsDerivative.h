#ifndef _CoefficientsDerivative_H
#define _CoefficientsDerivative_H

#include <Global.h>
#include <vector>

// �������
// Mach : ���
// Alpha : ����
// Beta : ����
// Delta_e : ����������Ƕȣ������ƣ�
// Delta_a : ����������Ƕȣ������ƣ�
// Delta_r : �����Ƕȣ������ƣ�
// C : ƽ�������ҳ�
// B : ��չ
// V : �ٶ�
// p : ��ת���ٶ�
// q : �������ٶ�
// r : ƫ�����ٶ�

// �������
// y : ��ά���������зֱ���� CD, CY, CL, Cl, Cm, Cn

// ��ȡ����ϵ���Է�����ƫ�������ϵ��
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
