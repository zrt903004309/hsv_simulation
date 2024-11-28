#ifndef _CoefficientsSixDoF_FandM_H
#define _CoefficientsSixDoF_FandM_H

#include "Global.h"
#include <vector>
using namespace std;

void CoefficientsSixDoF_FandM(const double& Mach, const double& Alpha, const double& Beta, const double& Delta_e, const double& Delta_a, const double& Delta_r, const double& C, const double& B, const double& Vel, const double& omega_x, const double& omega_y, const double& omega_z, vector<double>& y);
#endif
