#ifndef _Coefficients_H
#define _Coefficients_H

#include <Global.h>
#include <vector>
using namespace std;

void getCoefficients(const double& Mach, const double& Alpha, const double& Beta, const double& Delta_e, const double& Delta_a, const double& Delta_r, const double& C, const double& B, const double& Vel, const double& omega_x, const double& omega_y, const double& omega_z, vector<double>& y);

void getCoefficientsDetails(const double& Mach, const double& Alpha, const double& Beta, const double& Delta_e, const double& Delta_a, const double& Delta_r, const double& C, const double& B, const double& Vel, const double& omega_x, const double& omega_y, const double& omega_z, vector<vector<double>>& y);

void loadAeroCoefficients(const double& Mach_step, const double& Alpha_step);

#endif
