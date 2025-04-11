#include <Global.h>

void AtmoPara::GetAtmoPara(double h) //获取大气参数，输入为高度，输出为声速a，重力加速度g和大气密度rho
{
	// 包含各种常量，参考杨炳尉 标准大气参数的公式表示，误差在万分之一
	double rho0 = 1.225, P0 = 101325, Re = 6356.766, g0 = 9.8066;
	double W, T, P, H;
	double a, g, rho;

	h = fabs(h) / 1000.0;
	H = h / (1 + h / Re);

	g = g0 * pow(Re / (Re + h), 2);

	if (h <= 11.0191)
	{
		W = 1 - H / 44.3308;
		T = 288.15 * W;
		P = P0 * pow(W, 5.2559);
		rho = rho0 * pow(W, 4.2559);
	}
	else if (h <= 20.0631 && h > 11.0191)
	{
		W = exp((14.9647 - H) / 6.3416);
		T = 216.65;
		P = P0 * W * 0.11953;
		rho = rho0 * W * 0.15898;
	}
	else if (h <= 32.1619 && h > 20.0631)
	{
		W = 1 + (H - 24.9021) / 221.552;
		T = 221.552 * W;
		P = P0 * (2.5158 / 100) * pow(W, -34.1629);
		rho = rho0 * (3.2722 / 100) * pow(W, -35.1629);
	}
	else if (h <= 47.3501 && h > 32.1619)
	{
		W = 1 + (H - 39.7499) / 89.4107;
		T = 250.350 * W;
		P = P0 * (2.8338 / 1000) * pow(W, -12.2011);
		rho = rho0 * (3.2618 / 1000) * pow(W, -13.2011);
	}
	else if (h <= 51.4125 && h > 47.3501)
	{
		W = exp((48.6252 - H) / 7.9223);
		T = 270.650;
		P = P0 * (8.9155 / 10000) * W;
		rho = rho0 * (9.4920 / 10000) * W;
	}
	else if (h <= 71.8020 && h > 51.4125)
	{
		W = 1 - (H - 59.4390) / 88.2218;
		T = 247.021 * W;
		P = P0 * (2.1671 / 10000) * pow(W, 12.2011);
		rho = rho0 * (2.5280 / 10000) * pow(W, 11.2011);
	}
	else if (h <= 86.0 && h > 71.8020)
	{
		W = 1 - (H - 78.0303) / 100.295;
		T = 200.590 * W;
		P = P0 * (1.2274 / 100000) * pow(W, 17.0816);
		rho = rho0 * (1.7632 / 100000) * pow(W, 16.0816);
	}
	else if (h > 86.0)					// 这里不确定91km以上是否成立，原文中只到91km
	{
		W = exp((87.2848 - H) / 5.47);
		T = 186.87;
		P = P0 * ((2.2730 + 1.042 / 1000 * H) / 1000000) * W;
		rho = rho0 * (3.6411 / 1000000) * W;
	}

	a = 20.0468 * sqrt(T);
	this->a = a;
	this->g = g;
	this->rho = rho;

}