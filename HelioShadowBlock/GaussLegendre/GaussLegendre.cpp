#include "GaussLegendre.h"


void GaussLegendre::legendre(const double t, const double m, double&p, double& dp)
{
	double p0 = 1.0;
	double p1 = t;
	for (int k = 1; k < m; k++) {
		p = ((2.0*k + 1)*t*p1 - k*p0) / (1.0 + k);
		p0 = p1;
		p1 = p;
	}
	dp = m*(p0 - t*p1) / (1.0 - t*t);
}

///
//给定积分上下限x1和x2及阶数n，返回长度为n的x_list及w_list，
//其中分别存放n点gauss - legendre求积分公式的坐标点及权重
//	:param x1 : 积分下限
//	:param x2 : 积分上限
//	:param x : 求积分公式的坐标点
//	:param w : 求积分公式的权重
//	:param n : 高斯积分阶数
///
void GaussLegendre::CalcWeight(const int m, VectorXd& x, VectorXd&w, const double a, const double b)
{
	//int m = (n+1) / 2;
	//x.resize(n);
	//w.resize(n);
	//double xm = 0.5*(b + a);
	//double xl = 0.5*(b - a);
	//const double EPS = 3e-11;
	//double z, z1, pp, p1, p2, p3;
	//for (int i = 0; i < m; i++) {
	//	z = cos(PI*(i + 0.75) / (n + 0.5));
	//	z1 = -z;
	//	while (abs(z - z1) > EPS) {
	//		p1 = 1.0;
	//		p2 = 0.0;
	//		for (int j = 0; j < n; j++) {
	//			p3 = p2;
	//			p2 = p1;
	//			p1 = ((2 * j + 1)*z*p2 - j*p3) / (j + 1);
	//		}
	//		pp = n*(z*p1 - p2) / (z*z - 1);
	//		z1 = z;
	//		z = z1 - p1 / pp;
	//	}
	//	x(i) = xm - xl*z;
	//	x(n - 1 - i) = xm + xl*z;
	//	w(i) = 2 * xl / ((1 - z*z)*pp*pp);
	//	w(n - 1 - i) = w(i);
	//}
	x.resize(m);
	w.resize(m);
	int nRoots = int((m + 1) / 2);
	double p, dp, dt, t;
	for (int i = 0; i < nRoots; i++) {
		t = cos(PI*(i + 0.75) / (m + 0.5));
		while (true) {
			legendre(t, m, p, dp);
			dt = -p / dp;
			t += dt;
			if (abs(dt) < Epsilon) {
				x[i] = t;
				x[m - 1 - i] = -t;
				w[i] = 2.0 / (1.0 - t*t) / (dp*dp);
				w[m - i - 1] = w[i];
				break;
			}
		}
	}
}
