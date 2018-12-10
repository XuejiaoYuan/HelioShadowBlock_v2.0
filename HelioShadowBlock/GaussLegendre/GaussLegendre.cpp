#include "GaussLegendre.h"


double GaussLegendre::calcInte(const Vector4d & x, const Vector4d & y, const double sigma, const double ratio)
{
	Vector2d map_v;
	double sum = 0.0;
	for (int i = 0; i < weight[0].size(); i++) {
		for (int j = 0; j < weight[1].size(); j++) {
			map_v = map(x, y, node[0][i], node[1][j]);
			sum += weight[0][i] * weight[1][j] * jacobi(x, y, node[0](i), node[1](j))*flux_func(map_v.x(), map_v.y(), sigma, ratio);
		}
	}
	return sum;
}

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
//��������������x1��x2������n�����س���Ϊn��x_list��w_list��
//���зֱ���n��gauss - legendre����ֹ�ʽ������㼰Ȩ��
//	:param x1 : ��������
//	:param x2 : ��������
//	:param x : ����ֹ�ʽ�������
//	:param w : ����ֹ�ʽ��Ȩ��
//	:param n : ��˹���ֽ���
///
void GaussLegendre::CalcWeight(const int m, VectorXd& x, VectorXd&w, const double a, const double b)
{
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
				x[i] = -t;
				x[m - 1 - i] = t;
				w[i] = 2.0 / (1.0 - t*t) / (dp*dp);
				w[m - i - 1] = w[i];
				break;
			}
		}
	}
}

inline double GaussLegendre::jacobi(const Vector4d& x, const Vector4d& y, double s, double t) {
	double J00 = -(1.0 - t)*x(0) + (1.0 - t)*x(1) + (1.0 + t)*x(2) - (1.0 - t)*x(3);
	double J01 = -(1.0 - t)*y(0) + (1.0 - t)*y(1) + (1.0 + t)*y(2) - (1.0 - t)*y(3);
	double J10 = -(1.0 - s)*x(0) - (1.0 + s)*x(1) + (1.0 + s)*x(2) + (1.0 - s)*x(3);
	double J11 = -(1.0 - s)*y(0) - (1.0 + s)*y(1) + (1.0 + s)*y(2) + (1.0 - s)*y(3);
	return (J00*J11 - J01*J10) / 16.0;
}

inline Vector2d GaussLegendre::map(const Vector4d&x, const Vector4d&y, double s, double t) {
	Vector4d N;
	N(0) = (1.0 - s)*(1.0 - t) / 4.0;
	N(1) = (1.0 + s)*(1.0 - t) / 4.0;
	N(2) = (1.0 + s)*(1.0 + t) / 4.0;
	N(3) = (1.0 - s)*(1.0 + t) / 4.0;
	Vector2d map_v;
	map_v.x() = N.dot(x);
	map_v.y() = N.dot(y);
	return map_v;
}

