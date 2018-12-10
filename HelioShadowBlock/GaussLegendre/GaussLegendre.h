#pragma once
#include "../Common/CommonFunc.h"

class GaussLegendre {
public:
	void calcNodeWeight(const int N, const int M) {
		node.clear();
		weight.clear();
		node.resize(2);
		weight.resize(2);
		CalcWeight(N, node[0], weight[0]);
		CalcWeight(M, node[1], weight[1]);
	}

	vector<VectorXd> getX() { return node; }
	vector<VectorXd> getW() { return weight; }
	double calcInte(const Vector4d& x, const Vector4d& y, const double sigma, const double ratio);
	void legendre(const double t, const double m, double&p, double& dp);
	void CalcWeight(const int n, VectorXd& x, VectorXd&w, const double a = -1, const double b = 1);
	double jacobi(const Vector4d& x, const Vector4d& y, double s, double t);
	Vector2d map(const Vector4d&x, const Vector4d&y, double s, double t);
	double flux_func(double x, double y, const double sigma, const double ratio) {
		return exp(-0.5 / pow(sigma, 2)*(pow(x, 2) + pow(y * ratio, 2)));
	}

	vector<VectorXd> node;
	vector<VectorXd> weight;
};