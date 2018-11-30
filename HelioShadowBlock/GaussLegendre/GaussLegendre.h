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
private:
	void legendre(const double t, const double m, double&p, double& dp);
	void CalcWeight(const int n, VectorXd& x, VectorXd&w, const double a = -1, const double b = 1);

	vector<VectorXd> node;
	vector<VectorXd> weight;
};