#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Receiver.h"

class FED
{
public:
	FED(vector<Receiver*>& recvs, const Vector2d& helio_size, double rho, double sigma = 1e-3,
		bool calc_mAA = false, double k = 1.36,  double sigma_sun = 2.24*1e-3) {
		this->recvs = recvs;
		this->helio_size = helio_size;
		this->rho = rho;
		this->calc_mAA = calc_mAA;
		this->k = k;
		this->sigma = sigma;
		this->sigma_sun = sigma_sun;
	}

	double calcFED(const Vector3d& helio_pos);

	Vector2d helio_size;
	vector<Receiver*> recvs;
	bool calc_mAA;
	double k;
	double sigma_sun;
	double sigma;
	double rho;

private:
	double calcAA(const double& dis);
	double calcS(const double& dis);
	double HeightFieldCauchyConvolveLine(double x0, double y0, double xc, double a, double yc, double b, double s);
};

