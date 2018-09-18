#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Receiver.h"

class FED
{
public:
	FED(vector<Receiver*>& recvs, const Vector2f& helio_size, float rho, float sigma = 1e-3,
		bool calc_mAA = false, float k = 1.36,  float sigma_sun = 2.24*1e-3) {
		this->recvs = recvs;
		this->helio_size = helio_size;
		this->rho = rho;
		this->calc_mAA = calc_mAA;
		this->k = k;
		this->sigma = sigma;
		this->sigma_sun = sigma_sun;
	}

	float calcFED(const Vector3f& helio_pos);

	Vector2f helio_size;
	vector<Receiver*> recvs;
	bool calc_mAA;
	float k;
	float sigma_sun;
	float sigma;
	float rho;

private:
	float calcAA(const float& dis);
	float calcS(const float& dis);
	float HeightFieldCauchyConvolveLine(float x0, float y0, float xc, float a, float yc, float b, float s);
};

