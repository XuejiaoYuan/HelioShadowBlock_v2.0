#include "FED.h"

float FED::calcFED(const Vector3f & helio_pos)
{
	Vector4f xc(-helio_size.x() / 2, helio_size.x() / 2, helio_size.x() / 2, -helio_size.x() / 2);
	Vector4f yc(-helio_size.y() / 2, -helio_size.y() / 2, helio_size.y() / 2, helio_size.y() / 2);

	Vector3f recv2helio = helio_pos - recvs[0]->recv_pos;
	float recv2helio_dis = recv2helio.norm();
	recv2helio = recv2helio.normalized();
	float mcosine = 0;
	if (recvs.size() == 1)
		mcosine = recv2helio.transpose() * recvs[0]->recv_normal;
	else {
		for (auto& recv : recvs) {
			float tmp_cos = recv2helio.transpose() * recv->recv_normal;
			if (tmp_cos > mcosine) {
				mcosine = tmp_cos;
			}
		}
	}

	float mAA = 1.0;
	if (calc_mAA)
		mAA = calcAA(recv2helio_dis);

	float s = calcS(recv2helio_dis);
	float ratio = 1.0 / (PI*pow(s, 2));

	float mS = 0;
	for (int i = 0; i < 4; i++) {
		float val = HeightFieldCauchyConvolveLine(0, 0, xc(i), xc((i + 1) % 4) - xc(i), yc(i), yc((i + 1) % 4) - yc(i), s);
		mS += ratio * val;
	}

	return rho*mAA*mcosine*mS;
}


float FED::calcAA(const float& dis) {

	float mAA = 0;
	if (dis <= 1000)
		mAA = 0.99321 - 0.0001176*dis + 1.97 * 1e-8 * pow(dis, 2);
	else
		mAA = exp(-0.0001106 * dis);
	return mAA;
}


float FED::calcS(const float&dis) {

	return k / (dis*pow(sigma_sun, 2) + pow(2 * sigma, 0.5));
}


float FED::HeightFieldCauchyConvolveLine(float x0, float y0, float x1, float a, float y1, float b, float s)
{
	if (b == 0 || s == 0)
		return 0;

	float m = 1 / pow(s, 2);
	float S1 = pow(b, 2) * (m + pow(x0 - x1, 2)) + pow(a, 2) * (m + pow(y0 - y1, 2)) - 2 * a*b*(x0 - x1)*(y0 - y1);
	float S2 = pow(b, 2) * pow(x0 - x1, 2) + pow(a, 2) * (m + pow(y0 - y1, 2)) - 2 * a*b*(x0 - x1)*(y0 - y1);

	if (abs(S1) < 1e-7 || abs(S2) < 1e-7)
		return 0;

	float S3 = b*(x0 - x1) + a*(-y0 + y1);
	float S4 = b*(-x0 + x1) + a*(y0 - y1);

	float S8 = a*(-x0 + x1) + b*(-y0 + y1);
	float S9 = pow(a, 2) + a*(-x0 + x1) + b*(b - y0 + y1);
	float S10 = m + pow(x0 - x1, 2) + pow(y0 - y1, 2);
	float S11 = m + pow(a - (x0 - x1), 2) + pow(b - y0 + y1, 2);

	float F0 = 1 / 4 * b*(2 * S4 * sqrt(S1)* atan(S8 / sqrt(S1)) / (b*m*S2) + 
		(2 * b*S3*atan(S8 / sqrt(S1))) / (sqrt(S1)* S2) + 
		2 * (-y0 + y1)*atan((-x0 + x1) / sqrt(m + pow(-y0 + y1, 2))) / (b*m*sqrt(m + pow(-y0 + y1, 2))) + 
		a*log(S10) / S2 - a*log(pow(b, 2) * S10) / S2);

	float F1 = 1 / 4 * b*(2 * S4 * sqrt(S1)* atan(S9 / sqrt(S1)) / (b*m*S2) + 
		2 * b* S3 * atan(S9 / sqrt(S1)) / (sqrt(S1) * S2) +
		(2 * (b - y0 + y1)*atan((a - x0 + x1) / sqrt(m + pow(b - y0 + y1, 2)))) / (b*m*sqrt(m + pow(b - y0 + y1, 2))) +
		a*log(S11) / S2 - a*log(pow(b, 2) * S11) / S2);

	return F1 - F0;
}
