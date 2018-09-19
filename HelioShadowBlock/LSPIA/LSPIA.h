#pragma once
#ifndef LSPIA_H
#define LSPIA_H
#include "../Common/CommonFunc.h"

class LSPIA
{
public:
	LSPIA() {};
	void set_datas(vector<MatrixXf*> field_data, vector<MatrixXf*> sample_field_data, vector<MatrixXf*> sd_bk_res);
	vector<vector<MatrixXf*>>& LSPIA_surface(const vector<int>&ctrl_num, const float miu);

private:
	vector<MatrixXf*> field_data;			// 定日镜场所有定日镜x与z坐标
	vector<MatrixXf*> sample_field_data;		// 定日镜场采样定日镜x与z坐标
	vector<vector<MatrixXf*>> sd_bk_res;	// 不同时刻下定日镜场采样定日镜阴影遮挡率

	float BaseFunction(const int i, const int k, const float u, const vector<float>&knot);
	vector<float>& knot_vector(const int k, const vector<float>&param, const int N, const int M);
	float surface_fitting(const MatrixXf*D, MatrixXf*P,
		const vector<MatrixXf*>&Nik, const float miu, const float threashold);
	float surface_adjusting_control_points(const vector<MatrixXf*>&D_even_odd, MatrixXf*P,
		const vector <MatrixXf*>&Nik, const float miu);
};

#endif //LSPIA_H