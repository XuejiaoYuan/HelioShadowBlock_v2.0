#pragma once
#ifndef LSPIA_H
#define LSPIA_H
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"

class LSPIA
{
public:
	LSPIA() {};
	void set_datas(vector<MatrixXd*> field_data, vector<MatrixXd*> sample_field_data, vector<MatrixXd*> sd_bk_res);
	vector<vector<MatrixXd*>> LSPIA_surface(const vector<int>&ctrl_num, const double miu);
	void checkFittingData(vector<Heliostat*>& helios, MatrixXd* field_index, vector<vector<MatrixXd*>>& fitting_data);
private:
	vector<MatrixXd*> field_data;			// 定日镜场所有定日镜x与z坐标
	vector<MatrixXd*> sample_field_data;		// 定日镜场采样定日镜x与z坐标
	vector<MatrixXd*> sample_sd_bk_res;		// 不同时刻下定日镜场采样定日镜阴影遮挡率
	MatrixXd* delta_even_uv;
	MatrixXd* delta_odd_uv;
	MatrixXd* delta;

	double BaseFunction(const int i, const int k, const double u, const vector<double>&knot);
	vector<double> knot_vector(const int k, const vector<double>&param, const int N, const int M);
	double surface_fitting(const MatrixXd*D, MatrixXd&P,
		const vector<MatrixXd*>&Nik, const double miu, const double threashold);
	double surface_adjusting_control_points(const vector<MatrixXd*>&D_even_odd, MatrixXd&P,
		const vector <MatrixXd*>&Nik, const double miu);
};

#endif //LSPIA_H