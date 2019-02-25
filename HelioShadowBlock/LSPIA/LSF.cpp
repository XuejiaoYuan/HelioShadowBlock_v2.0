#include "LSF.h"

void LSF::LSF_surface(SolarScene *solar_scene, unordered_map<int, double>& sdbk)
{

	lsf_data.clear();
	Z.clear();
	lsf_param.clear();
	lsf_data.resize(seg_row, vector<MatrixXd>(seg_col, MatrixXd::Zero(6, 6)));
	Z.resize(seg_row, vector<MatrixXd>(seg_col, MatrixXd::Zero(6, 1)));
	lsf_param.resize(seg_row, vector<MatrixXd>(seg_col, MatrixXd::Zero(6, 1)));

	double row_length = solar_scene->layouts[0]->layout_size.z();
	double col_length = solar_scene->layouts[0]->layout_size.x();
	double row_start = solar_scene->layouts[0]->layout_bound_pos.z();
	double col_start = solar_scene->layouts[0]->layout_bound_pos.x();
	double row_gap = row_length / seg_row;
	double col_gap = col_length / seg_col;

	vector<Heliostat*>& helios = solar_scene->helios;

	for (auto& iter : sdbk) {
		int index = iter.first;
		int i = helios[index]->helio_pos.z() / row_gap;
		int j = helios[index]->helio_pos.x() / col_gap;
		for(int k=-1; k<2; ++k)
			for (int l = -1; l < 2; ++l) {
				if (0 <= i + k && i + k < seg_row && 0 <= j + l && j + l < seg_col) {
					lsf_param[i + k][j + l] += helios[index]->lsf_param_M;
					Z[i + k][j + l] += sdbk[index] * helios[index]->lsf_param_v.transpose();
				}
			}
	}

	for (int i = 0; i < seg_row; ++i)
		for (int j = 0; j < seg_col; ++j)
			lsf_param[i][j] = lsf_data[i][j].inverse()*Z[i][j];

}

// TODO: 计算结果，检验结果