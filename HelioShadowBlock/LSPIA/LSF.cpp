#include "LSF.h"

void LSF::LSF_surface(SolarScene *solar_scene)
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

	double sample_gap = (double)helios.size() / (sample_row*sample_col);

	for(int m=0; m<sample_row;++m)
		for (int n = 0; n < sample_col; ++n) {
			int index = (m*sample_col*sample_gap + n*sample_gap);
			int i = (helios[index]->helio_pos.z() - row_start) / row_gap;
			int j = (helios[index]->helio_pos.x() - col_start) / col_gap;
			for (int k = -1; k<2; ++k)
				for (int l = -1; l < 2; ++l) {
					if (0 <= i + k && i + k < seg_row && 0 <= j + l && j + l < seg_col) {
						lsf_data[i + k][j + l] += helios[index]->lsf_param_M;
						Z[i + k][j + l] += (1-helios[index]->sd_bk)*helios[index]->flux_sum * helios[index]->lsf_param_v.transpose();
					}
				}
		}

#pragma omp parallel for
	for (int i = 0; i < seg_row; ++i)
		for (int j = 0; j < seg_col; ++j)
			lsf_param[i][j] = lsf_data[i][j].inverse()*Z[i][j];

	calcFittingData(helios, row_gap, col_gap, row_start, col_start);
}

double LSF::calcFittingData(vector<Heliostat*>& helios,
					double row_gap, double col_gap, double row_start, double col_start)
{
	double total_flux = 0;
#pragma omp parallel for
	for (int k = 0; k < helios.size(); ++k) {
		if (!helios[k]->fluxCalc) {
			int i = (helios[k]->helio_pos.z() - row_start) / row_gap;
			int j = (helios[k]->helio_pos.x() - col_start) / col_gap;
			helios[k]->total_e = (helios[k]->lsf_param_v*lsf_param[i][j])(0, 0);
		}
#pragma omp atomic
		total_flux += helios[k]->total_e;
	}
	return total_flux;
}

