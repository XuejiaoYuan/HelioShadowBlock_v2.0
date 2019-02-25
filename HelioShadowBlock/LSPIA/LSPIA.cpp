#include "LSPIA.h"

//void LSPIA::set_datas(vector<MatrixXd*> field_data, vector<MatrixXd*> sample_field_data, vector<MatrixXd*> sample_sd_bk_res)
//{
//	this->field_data = field_data;
//	this->sample_sd_bk_res = sample_sd_bk_res;
//	this->sample_field_data = sample_field_data;
//}


void LSPIA::setPreDatas(FieldSegment *_field_seg, vector<int>& _ctrl_num, const double _miu)
{
	field_seg = _field_seg;
	ctrl_num_row = MatrixXi::Ones(field_seg->seg_row, field_seg->seg_col);
	ctrl_num_col = MatrixXi::Ones(field_seg->seg_row, field_seg->seg_col);
	int per_row = _ctrl_num[0] / field_seg->seg_row;
	int per_col = _ctrl_num[1] / field_seg->seg_col;
	int rows = field_seg->seg_row;
	int cols = field_seg->seg_col;
	segCnt = MatrixXd::Ones(field_seg->seg_row, field_seg->seg_col);
	
	if(field_seg->seg_row != 1 ||  field_seg->seg_col != 1)
		for (int i = 0; i < field_seg->seg_row; ++i)
			for (int j = 0; j < field_seg->seg_col; ++j) {
				if (i == 0 || i == field_seg->seg_row - 1) ctrl_num_row(i, j) = min(2, field_seg->seg_row) * per_row;
				else ctrl_num_row(i, j) = min(3, field_seg->seg_row) * per_row;
				if (j == 0 || j == field_seg->seg_col - 1) ctrl_num_col(i, j) = min(2, field_seg->seg_col) * per_col;
				else ctrl_num_col(i, j) = min(3, field_seg->seg_col) * per_col;

				if ((i == 0 && j == 0) || (i == 0 && j == cols - 1) || (i == rows - 1 && j == 0) || (i == rows - 1 && j == cols - 1)) segCnt(i, j) = 4;
				else if (i == 0 || i == rows - 1 || j == 0 || j == cols - 1) segCnt(i, j) = 6;
				else segCnt(i, j) = 9;
			}

	miu = _miu;
}

 
void LSPIA::LSPIA_surface()
{
	vector<vector<vector<Vector2i>>>& even_label = field_seg->even_segment_region_label;
	vector<vector<vector<Vector2i>>>& odd_label = field_seg->odd_segment_region_label;
	vector<vector<vector<Vector2i>>>& even_sample_label = field_seg->even_sample_segment_region_label;
	vector<vector<vector<Vector2i>>>& odd_sample_label = field_seg->odd_sample_segment_region_label;
	vector<vector<MatrixXd*>>& even_res = field_seg->even_res;
	vector<vector<MatrixXd*>>& odd_res = field_seg->odd_res;
	vector<vector<MatrixXd*>>& even_sample_res = field_seg->even_sample_res;
	vector<vector<MatrixXd*>>& odd_sample_res = field_seg->odd_sample_res;

	// for test
	vector<MatrixXd*>& even_field_index = field_seg->even_field_index;
	vector<MatrixXd*>& odd_field_index = field_seg->odd_field_index;


	int tCnt = even_sample_res.size();
	int rCnt = even_sample_label.size();
	int seg_row = field_seg->seg_row;
	int seg_col = field_seg->seg_col;
	for (int k = 0; k < rCnt; ++k) {								// 多个不连续区域: Fermat:3; Rect: 1; CrossRect: 1
		MatrixXd* even_sample_pos_x = field_seg->even_sample_pos_x[k];
		MatrixXd* even_sample_pos_y = field_seg->even_sample_pos_y[k];
		MatrixXd* odd_sample_pos_x = field_seg->odd_sample_pos_x[k];
		MatrixXd* odd_sample_pos_y = field_seg->odd_sample_pos_y[k];
		MatrixXd* even_pos_x = field_seg->even_pos_x[k];
		MatrixXd* even_pos_y = field_seg->even_pos_y[k];
		MatrixXd* odd_pos_x = field_seg->odd_pos_x[k];
		MatrixXd* odd_pos_y = field_seg->odd_pos_y[k];
		fstream outFile;
		//fstream outFile("odd_sample.txt", ios_base::out);
		//for (int i = 0; i < odd_sample_pos_x->rows(); ++i) {
		//	for (int j = 0; j < odd_sample_pos_x->cols(); ++j)
		//		outFile << "(" << (*odd_sample_pos_x)(i, j) << ", " << (*odd_sample_pos_y)(i, j) << ")\t";
		//	outFile << endl;
		//}
		//outFile.close();
		//outFile.open("even_sample.txt", ios_base::out);
		//for (int i = 0; i < even_sample_pos_x->rows(); ++i) {
		//	for (int j = 0; j < even_sample_pos_x->cols(); ++j)
		//		outFile << "(" << (*even_sample_pos_x)(i, j) << ", " << (*even_sample_pos_y)(i, j) << ")\t";
		//	outFile << endl;
		//}
		//outFile.close();

		MatrixXd k_even_res = MatrixXd::Zero(even_pos_x->rows(), even_pos_x->cols());
		MatrixXd k_odd_res = MatrixXd::Zero(odd_pos_x->rows(), odd_pos_x->cols());


#pragma omp parallel for
		for (int i = 0; i < seg_row; ++i) {				// 同一场景下的分割区域
			for (int j = 0; j < seg_col; ++j) {
				// 1. 采样子区域划分
				Vector2i even_sample_start = even_sample_label[k][i*seg_col + j][2];
				Vector2i even_sample_len = even_sample_label[k][i*seg_col + j][3];
				Vector2i odd_sample_start = odd_sample_label[k][i*seg_col + j][2];
				Vector2i odd_sample_len = odd_sample_label[k][i*seg_col + j][3];
				MatrixXd sub_even_sample_pos_x = even_sample_pos_x->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_even_sample_pos_y = even_sample_pos_y->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_odd_sample_pos_x = odd_sample_pos_x->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
				MatrixXd sub_odd_sample_pos_y = odd_sample_pos_y->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
				MatrixXd sub_even_sample_res = even_sample_res[0][k]->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_odd_sample_res = odd_sample_res[0][k]->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());

				// 2. 子区域划分
				Vector2i even_start = even_label[k][i*seg_col + j][2];
				Vector2i even_len = even_label[k][i*seg_col + j][3];
				Vector2i even_target_start = even_label[k][i*seg_col + j][4];
				Vector2i even_t_s = even_label[k][i*seg_col + j][0];
				Vector2i even_target_len = even_label[k][i*seg_col + j][1];
				Vector2i odd_start = odd_label[k][i*seg_col + j][2];
				Vector2i odd_len = odd_label[k][i*seg_col + j][3];
				Vector2i odd_target_start = odd_label[k][i*seg_col + j][4];
				Vector2i odd_t_s = odd_label[k][i*seg_col + j][0];
				Vector2i odd_target_len = odd_label[k][i*seg_col + j][1];

				MatrixXd sub_even_pos_x = even_pos_x->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				MatrixXd sub_even_pos_y = even_pos_y->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				MatrixXd sub_odd_pos_x = odd_pos_x->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());
				MatrixXd sub_odd_pos_y = odd_pos_y->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());

				// 3. 采样矩阵参数计算
				vector<int> ctrl_num = { ctrl_num_row(i,j), ctrl_num_col(i, j) };
				int even_sample_row = sub_even_sample_pos_x.rows();
				int sample_col = sub_even_sample_pos_x.cols();
				int odd_sample_row = sub_odd_sample_pos_x.rows();
				double start_x = min(sub_even_sample_pos_x(0, 0), sub_odd_sample_pos_x(0, 0));
				double start_y = min(sub_even_sample_pos_y(even_sample_row - 1, sample_col - 1), sub_odd_sample_pos_y(odd_sample_row - 1, sample_col - 1));
				double end_x = max(sub_even_sample_pos_x(even_sample_row - 1, sample_col - 1), sub_odd_sample_pos_x(odd_sample_row - 1, sample_col - 1));
				double end_y = max(sub_even_sample_pos_y(0, 0), sub_odd_sample_pos_y(0, 0));
				vector<vector<double>> knot_uv = initParameters(
					even_sample_row + odd_sample_row, sample_col, start_x, start_y, end_x, end_y, ctrl_num
				);


				vector<MatrixXd*> Nik = initBaseFunction(knot_uv, &sub_even_sample_pos_x, &sub_even_sample_pos_y, &sub_odd_sample_pos_x, &sub_odd_sample_pos_y, ctrl_num);
				MatrixXd ctrls = initCtrlPoints(&sub_even_sample_res, &sub_odd_sample_res, ctrl_num);

				// 4. 结果矩阵参数计算
				vector<MatrixXd*> f_Nik = initBaseFunction(knot_uv, &sub_even_pos_x, &sub_even_pos_y, &sub_odd_pos_x, &sub_odd_pos_y, ctrl_num);

				for (int t = 0; t < tCnt; ++t) {					// 同一区域不同时刻结果拟合					
					MatrixXd t_sub_even_sample_res = even_sample_res[t][k]->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
					MatrixXd t_sub_odd_sample_res = odd_sample_res[t][k]->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
					double sd_bk_error = surface_fitting(vector<MatrixXd*>{&t_sub_even_sample_res, &t_sub_odd_sample_res}, ctrls, Nik, miu, 0.1);

					// 2. Calculate the fitting surface
					MatrixXd calc_sd_bk_even = (*f_Nik[0]) * ctrls * f_Nik[2]->transpose();
					MatrixXd calc_sd_bk_odd = (*f_Nik[1]) *  ctrls * f_Nik[3]->transpose();


					k_even_res.block(even_t_s.x(), even_t_s.y(), even_target_len.x(), even_target_len.y()) =
						calc_sd_bk_even.block(even_target_start.x(), even_target_start.y(), even_target_len.x(), even_target_len.y());
					k_odd_res.block(odd_t_s.x(), odd_t_s.y(), odd_target_len.x(), odd_target_len.y()) =
						calc_sd_bk_odd.block(odd_target_start.x(), odd_target_start.y(), odd_target_len.x(), odd_target_len.y());
				}
			}
		}
		
		outFile.open("lspia_fitting.txt", ios_base::out);
		for (int i = 0; i < k_even_res.rows(); ++i)
			for (int j = 0; j < k_even_res.cols(); ++j) {
				int index = (*even_field_index[k])(i, j);
				auto h = field_seg->solar_scene->helios[index];
				double gt = h->flux_sum * (1 - h->sd_bk);
				double dis = gt - k_even_res(i, j);
				outFile << (*even_pos_x)(i, j) << ' ' << (*even_pos_y)(i, j) << ' ' << k_even_res(i, j) << ' ' << gt << ' ' <<  dis << ' ' << dis/gt << endl;
			}
		for (int i = 0; i < k_odd_res.rows(); ++i)
			for (int j = 0; j < k_odd_res.cols(); ++j) {
				int index = (*odd_field_index[k])(i, j);
				auto h = field_seg->solar_scene->helios[index];
				double gt = h->flux_sum * (1 - h->sd_bk);
				double dis = gt - k_odd_res(i, j);
				outFile << (*odd_pos_x)(i, j) << ' ' << (*odd_pos_y)(i, j) << ' ' << k_odd_res(i, j) << ' ' << gt << ' ' << dis << ' ' << dis/gt << endl;
			}

		vector<int>& exclude_helio_index = field_seg->solar_scene->layouts[0]->exclude_helio_index[k];
		vector<double>& exclude_helio_res = field_seg->solar_scene->layouts[0]->exclude_helio_res[0][k];
		vector<Heliostat*> helios = field_seg->solar_scene->helios;
		for (int i = 0; i < exclude_helio_res.size(); ++i) {
			Heliostat* helio = helios[exclude_helio_index[i]];
			outFile << helio->helio_pos.x() << ' ' << helio->helio_pos.z() << ' ' << exclude_helio_res[i] << ' ' <<  exclude_helio_res[i] << ' ' << 0 << ' ' << 0 << endl;
		}

		outFile.close();

	}


}


void LSPIA::LSF_surface()
{
	vector<vector<vector<Vector2i>>>& even_label = field_seg->even_segment_region_label;
	vector<vector<vector<Vector2i>>>& odd_label = field_seg->odd_segment_region_label;
	vector<vector<vector<Vector2i>>>& even_sample_label = field_seg->even_sample_segment_region_label;
	vector<vector<vector<Vector2i>>>& odd_sample_label = field_seg->odd_sample_segment_region_label;
	vector<vector<MatrixXd*>>& even_res = field_seg->even_res;
	vector<vector<MatrixXd*>>& odd_res = field_seg->odd_res;
	vector<vector<MatrixXd*>>& even_sample_res = field_seg->even_sample_res;
	vector<vector<MatrixXd*>>& odd_sample_res = field_seg->odd_sample_res;

	// for test
	vector<MatrixXd*>& even_field_index = field_seg->even_field_index;
	vector<MatrixXd*>& odd_field_index = field_seg->odd_field_index;

	int tCnt = even_sample_res.size();
	int rCnt = even_sample_label.size();
	int seg_row = field_seg->seg_row;
	int seg_col = field_seg->seg_col;
	vector<vector<MatrixXd>> res;
	for (int k = 0; k < rCnt; ++k) {								// 多个不连续区域: Fermat:3; Rect: 1; CrossRect: 1
		MatrixXd* even_sample_pos_x = field_seg->even_sample_pos_x[k];
		MatrixXd* even_sample_pos_y = field_seg->even_sample_pos_y[k];
		MatrixXd* odd_sample_pos_x = field_seg->odd_sample_pos_x[k];
		MatrixXd* odd_sample_pos_y = field_seg->odd_sample_pos_y[k];
		MatrixXd* even_pos_x = field_seg->even_pos_x[k];
		MatrixXd* even_pos_y = field_seg->even_pos_y[k];
		MatrixXd* odd_pos_x = field_seg->odd_pos_x[k];
		MatrixXd* odd_pos_y = field_seg->odd_pos_y[k];

		MatrixXd k_even_res = MatrixXd::Zero(even_pos_x->rows(), even_pos_x->cols());
		MatrixXd k_odd_res = MatrixXd::Zero(odd_pos_x->rows(), odd_pos_x->cols());

#pragma omp parallel for
		for (int i = 0; i < seg_row; ++i) {				// 同一场景下的分割区域
			for (int j = 0; j < seg_col; ++j) {
				// 1. 采样子区域划分
				Vector2i even_sample_start = even_sample_label[k][i*seg_col + j][2];
				Vector2i even_sample_len = even_sample_label[k][i*seg_col + j][3];
				Vector2i odd_sample_start = odd_sample_label[k][i*seg_col + j][2];
				Vector2i odd_sample_len = odd_sample_label[k][i*seg_col + j][3];
				MatrixXd sub_even_sample_pos_x = even_sample_pos_x->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_even_sample_pos_y = even_sample_pos_y->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_odd_sample_pos_x = odd_sample_pos_x->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
				MatrixXd sub_odd_sample_pos_y = odd_sample_pos_y->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
				MatrixXd sub_even_sample_res = even_sample_res[0][k]->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
				MatrixXd sub_odd_sample_res = odd_sample_res[0][k]->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());

				// 2. 子区域划分
				Vector2i e_t_s = even_label[k][i*seg_col + j][0];
				Vector2i even_start = even_label[k][i*seg_col + j][2];
				Vector2i even_len = even_label[k][i*seg_col + j][3];
				Vector2i even_target_start = even_label[k][i*seg_col + j][4];
				Vector2i even_target_len = even_label[k][i*seg_col + j][1];
				Vector2i o_t_s = odd_label[k][i*seg_col + j][0];
				Vector2i odd_start = odd_label[k][i*seg_col + j][2];
				Vector2i odd_len = odd_label[k][i*seg_col + j][3];
				Vector2i odd_target_start = odd_label[k][i*seg_col + j][4];
				Vector2i odd_target_len = odd_label[k][i*seg_col + j][1];

				//MatrixXd sub_even_pos_x = even_pos_x->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				//MatrixXd sub_even_pos_y = even_pos_y->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				//MatrixXd sub_odd_pos_x = odd_pos_x->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());
				//MatrixXd sub_odd_pos_y = odd_pos_y->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());
				MatrixXd sub_even_pos_x = even_pos_x->block(e_t_s.x(), e_t_s.y(), even_target_len.x(), even_target_len.y());
				MatrixXd sub_even_pos_y = even_pos_y->block(e_t_s.x(), e_t_s.y(), even_target_len.x(), even_target_len.y());
				MatrixXd sub_odd_pos_x = odd_pos_x->block(o_t_s.x(), o_t_s.y(), odd_target_len.x(), odd_target_len.y());
				MatrixXd sub_odd_pos_y = odd_pos_y->block(o_t_s.x(), o_t_s.y(), odd_target_len.x(), odd_target_len.y());


				for (int t = 0; t < tCnt; ++t) {					// 同一区域不同时刻结果拟合					
					MatrixXd t_sub_even_sample_res = even_sample_res[t][k]->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());

					MatrixXd param = LSF(sub_even_sample_pos_x, sub_even_sample_pos_y, sub_odd_sample_pos_x, sub_odd_sample_pos_y, sub_even_sample_res, sub_odd_sample_res);
					MatrixXd lsf_sub_even_target = LSF_res(param, sub_even_pos_x, sub_even_pos_y);
					MatrixXd lsf_sub_odd_target = LSF_res(param, sub_odd_pos_x, sub_odd_pos_y);

					k_even_res.block(e_t_s.x(), e_t_s.y(), even_target_len.x(), even_target_len.y()) = lsf_sub_even_target;
					k_odd_res.block(o_t_s.x(), o_t_s.y(), odd_target_len.x(), odd_target_len.y()) = lsf_sub_odd_target;
					//k_even_res.block(even_start.x(), even_start.y(), even_len.x(), even_len.y()) = lsf_sub_even_target;
					//k_odd_res.block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y()) = lsf_sub_odd_target;
				}
			}
		}
		//k_even_res.array() /= segCnt.array();
		//k_odd_res.array() /= segCnt.array();

		fstream outFile("lsf_fitting.txt", ios_base::out);
		for (int i = 0; i < k_even_res.rows(); ++i)
			for (int j = 0; j < k_even_res.cols(); ++j) {
				int index = (*even_field_index[k])(i, j);
				auto h = field_seg->solar_scene->helios[index];
				double gt = h->flux_sum * (1 - h->sd_bk);
				double dis = gt - k_even_res(i, j);
				outFile << (*even_pos_x)(i, j) << ' ' << (*even_pos_y)(i, j) << ' ' << k_even_res(i, j) << ' ' << gt << ' ' <<  dis << ' ' << dis/gt << endl;
			}
		for (int i = 0; i < k_odd_res.rows(); ++i)
			for (int j = 0; j < k_odd_res.cols(); ++j) {
				int index = (*odd_field_index[k])(i, j);
				auto h = field_seg->solar_scene->helios[index];
				double gt = h->flux_sum * (1 - h->sd_bk);
				double dis = gt - k_odd_res(i, j);
				outFile << (*odd_pos_x)(i, j) << ' ' << (*odd_pos_y)(i, j) << ' ' << k_odd_res(i, j) << ' ' << gt << ' ' << dis << ' ' << dis/gt << endl;
			}

		vector<int>& exclude_helio_index = field_seg->solar_scene->layouts[0]->exclude_helio_index[k];
		vector<double>& exclude_helio_res = field_seg->solar_scene->layouts[0]->exclude_helio_res[0][k];
		vector<Heliostat*> helios = field_seg->solar_scene->helios;
		for (int i = 0; i < exclude_helio_res.size(); ++i) {
			Heliostat* helio = helios[exclude_helio_index[i]];
			outFile << helio->helio_pos.x() << ' ' << helio->helio_pos.z() << ' ' << exclude_helio_res[i] << ' ' <<  exclude_helio_res[i] << ' ' << 0 << ' ' << 0 << endl;
		}

		outFile.close();

	}
}

void LSPIA::checkFittingData(vector<Heliostat*>& helios, MatrixXd * field_index, vector<vector<MatrixXd*>>& fitting_data)
{
	MatrixXd* data_even = fitting_data[0][0];
	MatrixXd* data_odd = fitting_data[0][1];
	MatrixXd* data_ptr;

	int tmp_col;
	int col = field_index->cols();
	fstream outFile("fitting_error.txt", ios_base::out);
	for (int i = 0; i < field_index->rows(); i++) {
		if (i % 2) {
			tmp_col = col - 1;
			data_ptr = data_odd;
		}
		else {
			tmp_col = col;
			data_ptr = data_even;
		}
		for (int j = 0; j < tmp_col; j++) {
			Heliostat* helio = helios[(*field_index)(i, j)];
			outFile << helio->helio_pos.x() << ' ' << helio->helio_pos.z() << ' ' << helio->sd_bk <<
				' ' << (*data_ptr)(i/2, j) << ' ' << helio->sd_bk - (*data_ptr)(i/2, j) << endl;
		}
	}
	outFile.close();
}


double LSPIA::BaseFunction(const int i, const int k, const double u, const vector<double>& knot)
{
	double Nik_u = 0;
	if (k == 1) {
		if (knot[i] <= u && u <= knot[i + 1])
			Nik_u = 1;
		else
			Nik_u = 0;
	}
	else {
		double length1 = knot[i + k - 1] - knot[i];
		double length2 = knot[i + k] - knot[i + 1];
		if (length1 == 0)
			length1 = 1.0;
		if (length2 == 0.0)
			length2 = 1.0;
		Nik_u = (u - knot[i]) / length1 * BaseFunction(i, k - 1, u, knot)
			+ (knot[i + k] - u) / length2 * BaseFunction(i + 1, k - 1, u, knot);
	}
	return Nik_u;
}

vector<double> LSPIA::knot_vector(const int k, const VectorXd& param, const int N, const int M)
{
	int m = N + k;
	vector<double> knot(m + 1, 0);
	for (int i = 0; i < k + 1; i++)
		knot[i] = param[0];
	for (int i = m - k; i < m + 1; i++)
		knot[i] = param[param.rows()-1];
	for (int i = k + 1; i < m - k; i++) {
		int j = i - 3;
		double jd = j*M / (N - 3.0);
		int n = int(jd);
		double alpha = jd - n;
		knot[i] = (1 - alpha)*param[n - 1] + alpha*param[n];
	}
	return knot;
}

double LSPIA::surface_fitting(const vector<MatrixXd*>& D_even_odd, MatrixXd& P, const vector<MatrixXd*>& Nik, const double miu, const double threashold)
{

	int cnt = 0;
	double pre_error, cur_error;
	double tmp_miu = miu;
	MatrixXd preP = P;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
	pre_error = cur_error;
	//cout << "iteration: " << ++cnt << " error: " << cur_error << endl;

	preP = P;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
	//cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
	while (abs(pre_error - cur_error) > threashold) {
		pre_error = cur_error;
		cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
		//cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
		
	}
	//cout << sqrt(cur_error / (D_even_odd[0]->rows()*D_even_odd[0]->cols() + D_even_odd[1]->rows() * D_even_odd[1]->cols())) << endl;
	return cur_error;
}

double LSPIA::surface_adjusting_control_points(const vector<MatrixXd*>& D_even_odd, MatrixXd& P, 
	const vector<MatrixXd*>& Nik, const double miu)
{
	MatrixXd delta_even_uv = (*D_even_odd[0]) - (*Nik[0]) * P * Nik[2]->transpose();
	MatrixXd delta_odd_uv = (*D_even_odd[1]) - (*Nik[1]) * P * Nik[3]->transpose();
	auto error = delta_even_uv.squaredNorm() + delta_odd_uv.squaredNorm();
	
	MatrixXd delta = miu * (Nik[0]->transpose() * (delta_even_uv * (*Nik[2])) + Nik[1]->transpose() * (delta_odd_uv * (*Nik[3])));
	P += delta;
	return error;
}

vector<vector<double>> LSPIA::initParameters(const int rows, const int cols, const double start_x, const double start_y ,const double end_x, const double end_y, const vector<int>& ctrl_num)
{
	// Calculate the parameters
	VectorXd param_u(rows), param_v(cols);
	param_u.setLinSpaced(start_y, end_y);
	param_v.setLinSpaced(start_x, end_x);

	int p = 3;
	int q = 3;
	vector<vector<double>> knot_uv(2);
	knot_uv[0] = knot_vector(p, param_u, ctrl_num[0], rows);
	knot_uv[1] = knot_vector(q, param_v, ctrl_num[1], cols);
	return knot_uv;
}

vector<MatrixXd*> LSPIA::initBaseFunction(vector<vector<double>>& knot_uv, MatrixXd *even_x, MatrixXd *even_y, MatrixXd *odd_x, MatrixXd*odd_y, const vector<int>& ctrl_num)
{
	// Calculate b-spline blending basis
	int p = 3;
	int q = 3;
	int row_even = even_x->rows();
	int row_odd = odd_x->rows();
	int col = even_x->cols();

	MatrixXd* Nik_u_even = new MatrixXd(row_even, ctrl_num[0]);
	MatrixXd* Nik_u_odd = new MatrixXd(row_odd, ctrl_num[0]);
	MatrixXd* Nik_v_even = new MatrixXd(col, ctrl_num[1]);
	MatrixXd* Nik_v_odd = new MatrixXd(col, ctrl_num[1]);
	
	for (int i = row_even - 1; i > -1; --i) 
		for (int j = 0; j < ctrl_num[0]; ++j) 
			(*Nik_u_even)(i, j) = BaseFunction(j, p + 1, (*even_y)(i, 0), knot_uv[0]);

	for (int i = row_odd - 1; i > -1; --i) 
		for (int j = 0; j < ctrl_num[0]; ++j) 
			(*Nik_u_odd)(i, j) = BaseFunction(j, p + 1, (*odd_y)(i, 0), knot_uv[0]);

	for (int i = 0; i < col; ++i) 
		for (int j = 0; j < ctrl_num[1]; ++j) {
			(*Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*even_x)(0, i), knot_uv[1]);
			(*Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*odd_x)(0, i), knot_uv[1]);
		}

	vector<MatrixXd*> Nik = { Nik_u_even, Nik_u_odd, Nik_v_even, Nik_v_odd };

	return Nik;
}

MatrixXd LSPIA::initCtrlPoints(MatrixXd * even_res, MatrixXd * odd_res, const vector<int>& ctrl_num)
{
	int rows = even_res->rows();
	int cols = even_res->cols();
	MatrixXd ctrls(ctrl_num[0], ctrl_num[1]);
	for (int i = 0; i < ctrl_num[0]; ++i) {
		int f_i = int(rows*i / ctrl_num[0]);
		for (int j = 0; j < ctrl_num[1]; ++j) {
			int f_j = int(cols*j / ctrl_num[1]);
			ctrls(i, j) = (*even_res)(f_i, f_j);
		}
	}
	return ctrls;
}


MatrixXd LSPIA::LSF( MatrixXd& even_pos_x, MatrixXd& even_pos_y, MatrixXd& odd_pos_x, MatrixXd& odd_pos_y, MatrixXd& even_res, MatrixXd& odd_res) {
	MatrixXd M(6, 6), Z(6, 1);
	M.setZero();
	Z.setZero();
	accumulate(M, Z, even_pos_x, even_pos_y, even_res);
	accumulate(M, Z, odd_pos_x, odd_pos_y, odd_res);
	MatrixXd A = M.inverse()*Z;

	return A;
}


void LSPIA::accumulate(MatrixXd& M, MatrixXd& Z, MatrixXd& pos_x, MatrixXd& pos_y, MatrixXd& res) {
	int rows = pos_x.rows();
	int cols = pos_x.cols();
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			VectorXd v(6);
			v(0) = 1;
			v(1) = pos_x(i, j);
			v(2) = pos_y(i, j);
			v(3) = pos_x(i, j)*pos_x(i, j);
			v(4) = pos_x(i, j)*pos_y(i, j);
			v(5) = pos_y(i, j)*pos_y(i, j);
			
			M.row(0) += v;
			M.row(1) += pos_x(i, j)*v;
			M.row(2) += pos_y(i, j)*v;
			M.row(3) += pos_x(i, j)*pos_x(i, j)*v;
			M.row(4) += pos_x(i, j)*pos_y(i, j)*v;
			M.row(5) += pos_y(i, j)*pos_y(i, j)*v;
			
			Z.col(0) += res(i, j)*v.transpose();
		}
	}
}

MatrixXd LSPIA::LSF_res(MatrixXd& param, MatrixXd& pos_x, MatrixXd& pos_y) 
{
	int rows = pos_x.rows();
	int cols = pos_x.cols();
	MatrixXd m(1, 6);
	m.setOnes();
	MatrixXd res(rows, cols);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			m(0, 1) = pos_x(i, j);
			m(0, 2) = pos_y(i, j);
			m(0, 3) = pos_x(i, j)*pos_x(i, j);
			m(0, 4) = pos_x(i, j)*pos_y(i, j);
			m(0, 5) = pos_y(i, j)*pos_y(i, j);
			res(i, j) = (m*param)(0, 0);
		}
	}
	return res;
}

