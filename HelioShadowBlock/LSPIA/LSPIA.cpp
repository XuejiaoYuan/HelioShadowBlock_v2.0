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
	if(field_seg->seg_row != 1 ||  field_seg->seg_col != 1)
		for (int i = 0; i < field_seg->seg_row; ++i)
			for (int j = 0; j < field_seg->seg_col; ++j) {
				if (i == 0 || i == field_seg->seg_row - 1) ctrl_num_row(i, j) = min(2, field_seg->seg_row) * per_row;
				else ctrl_num_row(i, j) = min(3, field_seg->seg_row) * per_row;
				if (j == 0 || j == field_seg->seg_col - 1) ctrl_num_col(i, j) = min(2, field_seg->seg_col) * per_col;
				else ctrl_num_col(i, j) = min(3, field_seg->seg_col) * per_col;
				cout << ctrl_num_row(i, j) << ' ' << ctrl_num_col(i, j) << endl;
			}

	miu = _miu;
}

//	LSPIA_surface
//	计算所有时刻下阴影遮挡曲面的控制顶点
// :param field_data: 定日镜场x与z坐标
// :param sample_filed_data: 定日镜采样点x与z坐标
// :param sample_data: 定日镜采样点阴影遮挡结果
//vector<vector<MatrixXd*>> LSPIA::LSPIA_surface(const vector<int>&ctrl_num, const double miu)
//{
//	// Step0. Init parameters
//	int row = sample_field_data[0]->rows();
//	int col = sample_field_data[0]->cols();
//	int row_even = int(row / 2) + row % 2;
//	int row_odd = int(row / 2);
//	int p = 3;
//	int q = 3;
//	MatrixXd* sample_x = sample_field_data[0];
//	MatrixXd* sample_y = sample_field_data[1];
//	delta_even_uv = new MatrixXd(row / 2, col);
//	delta_odd_uv = new MatrixXd(row / 2, col);
//	delta = new MatrixXd(ctrl_num[0], ctrl_num[1]);
//	
//	// Step1. Calculate the parameters
//	vector<double> param_u, param_v;
//	for (int i = row - 1; i > -1; i--)
//		param_u.push_back((*sample_y)(i, 0));
//	double d = ((*sample_x)(0, col - 1) - (*sample_x)(0, 0)) / (col-1);
//	for (int i = 0; i < col; i++)
//		param_v.push_back((*sample_x)(0, 0) + i*d);
//
//	// Step2. Calculate the knot vectors
//	vector<vector<double>> knot_uv(2);
//	knot_uv[0] = knot_vector(p, param_u, ctrl_num[0], row);
//	knot_uv[1] = knot_vector(q, param_v, ctrl_num[1], col);
//
//	// Step3. Calculate b-spline blending basis
//	MatrixXd* Nik_u_even = new MatrixXd(row_even, ctrl_num[0]);
//	MatrixXd* Nik_u_odd = new MatrixXd(row_odd, ctrl_num[0]);
//	MatrixXd* Nik_v_even = new MatrixXd(col, ctrl_num[1]);
//	MatrixXd* Nik_v_odd = new MatrixXd(col, ctrl_num[1]);
//	for (int i = row - 1; i > -1; i--) {
//		for (int j = 0; j < ctrl_num[0]; j++)
//			if (i % 2) {
//				(*Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
//			}
//			else {
//				(*Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
//				//cout << (*Nik_u_even)(int(i / 2), j) << ' ';
//			}
//		//cout << endl;
//	}
//	for (int i = 0; i < col; i++) {
//		for (int j = 0; j < ctrl_num[1]; j++) {
//			(*Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*sample_x)(0, i), knot_uv[1]);
//			(*Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*sample_x)(1, i), knot_uv[1]);
//		}
//	}
//	vector<MatrixXd*> Nik = { Nik_u_even, Nik_u_odd, Nik_v_even, Nik_v_odd };
//
//	// Step4. Calculate helio position blending basis
//	int f_row = field_data[0]->rows();
//	int f_col = field_data[0]->cols();
//	int f_row_even = f_row / 2 + f_row % 2;
//	int f_row_odd = f_row / 2;
//	MatrixXd* f_Nik_u_even = new MatrixXd(f_row_even, ctrl_num[0]);
//	MatrixXd* f_Nik_u_odd = new MatrixXd(f_row_odd, ctrl_num[0]);
//	MatrixXd* f_Nik_v_even = new MatrixXd(f_col, ctrl_num[1]);
//	MatrixXd* f_Nik_v_odd = new MatrixXd(f_col - 1, ctrl_num[1]);
//	for (int i = f_row - 1; i > -1; i--)
//		for (int j = 0; j < ctrl_num[0]; j++)
//			if (i % 2)
//				(*f_Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);
//			else
//				(*f_Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);
//
//	for (int i = 0; i < f_col-1; i++)
//		for (int j = 0; j < ctrl_num[1]; j++) {
//			(*f_Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(0, i), knot_uv[1]);
//			(*f_Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(1, i), knot_uv[1]);
//		}
//	for (int j = 0; j < ctrl_num[1]; j++)
//		(*f_Nik_v_even)(f_col - 1, j) = BaseFunction(j, q + 1, (*field_data[0])(0, f_col - 1), knot_uv[1]);
//
//	vector<MatrixXd*> f_Nik = { f_Nik_u_even, f_Nik_u_odd, f_Nik_v_even, f_Nik_v_odd };
//
//	// Step5. Calculate the fitting error and control points
//	MatrixXd* ctrl_sd_bk = new MatrixXd(ctrl_num[0], ctrl_num[1]);
//	// 1. Select the initial points of shadow and block 
//	for (int i = 0; i < ctrl_num[0] - 1; i++) {
//		int f_i = row*i / double(ctrl_num[0]);
//		for (int j = 0; j < ctrl_num[1] - 1; j++) {
//			int f_j = col*j / double(ctrl_num[1]);
//			(*ctrl_sd_bk)(i, j) = (*sample_sd_bk_res[0])(f_i, f_j);;
//		}
//		(*ctrl_sd_bk)(i, ctrl_num[1] - 1) = (*sample_sd_bk_res[0])(f_i, sample_sd_bk_res[0]->cols() - 1);
//	}
//	for (int j = 0; j < ctrl_num[1] - 1; j++) {
//		int f_j = col*j / double(ctrl_num[1]);
//		(*ctrl_sd_bk)(ctrl_num[0] - 1, j) = (*sample_sd_bk_res[0])(sample_sd_bk_res[0]->rows() - 1, f_j);
//	}
//	ctrl_sd_bk->bottomRightCorner<1, 1>() = sample_sd_bk_res[0]->bottomRightCorner<1, 1>();
//	vector<vector<MatrixXd*>> calc(sample_sd_bk_res.size());
//
//#pragma omp parallel for
//	for (int cnt = 0; cnt < sample_sd_bk_res.size(); cnt++) {
//		MatrixXd tmp_ctrl_sd_bk = *ctrl_sd_bk;
//		// 1. fitting shadow&block surface
//		double sd_bk_error = surface_fitting(sample_sd_bk_res[cnt], tmp_ctrl_sd_bk, Nik, miu, 1e-4);
//
//		// 2. Calculate the fitting surface
//		MatrixXd* calc_sd_bk_even = new MatrixXd((*f_Nik[0]) * tmp_ctrl_sd_bk * f_Nik[2]->transpose());
//		MatrixXd* calc_sd_bk_odd = new MatrixXd((*f_Nik[1]) * tmp_ctrl_sd_bk * f_Nik[3]->transpose());
//
//		calc[cnt].push_back(calc_sd_bk_even);
//		calc[cnt].push_back(calc_sd_bk_odd);
//	}
//
//	delete ctrl_sd_bk;
//	delete delta_even_uv;
//	delete delta_odd_uv;
//	delete delta;
//	
//	return calc;
//}

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
		fstream outFile("odd_sample.txt", ios_base::out);
		for (int i = 0; i < odd_sample_pos_x->rows(); ++i) {
			for (int j = 0; j < odd_sample_pos_x->cols(); ++j)
				outFile << "(" << (*odd_sample_pos_x)(i, j) << ", " << (*odd_sample_pos_y)(i, j) << ")\t";
			outFile << endl;
		}
		outFile.close();
		outFile.open("even_sample.txt", ios_base::out);
		for (int i = 0; i < even_sample_pos_x->rows(); ++i) {
			for (int j = 0; j < even_sample_pos_x->cols(); ++j)
				outFile << "(" << (*even_sample_pos_x)(i, j) << ", " << (*even_sample_pos_y)(i, j) << ")\t";
			outFile << endl;
		}
		outFile.close();
		for (int i = 0; i < seg_row; ++i) {				// 同一场景下的分割区域
			for (int j = 0; j < seg_col; ++j) {
				// 1. 采样子区域划分
				Vector2i even_sample_start = even_sample_label[k][i*seg_col + j][2];
				Vector2i even_sample_len = even_sample_label[k][i*seg_col + j][3];
				Vector2i odd_sample_start = odd_sample_label[k][i*seg_col + j][2];
				Vector2i odd_sample_len = odd_sample_label[k][i*seg_col + j][3];
				cout << even_sample_start.x() << ' ' << even_sample_start.y() << ' ' << even_sample_len.x() << ' ' << even_sample_len.y() << endl;
				cout << odd_sample_start.x() << ' ' << odd_sample_start.y() << ' ' << odd_sample_len.x() << ' ' << odd_sample_len.y() << endl;
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
				Vector2i even_target_len = even_label[k][i*seg_col + j][1];
				Vector2i odd_start = odd_label[k][i*seg_col + j][2];
				Vector2i odd_len = odd_label[k][i*seg_col + j][3];
				Vector2i odd_target_start = odd_label[k][i*seg_col + j][4];
				Vector2i odd_target_len = odd_label[k][i*seg_col + j][1];

				MatrixXd sub_even_pos_x = even_pos_x->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				MatrixXd sub_even_pos_y = even_pos_y->block(even_start.x(), even_start.y(), even_len.x(), even_len.y());
				MatrixXd sub_odd_pos_x = odd_pos_x->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());
				MatrixXd sub_odd_pos_y = odd_pos_y->block(odd_start.x(), odd_start.y(), odd_len.x(), odd_len.y());
				fstream outFile("odd_pos_y.txt", ios_base::out);
				for (int i = 0; i < odd_pos_y->rows(); ++i) {
					for (int j = 0; j < odd_pos_y->cols(); ++j)
						outFile << (*odd_pos_y)(i, j) << ' ';
					outFile << endl;
				}
				outFile.close();

				// 3. 采样矩阵参数计算
				vector<int> ctrl_num = { ctrl_num_row(i,j), ctrl_num_col(i, j) };
				int even_sample_row = sub_even_sample_pos_x.rows();
				int sample_col = sub_even_sample_pos_x.cols();
				int odd_sample_row = sub_odd_sample_pos_x.rows();
				double start_x = min(sub_even_sample_pos_x(0, 0), sub_odd_sample_pos_x(0, 0));
				double start_y = min(sub_even_sample_pos_y(even_sample_row - 1, sample_col - 1), sub_odd_sample_pos_y(odd_sample_row - 1, sample_col - 1)); 
				double end_x = max(sub_even_sample_pos_x(even_sample_row - 1, sample_col - 1), sub_odd_sample_pos_x(odd_sample_row - 1, sample_col - 1));
				cout << sub_even_sample_pos_x(even_sample_row - 1, sample_col - 1) << ' ' << sub_odd_sample_pos_x(odd_sample_row - 1, sample_col - 1) << endl;
				double end_y = max(sub_even_sample_pos_y(0, 0), sub_odd_sample_pos_y(0, 0));
				cout << start_x << ' ' << start_y << ' ' << end_x << ' ' << end_y << endl;
				vector<vector<double>> knot_uv = initParameters(
					even_sample_row + odd_sample_row, sample_col, start_x, start_y, end_x, end_y, ctrl_num
				);
				
				outFile.open("knot_uv.txt", ios_base::out);
				for (int i = 0; i < knot_uv.size(); i++) {
					for (int j = 0; j < knot_uv[i].size(); ++j)
						outFile << knot_uv[i][j] << ' ';
					outFile << endl;
				}
				outFile.close();

				outFile.open("sub_odd_pos_y.txt", ios_base::out);
				for (int i = 0; i < sub_odd_pos_y.rows(); ++i) {
					for (int j = 0; j < sub_odd_pos_y.cols(); ++j)
						outFile << sub_odd_pos_y(i, j) << ' ';
					outFile << endl;
				}
				outFile.close();

				outFile.open("sub_sample_odd_pos_y.txt", ios_base::out);
				for (int i = 0; i < sub_odd_sample_pos_y.rows(); ++i) {
					for (int j = 0; j < sub_odd_sample_pos_y.cols(); ++j)
						outFile << sub_odd_sample_pos_y(i, j) << ' ';
					outFile << endl;
				}
				outFile.close();


				vector<MatrixXd*> Nik = initBaseFunction(knot_uv, &sub_even_sample_pos_x, &sub_even_sample_pos_y, &sub_odd_sample_pos_x, &sub_odd_sample_pos_y, ctrl_num);
				MatrixXd ctrls = initCtrlPoints(&sub_even_sample_res, &sub_odd_sample_res, ctrl_num);

				// 4. 结果矩阵参数计算
				vector<MatrixXd*> f_Nik = initBaseFunction(knot_uv, &sub_even_pos_x, &sub_even_pos_y, &sub_odd_pos_x, &sub_odd_pos_y, ctrl_num);
				
				for (int t = 0; t < tCnt; ++t) {					// 同一区域不同时刻结果拟合
					//MatrixXd sub_even_target = even_res[k][t]->block(even_target_start.x(), even_target_start.y(), even_end.x(), even_end.y());
					//MatrixXd sub_odd_target = odd_res[k][t]->block(odd_target_start.x(), odd_target_start.y(), odd_target_end.x(), odd_target_end.y());
					MatrixXd t_sub_even_sample_res = even_sample_res[t][k]->block(even_sample_start.x(), even_sample_start.y(), even_sample_len.x(), even_sample_len.y());
					MatrixXd t_sub_odd_sample_res = odd_sample_res[t][k]->block(odd_sample_start.x(), odd_sample_start.y(), odd_sample_len.x(), odd_sample_len.y());
					double sd_bk_error = surface_fitting(vector<MatrixXd*>{&t_sub_even_sample_res, &t_sub_odd_sample_res}, ctrls, Nik, miu, 1e-4);

					outFile.open("sub_even_sample_res.txt", ios_base::out);
					for (int i = 0; i < t_sub_even_sample_res.rows(); ++i) {
						for (int j = 0; j < t_sub_even_sample_res.cols(); ++j)
							outFile << t_sub_even_sample_res(i, j) << ' ';
						outFile << endl;
					}
					outFile.close();

					// 2. Calculate the fitting surface
					MatrixXd calc_sd_bk_even = (*f_Nik[0]) * ctrls * f_Nik[2]->transpose();
					MatrixXd calc_sd_bk_odd = (*f_Nik[1]) *  ctrls * f_Nik[3]->transpose();
					cout << even_target_start.x() << ' ' << even_target_start.y() << ' ' << even_target_len.x() << ' ' << even_target_len.y() << endl;

					MatrixXd sub_even_target = calc_sd_bk_even.block(even_target_start.x(), even_target_start.y(), even_target_len.x(), even_target_len.y());
					MatrixXd sub_odd_target = calc_sd_bk_odd.block(odd_target_start.x(), odd_target_start.y(), odd_target_len.x(), odd_target_len.y());

					cout << even_target_start.x() << ' ' << even_target_start.y() << ' ' << even_target_len.x() << ' ' << even_target_len.y() << endl;
					cout << sub_even_target.rows() << ' ' << sub_even_target.cols() << endl;

					// for test
					fstream outFile("fitting_t" + to_string(t) + "_row" + to_string(i) + "_col" + to_string(j) + ".txt", ios_base::out);
					MatrixXd sub_even_target_pos_x = sub_even_pos_x.block(even_target_start.x(), even_target_start.y(), even_target_len.x(), even_target_len.y());
					MatrixXd sub_even_target_pos_y = sub_even_pos_y.block(even_target_start.x(), even_target_start.y(), even_target_len.x(), even_target_len.y());
					MatrixXd sub_odd_target_pos_x = sub_odd_pos_x.block(odd_target_start.x(), odd_target_start.y(), odd_target_len.x(), odd_target_len.y());
					MatrixXd sub_odd_target_pos_y = sub_odd_pos_y.block(odd_target_start.x(), odd_target_start.y(), odd_target_len.x(), odd_target_len.y());
					for (int i = 0; i < sub_even_target.rows(); ++i)
						for (int j = 0; j < sub_even_target.cols(); ++j)
							outFile << sub_even_target_pos_x(i, j) << ' ' << sub_even_target_pos_y(i, j) << ' ' << sub_even_target(i, j) << endl;

					for (int i = 0; i < sub_odd_target.rows(); ++i)
						for (int j = 0; j < sub_odd_target.cols(); ++j)
							outFile << sub_odd_target_pos_x(i, j) << ' ' << sub_odd_target_pos_y(i, j) << ' ' << sub_odd_target(i, j) << endl;

					vector<int>& exclude_helio_index = field_seg->solar_scene->layouts[0]->exclude_helio_index[k];
					vector<double>& exclude_helio_res = field_seg->solar_scene->layouts[0]->exclude_helio_res[t][k];
					vector<Heliostat*> helios = field_seg->solar_scene->helios;
					for (int i = 0; i < exclude_helio_res.size(); ++i) {
						Heliostat* helio = helios[exclude_helio_index[i]];
						outFile << helio->helio_pos.x() << ' ' << helio->helio_pos.z() << ' ' << exclude_helio_res[i] << endl;

					}
						
					outFile.close();
				}
			}
		}

	}

	fstream outFile("gt_helio.txt", ios_base::out);
	for (auto& h : field_seg->solar_scene->helios) 
		outFile << h->helio_pos.x() << ' ' << h->helio_pos.z() << ' ' << h->flux_sum << endl;
	outFile.close();
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

//void LSPIA::_lspia_surface(FieldSegment* field_seg, const vector<int>& ctrl_num, const double miu)
//{
//	// Step0. Init parameters
//	int row = sample_field_data[0]->rows();
//	int col = sample_field_data[0]->cols();
//	int row_even = int(row / 2) + row % 2;
//	int row_odd = int(row / 2);
//	int p = 3;
//	int q = 3;
//	MatrixXd* sample_x = sample_field_data[0];
//	MatrixXd* sample_y = sample_field_data[1];
//	//delta_even_uv = new MatrixXd(row / 2, col);
//	//delta_odd_uv = new MatrixXd(row / 2, col);
//	//delta = new MatrixXd(ctrl_num[0], ctrl_num[1]);
//
//	// Step1. Calculate the parameters
//	vector<double> param_u, param_v;
//	for (int i = row - 1; i > -1; i--)
//		param_u.push_back((*sample_y)(i, 0));
//	double d = ((*sample_x)(0, col - 1) - (*sample_x)(0, 0)) / (col - 1);
//	for (int i = 0; i < col; i++)
//		param_v.push_back((*sample_x)(0, 0) + i*d);
//
//	// Step2. Calculate the knot vectors
//	vector<vector<double>> knot_uv(2);
//	knot_uv[0] = knot_vector(p, param_u, ctrl_num[0], row);
//	knot_uv[1] = knot_vector(q, param_v, ctrl_num[1], col);
//
//	// Step3. Calculate b-spline blending basis
//	MatrixXd* Nik_u_even = new MatrixXd(row_even, ctrl_num[0]);
//	MatrixXd* Nik_u_odd = new MatrixXd(row_odd, ctrl_num[0]);
//	MatrixXd* Nik_v_even = new MatrixXd(col, ctrl_num[1]);
//	MatrixXd* Nik_v_odd = new MatrixXd(col, ctrl_num[1]);
//	for (int i = row - 1; i > -1; i--) {
//		for (int j = 0; j < ctrl_num[0]; j++)
//			if (i % 2) {
//				(*Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
//			}
//			else {
//				(*Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
//				//cout << (*Nik_u_even)(int(i / 2), j) << ' ';
//			}
//			//cout << endl;
//	}
//	for (int i = 0; i < col; i++) {
//		for (int j = 0; j < ctrl_num[1]; j++) {
//			(*Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*sample_x)(0, i), knot_uv[1]);
//			(*Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*sample_x)(1, i), knot_uv[1]);
//		}
//	}
//	vector<MatrixXd*> Nik = { Nik_u_even, Nik_u_odd, Nik_v_even, Nik_v_odd };
//
//	// Step4. Calculate helio position blending basis
//	int f_row = field_data[0]->rows();
//	int f_col = field_data[0]->cols();
//	int f_row_even = f_row / 2 + f_row % 2;
//	int f_row_odd = f_row / 2;
//	MatrixXd* f_Nik_u_even = new MatrixXd(f_row_even, ctrl_num[0]);
//	MatrixXd* f_Nik_u_odd = new MatrixXd(f_row_odd, ctrl_num[0]);
//	MatrixXd* f_Nik_v_even = new MatrixXd(f_col, ctrl_num[1]);
//	MatrixXd* f_Nik_v_odd = new MatrixXd(f_col - 1, ctrl_num[1]);
//	for (int i = f_row - 1; i > -1; i--)
//		for (int j = 0; j < ctrl_num[0]; j++)
//			if (i % 2)
//				(*f_Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);
//			else
//				(*f_Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);
//
//	for (int i = 0; i < f_col - 1; i++)
//		for (int j = 0; j < ctrl_num[1]; j++) {
//			(*f_Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(0, i), knot_uv[1]);
//			(*f_Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(1, i), knot_uv[1]);
//		}
//	for (int j = 0; j < ctrl_num[1]; j++)
//		(*f_Nik_v_even)(f_col - 1, j) = BaseFunction(j, q + 1, (*field_data[0])(0, f_col - 1), knot_uv[1]);
//
//	vector<MatrixXd*> f_Nik = { f_Nik_u_even, f_Nik_u_odd, f_Nik_v_even, f_Nik_v_odd };
//
//	// Step5. Calculate the fitting error and control points
//	MatrixXd* ctrl_sd_bk = new MatrixXd(ctrl_num[0], ctrl_num[1]);
//	// 1. Select the initial points of shadow and block 
//	for (int i = 0; i < ctrl_num[0] - 1; i++) {
//		int f_i = row*i / double(ctrl_num[0]);
//		for (int j = 0; j < ctrl_num[1] - 1; j++) {
//			int f_j = col*j / double(ctrl_num[1]);
//			(*ctrl_sd_bk)(i, j) = (*sample_sd_bk_res[0])(f_i, f_j);;
//		}
//		(*ctrl_sd_bk)(i, ctrl_num[1] - 1) = (*sample_sd_bk_res[0])(f_i, sample_sd_bk_res[0]->cols() - 1);
//	}
//	for (int j = 0; j < ctrl_num[1] - 1; j++) {
//		int f_j = col*j / double(ctrl_num[1]);
//		(*ctrl_sd_bk)(ctrl_num[0] - 1, j) = (*sample_sd_bk_res[0])(sample_sd_bk_res[0]->rows() - 1, f_j);
//	}
//	ctrl_sd_bk->bottomRightCorner<1, 1>() = sample_sd_bk_res[0]->bottomRightCorner<1, 1>();
//	vector<vector<MatrixXd*>> calc(sample_sd_bk_res.size());
//#pragma omp parallel for
//	for (int cnt = 0; cnt < sample_sd_bk_res.size(); cnt++) {
//		MatrixXd tmp_ctrl_sd_bk = *ctrl_sd_bk;
//		// 1. fitting shadow&block surface
//		double sd_bk_error = surface_fitting(sample_sd_bk_res[cnt], tmp_ctrl_sd_bk, Nik, miu, 1e-4);
//
//		// 2. Calculate the fitting surface
//		MatrixXd* calc_sd_bk_even = new MatrixXd((*f_Nik[0]) * tmp_ctrl_sd_bk * f_Nik[2]->transpose());
//		MatrixXd* calc_sd_bk_odd = new MatrixXd((*f_Nik[1]) * tmp_ctrl_sd_bk * f_Nik[3]->transpose());
//
//		calc[cnt].push_back(calc_sd_bk_even);
//		calc[cnt].push_back(calc_sd_bk_odd);
//	}
//
//	delete ctrl_sd_bk;
//	delete delta_even_uv;
//	delete delta_odd_uv;
//	delete delta;
//
//}

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
	//vector<MatrixXd*> D_even_odd;
	//MatrixXd* D_even = new MatrixXd(D->rows() / 2 + D->rows() % 2, D->cols());
	//MatrixXd* D_odd = new MatrixXd(D->rows() / 2, D->cols());
	//for (int i = 0; i < D->rows(); i++) {
	//	if (i % 2)
	//		D_odd->row(int(i / 2)) = D->row(i);
	//	else
	//		D_even->row(int(i / 2)) = D->row(i);
	//}
	//D_even_odd.push_back(D_even);
	//D_even_odd.push_back(D_odd);

	int cnt = 0;
	double pre_error, cur_error;
	double tmp_miu = miu;
	MatrixXd preP = P;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
	pre_error = cur_error;
	cout << "iteration: " << ++cnt << " error: " << cur_error << endl;

	preP = P;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
	cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
	while (abs(pre_error - cur_error) > threashold) {
		//if (cur_error > pre_error) {
		//	P = preP;
		//	tmp_miu -= 2* miu / 5;
		//	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
		//	cout << "iteration: " << cnt << " error: " << cur_error << endl;
		//	//preP = P;
		//}
		//else {
		//	preP = P;
		//	if (!(cnt % 5))
		//		tmp_miu += miu / 5;
		//}
		pre_error = cur_error;
		cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, tmp_miu);
		cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
		
	}
	//delete D_even;
	//delete D_odd;
	cout << sqrt(cur_error / (D_even_odd[0]->rows()*D_even_odd[0]->cols() + D_even_odd[1]->rows() * D_even_odd[1]->cols())) << endl;
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

	fstream outFile("Nik_u_even.txt", ios_base::out);
	for (int i = row_even - 1; i > -1; --i) {
		for (int j = 0; j < ctrl_num[0]; ++j)
			outFile << (*Nik_u_even)(i, j) << ' ';
		outFile << endl;
	}
	outFile.close();

	outFile.open("Nik_u_odd.txt", ios_base::out);
	for (int i = row_odd - 1; i > -1; --i) {
		for (int j = 0; j < ctrl_num[0]; ++j)
			outFile << (*Nik_u_odd)(i, j) << ' ';
		outFile << endl;
	}
	outFile.close();

	outFile.open("Nik_v_even.txt", ios_base::out);
	for (int i = 0; i < col; ++i) {
		for (int j = 0; j < ctrl_num[1]; ++j) 
			outFile << (*Nik_v_even)(i, j) << ' ';
		outFile << endl;
	}
	outFile.close();


	outFile.open("Nik_v_odd.txt", ios_base::out);
	for (int i = 0; i < col; ++i) {
		for (int j = 0; j < ctrl_num[1]; ++j)
			outFile << (*Nik_v_odd)(i, j) << ' ';
		outFile << endl;
	}
	outFile.close();
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



