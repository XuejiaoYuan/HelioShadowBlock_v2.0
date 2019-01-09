#include "LSPIA.h"

void LSPIA::set_datas(vector<MatrixXd*> field_data, vector<MatrixXd*> sample_field_data, vector<MatrixXd*> sample_sd_bk_res)
{
	this->field_data = field_data;
	this->sample_sd_bk_res = sample_sd_bk_res;
	this->sample_field_data = sample_field_data;
}

//	LSPIA_surface
//	计算所有时刻下阴影遮挡曲面的控制顶点
// :param field_data: 定日镜场x与z坐标
// :param sample_filed_data: 定日镜采样点x与z坐标
// :param sample_data: 定日镜采样点阴影遮挡结果
vector<vector<MatrixXd*>> LSPIA::LSPIA_surface(const vector<int>&ctrl_num, const double miu)
{
	// Step0. Init parameters
	int row = sample_field_data[0]->rows();
	int col = sample_field_data[0]->cols();
	int row_even = int(row / 2) + row % 2;
	int row_odd = int(row / 2);
	int p = 3;
	int q = 3;
	MatrixXd* sample_x = sample_field_data[0];
	MatrixXd* sample_y = sample_field_data[1];
	delta_even_uv = new MatrixXd(row / 2, col);
	delta_odd_uv = new MatrixXd(row / 2, col);
	delta = new MatrixXd(ctrl_num[0], ctrl_num[1]);
	
	// Step1. Calculate the parameters
	vector<double> param_u, param_v;
	for (int i = row - 1; i > -1; i--)
		param_u.push_back((*sample_y)(i, 0));
	double d = ((*sample_x)(0, col - 1) - (*sample_x)(0, 0)) / (col-1);
	for (int i = 0; i < col; i++)
		param_v.push_back((*sample_x)(0, 0) + i*d);

	// Step2. Calculate the knot vectors
	vector<vector<double>> knot_uv(2);
	knot_uv[0] = knot_vector(p, param_u, ctrl_num[0], row);
	knot_uv[1] = knot_vector(q, param_v, ctrl_num[1], col);

	// Step3. Calculate b-spline blending basis
	MatrixXd* Nik_u_even = new MatrixXd(row_even, ctrl_num[0]);
	MatrixXd* Nik_u_odd = new MatrixXd(row_odd, ctrl_num[0]);
	MatrixXd* Nik_v_even = new MatrixXd(col, ctrl_num[1]);
	MatrixXd* Nik_v_odd = new MatrixXd(col, ctrl_num[1]);
	for (int i = row - 1; i > -1; i--) {
		for (int j = 0; j < ctrl_num[0]; j++)
			if (i % 2) {
				(*Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
			}
			else {
				(*Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
				//cout << (*Nik_u_even)(int(i / 2), j) << ' ';
			}
		//cout << endl;
	}
	for (int i = 0; i < col; i++) {
		for (int j = 0; j < ctrl_num[1]; j++) {
			(*Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*sample_x)(0, i), knot_uv[1]);
			(*Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*sample_x)(1, i), knot_uv[1]);
		}
	}
	vector<MatrixXd*> Nik = { Nik_u_even, Nik_u_odd, Nik_v_even, Nik_v_odd };

	// Step4. Calculate helio position blending basis
	int f_row = field_data[0]->rows();
	int f_col = field_data[0]->cols();
	int f_row_even = f_row / 2 + f_row % 2;
	int f_row_odd = f_row / 2;
	MatrixXd* f_Nik_u_even = new MatrixXd(f_row_even, ctrl_num[0]);
	MatrixXd* f_Nik_u_odd = new MatrixXd(f_row_odd, ctrl_num[0]);
	MatrixXd* f_Nik_v_even = new MatrixXd(f_col, ctrl_num[1]);
	MatrixXd* f_Nik_v_odd = new MatrixXd(f_col - 1, ctrl_num[1]);
	for (int i = f_row - 1; i > -1; i--)
		for (int j = 0; j < ctrl_num[0]; j++)
			if (i % 2)
				(*f_Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);
			else
				(*f_Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, (*field_data[1])(f_row - 1 - i, 0), knot_uv[0]);

	for (int i = 0; i < f_col-1; i++)
		for (int j = 0; j < ctrl_num[1]; j++) {
			(*f_Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(0, i), knot_uv[1]);
			(*f_Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*field_data[0])(1, i), knot_uv[1]);
		}
	for (int j = 0; j < ctrl_num[1]; j++)
		(*f_Nik_v_even)(f_col - 1, j) = BaseFunction(j, q + 1, (*field_data[0])(0, f_col - 1), knot_uv[1]);

	vector<MatrixXd*> f_Nik = { f_Nik_u_even, f_Nik_u_odd, f_Nik_v_even, f_Nik_v_odd };

	// Step5. Calculate the fitting error and control points
	MatrixXd* ctrl_sd_bk = new MatrixXd(ctrl_num[0], ctrl_num[1]);
	// 1. Select the initial points of shadow and block 
	for (int i = 0; i < ctrl_num[0] - 1; i++) {
		int f_i = row*i / double(ctrl_num[0]);
		for (int j = 0; j < ctrl_num[1] - 1; j++) {
			int f_j = col*j / double(ctrl_num[1]);
			(*ctrl_sd_bk)(i, j) = (*sample_sd_bk_res[0])(f_i, f_j);;
		}
		(*ctrl_sd_bk)(i, ctrl_num[1] - 1) = (*sample_sd_bk_res[0])(f_i, sample_sd_bk_res[0]->cols() - 1);
	}
	for (int j = 0; j < ctrl_num[1] - 1; j++) {
		int f_j = col*j / double(ctrl_num[1]);
		(*ctrl_sd_bk)(ctrl_num[0] - 1, j) = (*sample_sd_bk_res[0])(sample_sd_bk_res[0]->rows() - 1, f_j);
	}
	ctrl_sd_bk->bottomRightCorner<1, 1>() = sample_sd_bk_res[0]->bottomRightCorner<1, 1>();
	vector<vector<MatrixXd*>> calc(sample_sd_bk_res.size());
#pragma omp parallel for
	for (int cnt = 0; cnt < sample_sd_bk_res.size(); cnt++) {
		MatrixXd tmp_ctrl_sd_bk = *ctrl_sd_bk;
		// 1. fitting shadow&block surface
		double sd_bk_error = surface_fitting(sample_sd_bk_res[cnt], tmp_ctrl_sd_bk, Nik, miu, 1e-4);

		// 2. Calculate the fitting surface
		MatrixXd* calc_sd_bk_even = new MatrixXd((*f_Nik[0]) * tmp_ctrl_sd_bk * f_Nik[2]->transpose());
		MatrixXd* calc_sd_bk_odd = new MatrixXd((*f_Nik[1]) * tmp_ctrl_sd_bk * f_Nik[3]->transpose());

		calc[cnt].push_back(calc_sd_bk_even);
		calc[cnt].push_back(calc_sd_bk_odd);
	}

	delete ctrl_sd_bk;
	delete delta_even_uv;
	delete delta_odd_uv;
	delete delta;
	
	return calc;
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

vector<double> LSPIA::knot_vector(const int k, const vector<double>& param, const int N, const int M)
{
	int m = N + k;
	vector<double> knot(m + 1, 0);
	for (int i = 0; i < k + 1; i++)
		knot[i] = param[0];
	for (int i = m - k; i < m + 1; i++)
		knot[i] = param[param.size() - 1];
	for (int i = k + 1; i < m - k; i++) {
		int j = i - 3;
		double jd = j*M / (N - 3.0);
		int n = int(jd);
		double alpha = jd - n;
		knot[i] = (1 - alpha)*param[n - 1] + alpha*param[n];
	}
	return knot;
}

double LSPIA::surface_fitting(const MatrixXd* D, MatrixXd& P, const vector<MatrixXd*>& Nik, const double miu, const double threashold)
{
	vector<MatrixXd*> D_even_odd;
	MatrixXd* D_even = new MatrixXd(D->rows() / 2 + D->rows() % 2, D->cols());
	MatrixXd* D_odd = new MatrixXd(D->rows() / 2, D->cols());
	for (int i = 0; i < D->rows(); i++) {
		if (i % 2)
			D_odd->row(int(i / 2)) = D->row(i);
		else
			D_even->row(int(i / 2)) = D->row(i);
	}
	D_even_odd.push_back(D_even);
	D_even_odd.push_back(D_odd);

	double pre_error, cur_error;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
	pre_error = cur_error;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
	int cnt = 0;
	while (abs(pre_error - cur_error) > threashold) {
		pre_error = cur_error;
		cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
		cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
	}
	delete D_even;
	delete D_odd;
	return cur_error;
}

double LSPIA::surface_adjusting_control_points(const vector<MatrixXd*>& D_even_odd, MatrixXd& P, const vector<MatrixXd*>& Nik, const double miu)
{
	*delta_even_uv = (*D_even_odd[0]) - (*Nik[0]) * P * Nik[2]->transpose();
	*delta_odd_uv = (*D_even_odd[1]) - (*Nik[1]) * P * Nik[3]->transpose();
	*delta = miu * (Nik[0]->transpose() * (*delta_even_uv * (*Nik[2])) + Nik[1]->transpose() * (*delta_odd_uv * (*Nik[3])));
	P += *delta;
	auto error = delta_even_uv->squaredNorm() + delta_odd_uv->squaredNorm();
	return error;
}



