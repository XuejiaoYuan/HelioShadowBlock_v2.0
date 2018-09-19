#include "LSPIA.h"

void LSPIA::set_datas(vector<MatrixXf*> field_data, vector<MatrixXf*> sample_field_data, vector<MatrixXf*> sd_bk_res)
{
	this->field_data = field_data;
	this->sd_bk_res = sd_bk_res;
	this->sample_field_data = sample_field_data;
}

//	LSPIA_surface
//	计算所有时刻下阴影遮挡曲面的控制顶点
// :param field_data: 定日镜场x与z坐标
// :param sample_filed_data: 定日镜采样点x与z坐标
// :param sample_data: 定日镜采样点阴影遮挡结果
vector<vector<MatrixXf*>> & LSPIA::LSPIA_surface(const vector<int>&ctrl_num, const float miu)
{
	// Step0. Init parameters
	int row = sample_field_data[0]->rows();
	int col = sample_field_data[0]->cols();
	int row_even = int(row / 2) + row % 2;
	int row_odd = int(row / 2);
	int p = 3;
	int q = 3;
	MatrixXf* sample_x = sample_field_data[0];
	MatrixXf* sample_y = sample_field_data[1];

	// Step1. Calculate the parameters
	vector<float> param_u, param_v;
	for (int i = row - 1; i > -1; i--)
		param_u.push_back((*sample_y)(i, 0));
	float delta = ((*sample_x)(0, col - 1) - (*sample_x)(0, 0)) / col;
	for (int i = 0; i < col; i++)
		param_v.push_back((*sample_x)(0, 0) + i*delta);

	// Step2. Calculate the knot vectors
	vector<vector<float>> knot_uv(2);
	knot_uv[0] = knot_vector(p, param_u, ctrl_num[0], row);
	knot_uv[1] = knot_vector(q, param_v, ctrl_num[1], col);

	// Step3. Calculate b-spline blending basis
	MatrixXf* Nik_u_even = new MatrixXf(row_even, ctrl_num[0]);
	MatrixXf* Nik_u_odd = new MatrixXf(row_odd, ctrl_num[0]);
	MatrixXf* Nik_v_even = new MatrixXf(col, ctrl_num[1]);
	MatrixXf* Nik_v_odd = new MatrixXf(col, ctrl_num[1]);
	for (int i = row - 1; i > -1; i--)
		for (int j = 0; j < ctrl_num[0]; j++)
			if (i % 2)
				(*Nik_u_odd)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
			else
				(*Nik_u_even)(int(i / 2), j) = BaseFunction(j, p + 1, param_u[i], knot_uv[0]);
	for (int i = 0; i < col; i++) {
		for (int j = 0; j < ctrl_num[1]; i++) {
			(*Nik_v_even)(i, j) = BaseFunction(j, q + 1, (*sample_x)(0, i), knot_uv[1]);
			(*Nik_v_odd)(i, j) = BaseFunction(j, q + 1, (*sample_x)(1, i), knot_uv[1]);
		}
	}
	vector<MatrixXf*> Nik = { Nik_u_even, Nik_u_odd, Nik_v_even, Nik_v_odd };

	// Step4. Calculate helio position blending basis
	int f_row = field_data[0]->rows();
	int f_col = field_data[0]->cols();
	int f_row_even = f_row / 2 + f_row % 2;
	int f_row_odd = f_row / 2;
	MatrixXf* f_Nik_u_even = new MatrixXf(f_row_even, ctrl_num[0]);
	MatrixXf* f_Nik_u_odd = new MatrixXf(f_row_odd, ctrl_num[0]);
	MatrixXf* f_Nik_v_even = new MatrixXf(f_col, ctrl_num[1]);
	MatrixXf* f_Nik_v_odd = new MatrixXf(f_col - 1, ctrl_num[1]);
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
	vector<MatrixXf*> f_Nik = { f_Nik_u_even, f_Nik_u_odd, f_Nik_v_even, f_Nik_v_odd };

	// Step5. Calculate the fitting error and control points
	vector<MatrixXf*> ctrl_sd(sa
		mple_data.size(), new MatrixXf(ctrl_num[0], ctrl_num[1]));
	vector<MatrixXf*> ctrl_bk(sample_data.size(), new MatrixXf(ctrl_num[0], ctrl_num[1]));
	vector<vector<MatrixXf*>> calc(sample_data.size());
	for (int cnt = 0; cnt < sample_data.size(); cnt++) {
		// 1. Select the initial points of shadow and block 
		for (int i = 0; i < ctrl_num[0] - 1; i++) {
			int f_i = row*i / float(ctrl_num[0]);
			for (int j = 0; j < ctrl_num[1] - 1; j++) {
				int f_j = col*j / float(ctrl_num[1]);
				(*ctrl_sd[cnt]) << (*sample_data[cnt][2])(f_i, f_j);
				(*ctrl_bk[cnt]) << (*sample_data[cnt][3])(f_i, f_j);
			}
			(*ctrl_sd[cnt]) << (*sample_data[cnt][2])(f_i, sample_data[cnt][2]->cols() - 1);
			(*ctrl_bk[cnt]) << (*sample_data[cnt][3])(f_i, sample_data[cnt][3]->cols() - 1);
		}
		for (int j = 0; j < ctrl_num[1] - 1; j++) {
			int f_j = col*j / float(ctrl_num[1]);
			(*ctrl_sd[cnt]) << (*sample_data[cnt][2])(sample_data[cnt][2]->rows() - 1, f_j);
			(*ctrl_bk[cnt]) << (*sample_data[cnt][3])(sample_data[cnt][3]->rows() - 1, f_j);
		}
		(*ctrl_sd[cnt]) << sample_data[cnt][2]->bottomRightCorner<1, 1>();
		(*ctrl_bk[cnt]) << sample_data[cnt][3]->bottomRightCorner<1, 1>();

		// 1.1 fitting shadow surface
		float sd_error = surface_fitting(sample_data[cnt][0], ctrl_sd[cnt], Nik, miu, 0.3);
		cout << "sd_error: " << sd_error << endl;

		// 1.2 fitting block surface
		float bk_error = surface_fitting(sample_data[cnt][1], ctrl_bk[cnt], Nik, miu, 0.3);
		cout << "bk_error: " << bk_error << endl;

		// 2. Calculate the fitting surface	
		MatrixXf* calc_sd_even = new MatrixXf((*Nik[0]) * (*ctrl_sd[cnt]) * Nik[2]->transpose());
		MatrixXf* calc_sd_odd = new MatrixXf((*Nik[1]) * (*ctrl_sd[cnt]) * Nik[3]->transpose());

		MatrixXf* calc_bk_even = new MatrixXf((*Nik[0]) * (*ctrl_sd[cnt]) * Nik[2]->transpose());
		MatrixXf* calc_bk_odd = new MatrixXf((*Nik[1]) * (*ctrl_sd[cnt]) * Nik[3]->transpose());
		calc[cnt].push_back(calc_sd_even);
		calc[cnt].push_back(calc_sd_odd);
		calc[cnt].push_back(calc_bk_even);
		calc[cnt].push_back(calc_bk_odd);
	}
	
	return calc;
}

float LSPIA::BaseFunction(const int i, const int k, const float u, const vector<float>& knot)
{
	float Nik_u = 0;
	if (k == 1) {
		if (knot[i] <= u && u <= knot[i + 1])
			Nik_u = 1;
		else
			Nik_u = 0;
	}
	else {
		float length1 = knot[i + k - 1] - knot[i];
		float length2 = knot[i + k] - knot[i + 1];
		if (length1 == 0)
			length1 = 1.0;
		if (length2 == 0.0)
			length2 = 1.0;
		Nik_u = (u - knot[i]) / length1 * BaseFunction(i, k - 1, u, knot)
			+ (knot[i + k] - u) / length2 * BaseFunction(i + 1, k - 1, u, knot);
	}
	return Nik_u;
}

vector<float>& LSPIA::knot_vector(const int k, const vector<float>& param, const int N, const int M)
{
	int m = N + k;
	vector<float> knot(m + 1, 0);
	for (int i = 0; i < k + 1; i++)
		knot[i] = param[0];
	for (int i = m - k; i < m + 1; i++)
		knot[i] = param[param.size() - 1];
	for (int i = k + 1; i < m - k; i++) {
		int j = i - 3;
		float jd = j*M / (N - 3.0);
		int n = int(jd);
		float alpha = jd - n;
		knot[i] = (1 - alpha)*param[n - 1] + alpha*param[n];
	}
	return knot;
}

float LSPIA::surface_fitting(const MatrixXf* D, MatrixXf* P, const vector<MatrixXf*>& Nik, const float miu, const float threashold)
{
	vector<MatrixXf*> D_even_odd(2);
	for (int i = 0; i < D->rows(); i++) {
		if (i % 2)
			D_even_odd[1]->row(int(i / 2)) = D->row(i);
		else
			D_even_odd[0]->row(int(i / 2)) = D->row(i);
	}
	float pre_error, cur_error;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
	pre_error = cur_error;
	cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
	int cnt = 0;
	while (abs(pre_error - cur_error) > threashold) {
		pre_error = cur_error;
		cur_error = surface_adjusting_control_points(D_even_odd, P, Nik, miu);
		cout << "iteration: " << ++cnt << " error: " << cur_error << endl;
	}
	return cur_error;
}

float LSPIA::surface_adjusting_control_points(const vector<MatrixXf*>& D_even_odd, MatrixXf* P, const vector<MatrixXf*>& Nik, const float miu)
{
	float error = 0;
	MatrixXf* delta_even_uv = new MatrixXf((*D_even_odd[0]) -(*Nik[0]) * (*P) * (*Nik[2]).transpose());
	MatrixXf* delta_odd_uv = new MatrixXf((*D_even_odd[1]) - (*Nik[1]) * (*P) * (*Nik[3]).transpose());
	MatrixXf* delta = new MatrixXf(miu*(
		(*Nik[0]).transpose() * (*delta_even_uv) * (*Nik[2]) +
		(*Nik[1]).transpose() * (*delta_odd_uv) * (*Nik[3])
		));
	*P += *delta;
	error += delta_even_uv->squaredNorm() + delta_odd_uv->squaredNorm();
	return error;
}



