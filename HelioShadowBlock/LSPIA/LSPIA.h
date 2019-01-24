#pragma once
#ifndef LSPIA_H
#define LSPIA_H
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"
#include "../DataStructure/FieldSegment.h"

class LSPIA
{
public:
	LSPIA() {};
	void setPreDatas(FieldSegment *_field_seg) { field_seg = _field_seg; }
	void setPreDatas(FieldSegment *field_seg, vector<int>& ctrl_num, const double miu);
	void LSPIA_surface();
	void LSF_surface();
	void checkFittingData(vector<Heliostat*>& helios, MatrixXd* field_index, vector<vector<MatrixXd*>>& fitting_data);

private:
	FieldSegment *field_seg;
	MatrixXi ctrl_num_row, ctrl_num_col;
	MatrixXd segCnt;
	double miu;

	double BaseFunction(const int i, const int k, const double u, const vector<double>&knot);
	vector<double> knot_vector(const int k, const VectorXd&param, const int N, const int M);
	double surface_fitting(const vector<MatrixXd*>& D_even_odd, MatrixXd&P,
		const vector<MatrixXd*>&Nik, const double miu, const double threashold);
	double surface_adjusting_control_points(const vector<MatrixXd*>&D_even_odd, MatrixXd&P,
		const vector <MatrixXd*>&Nik, const double miu);
	vector<vector<double>> LSPIA::initParameters(const int rows, const int cols, const double start_x, const double start_y, const double end_x, const double end_y, const vector<int>& ctrl_num);
	vector<MatrixXd*> initBaseFunction(vector<vector<double>>& knot_uv, MatrixXd *even_x, MatrixXd *even_y, MatrixXd *odd_x, MatrixXd *odd_y, const vector<int>& ctrl_num);
	MatrixXd initCtrlPoints(MatrixXd* even_res, MatrixXd* odd_res, const vector<int>& ctrl_num);
	MatrixXd LSF(MatrixXd& even_pos_x, MatrixXd& even_pos_y, MatrixXd& odd_pos_x, MatrixXd& odd_pos_y, MatrixXd& even_res, MatrixXd& odd_res);
	void accumulate(MatrixXd& M, MatrixXd& Z, MatrixXd& pos_x, MatrixXd& pos_y, MatrixXd& res);
	MatrixXd LSF_res(MatrixXd& param, MatrixXd& pos_x, MatrixXd& pos_y);
};

#endif //LSPIA_H