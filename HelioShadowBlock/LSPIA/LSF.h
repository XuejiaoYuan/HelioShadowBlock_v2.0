#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"
#include "../DataStructure/LsfFieldSegment.h"

class LSF {
public:
	LSF(int _seg_row, int _seg_col, int _sample_row, int _sample_col)
		: seg_row(_seg_row), seg_col(_seg_col), sample_row(_sample_row), sample_col(_sample_col) {};
	void LSF_surface(SolarScene* solar_scene);

private:
	//LsfFieldSegment field_seg;
	double calcFittingData(vector<Heliostat*>& helios,
					double row_gap, double col_gap, double row_start, double col_start);

	vector<vector<MatrixXd>> lsf_data, Z, lsf_param;
	int seg_row, seg_col;
	int sample_row, sample_col;
};