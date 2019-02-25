#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"
#include "../DataStructure/LsfFieldSegment.h"

class LSF {
public:
	LSF(int _seg_row, int _seg_col): seg_row(_seg_row), seg_col(_seg_col) {};
	void LSF_surface(SolarScene* solar_scene, unordered_map<int, double>& sdbk);
	void checkFittingData(vector<Heliostat*>& helios);

private:
	//LsfFieldSegment field_seg;
	vector<vector<MatrixXd>> lsf_data, Z, lsf_param;
	int seg_row, seg_col;
};