#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"

class LsfFieldSegment {
public:
	LsfFieldSegment(SolarScene* _solar_scene) : solar_scene(_solar_scene){}
	void setSegmentParam(int _seg_row, int _seg_col) {
		seg_row = _seg_row;
		seg_col = _seg_col;
	}
	void initFieldSegment();
	//unordered_map<int, MatrixXd> lsf_param;
	vector<vector<Vector2d>> boundPos;
	unordered_map<int, double> sample_res;

	SolarScene *solar_scene;
	int seg_row, seg_col;						// the rows and cols of field segmentation
	int sample_row, sample_col;

};