#include "LsfFieldSegment.h"

void LsfFieldSegment::initFieldSegment()
{
	//lsf_param.clear();
	boundPos.clear();
	//vector<Heliostat*>& helios = solar_scene->helios;
	//MatrixXd M(6, 6), v(6);
	//v(0) = 1;

	//for (auto&h : helios) {
	//	double pos_x = h->helio_pos.x();
	//	double pos_y = h->helio_pos.y();
	//	v(1) = pos_x;
	//	v(2) = pos_y;
	//	v(3) = pos_x*pos_x;
	//	v(4) = pos_x*pos_y;
	//	v(5) = pos_y*pos_y;

	//	M.row(0) = v;
	//	M.row(1) = pos_x*v;
	//	M.row(2) = pos_y*v;
	//	M.row(3) = pos_x*pos_x*v;
	//	M.row(4) = pos_x*pos_y*v;
	//	M.row(5) = pos_y*pos_y*v;
	//	lsf_param[h->helio_index] = M;
	//}

	double row_gap = solar_scene->layouts[0]->layout_size.z()/seg_row;
	double col_gap = solar_scene->layouts[0]->layout_size.x()/seg_col;
	Vector2d startPos = Vector2d(solar_scene->layouts[0]->layout_bound_pos.z(), 
								 solar_scene->layouts[0]->layout_bound_pos.x());
	boundPos.resize(seg_row);
	for (int i = 0; i <= seg_row; ++i) {
		for (int j = 0; j <= seg_col; ++j) {
			boundPos[i].push_back(Vector2d(startPos.x() + i*row_gap, startPos.y() + j*col_gap));
		}
	}
}

