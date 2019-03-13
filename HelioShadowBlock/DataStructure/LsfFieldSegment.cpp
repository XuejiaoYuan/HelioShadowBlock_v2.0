#include "LsfFieldSegment.h"

// [�������ָ�]
//		����������򻮷ֳ�С���򣬱�֤����������
//
void LsfFieldSegment::initFieldSegment()
{
	boundPos.clear();

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

