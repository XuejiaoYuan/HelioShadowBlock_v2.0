#include "FieldSegment.h"


void FieldSegment::initFieldMatrix(SolarScene * _solar_scene)
{
	solar_scene = _solar_scene;
	vector<Heliostat*>& helios = solar_scene->helios;
	vector<Layout*>& layouts = solar_scene->layouts;
	int row = layouts[0]->layout_row_col.x();
	int col = layouts[0]->layout_row_col.y();
	if (!field_seg.empty()) {
		for (auto&seg : field_seg) {
			for(int i=0; i<seg[2]->rows(); i++)
				for (int j = 0; j < seg[2]->cols(); j++) {
					int index = (*seg[2])(i, j);
					(*seg[0])(i, j) = helios[index]->helio_pos.x();
					(*seg[1])(i, j) = helios[index]->helio_pos.z();
				}
		}
	}
	else {
		MatrixXd *m_x, *m_z, *m_index;
		int total_row = layouts[0]->layout_row_col.x();			// 镜场中所有定日镜行数
		int total_col = layouts[0]->layout_row_col.y();			// 镜场中所有定日镜列数
		int m_row = total_row / seg_row;						// 分割区域中定日镜行数
		int m_col = total_col / seg_col;						// 分割区域中定日镜列数
		int tmp_row, tmp_col;
		int start_index = 0;
		for (int i = 0; i < seg_row; i++) {
			if (i == seg_row - 1) 
				tmp_row = total_row - (seg_row - 1)*m_row;
			else 
				tmp_row = m_row;
			for (int j = 0; j < seg_col; j++) {
				if (j == seg_col - 1) 
					tmp_col = total_col - (seg_col - 1)*m_col;
				else 
					tmp_col = m_col;
				m_x = new MatrixXd(tmp_row, tmp_col);
				m_z = new MatrixXd(tmp_row, tmp_col);
				m_index = new MatrixXd(tmp_row, tmp_col);

				start_index = i * m_row *  total_col + j * tmp_col;
				for (int m_i = 0; m_i < tmp_row; m_i++) {
					for (int m_j = 0; m_j < tmp_col; m_j++) {
						int h_index = start_index + m_i*tmp_col + m_j;
						(*m_x)(m_i, m_j) = helios[h_index]->helio_pos.x();
						(*m_z)(m_i, m_j) = helios[h_index]->helio_pos.z();
						(*m_index)(m_i, m_j) = helios[h_index]->helio_index;
					}
				}
				field_seg.push_back({ m_x, m_z, m_index });
				m_x = nullptr;
				m_z = nullptr;
				m_index = nullptr;
			}
		}
	}

}


void CrossFieldSegment::initFieldMatrix(SolarScene * _solar_scene) {
	solar_scene = _solar_scene;
	vector<Heliostat*>& helios = solar_scene->helios;
	vector<Layout*>& layouts = solar_scene->layouts;
	int row = layouts[0]->layout_row_col.x();
	int col = layouts[0]->layout_row_col.y();
	if (!field_seg.empty()) {
		for (auto&seg : field_seg) {
			for (int i = 0; i<seg[2]->rows(); i++)
				for (int j = 0; j < seg[2]->cols(); j++) {
					int index = (*seg[2])(i, j);
					(*seg[0])(i, j) = helios[index]->helio_pos.x();
					(*seg[1])(i, j) = helios[index]->helio_pos.z();
				}
		}
	}
	else {
		MatrixXd *m_x, *m_z, *m_index;
		int total_row = layouts[0]->layout_row_col.x();			// 镜场中所有定日镜行数
		int total_col = layouts[0]->layout_row_col.y();			// 镜场中所有定日镜列数
		int m_row = total_row / seg_row;						// 分割区域中定日镜行数
		int m_col = total_col / seg_col;						// 分割区域中定日镜列数
		int tmp_row, tmp_col;
		int start_index = 0;
		

		for (int i = 0; i < seg_row; i++) {
			if (i == seg_row - 1)
				tmp_row = total_row - (seg_row - 1)*m_row;
			else
				tmp_row = m_row;
			for (int j = 0; j < seg_col; j++) {
				if (j == seg_col - 1)
					tmp_col = total_col - (seg_col - 1)*m_col;
				else
					tmp_col = m_col;
				m_x = new MatrixXd(tmp_row, tmp_col);
				m_z = new MatrixXd(tmp_row, tmp_col);
				m_index = new MatrixXd(tmp_row, tmp_col);

				start_index = i * m_row *  total_col + j * tmp_col;
				for (int m_i = 0; m_i < tmp_row; m_i++) {
					for (int m_j = 0; m_j < tmp_col; m_j++) {
						int h_index = start_index + m_i*tmp_col + m_j;
						(*m_x)(m_i, m_j) = helios[h_index]->helio_pos.x();
						(*m_z)(m_i, m_j) = helios[h_index]->helio_pos.z();
						(*m_index)(m_i, m_j) = helios[h_index]->helio_index;
					}
				}
				field_seg.push_back({ m_x, m_z, m_index });
				m_x = nullptr;
				m_z = nullptr;
				m_index = nullptr;
			}
		}
	}
}