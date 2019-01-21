#include "FieldSegment.h"

void FieldSegment::initFieldSegment()
{
	vector<MatrixXd*>& field = solar_scene->layouts[0]->helio_index_store;
	vector<MatrixXd*>& m_x = solar_scene->layouts[0]->m_helio_x;
	vector<MatrixXd*>& m_y = solar_scene->layouts[0]->m_helio_y;
	MatrixXd *e_s_f_i, *e_s_p_x, *e_s_p_y;
	MatrixXd *o_s_f_i, *o_s_p_x, *o_s_p_y;

	for (int i = 0; i < field.size(); ++i) {
		e_s_f_i = new MatrixXd(sample_row, sample_col);
		e_s_p_x = new MatrixXd(sample_row, sample_col);
		e_s_p_y = new MatrixXd(sample_row, sample_col);

		o_s_f_i = new MatrixXd(sample_row, sample_col);
		o_s_p_x = new MatrixXd(sample_row, sample_col);
		o_s_p_y = new MatrixXd(sample_row, sample_col);

		matrixSeperate(field[i], m_x[i], m_y[i]);
		initFieldMatrix(even_field_index.back(), even_segment_region_label);
		initFieldMatrix(odd_field_index.back(), odd_segment_region_label);
		initSample(sample_row, sample_col, even_field_index.back(), even_segment_region_label.back(), e_s_f_i, e_s_p_x, e_s_p_y, even_sample_segment_region_label);
		initSample(sample_row, sample_col, odd_field_index.back(), odd_segment_region_label.back(), o_s_f_i, o_s_p_x, o_s_p_y, odd_sample_segment_region_label);
		even_sample_field_index.push_back(e_s_f_i);
		even_sample_pos_x.push_back(e_s_p_x);
		even_sample_pos_y.push_back(e_s_p_y);
		odd_sample_field_index.push_back(o_s_f_i);
		odd_sample_pos_x.push_back(o_s_p_x);
		odd_sample_pos_y.push_back(o_s_p_y);
		e_s_f_i = e_s_p_x = e_s_p_y = nullptr;
		o_s_f_i = o_s_p_x = o_s_p_y = nullptr;
	}
}

void FieldSegment::updateFieldSegment()
{
	LayoutType layout_type = solar_scene->layouts[0]->layout_type;
	switch (layout_type)
	{
	case RectLayoutType:
	case CrossRectLayoutType:
		_updateFieldSegment();
		break;
	case FermatLayoutType:
		clearParameters();
		initFieldSegment();
		break;
	default:
		break;
	}
}


void FieldSegment::matrixSeperate(MatrixXd * field, MatrixXd* m_x, MatrixXd * m_y)
{
	//Step 1. 将镜场矩阵按奇偶行分开，用于后续LSPIA拟合
	int rows = field->rows();
	int cols = field->cols();
	int even_rows = rows / 2 + rows % 2;
	int odd_rows = rows / 2;
	MatrixXd *m_even_field_index = new MatrixXd(even_rows, cols);
	MatrixXd *m_odd_field_index = new MatrixXd(odd_rows, cols);
	MatrixXd *m_even_x = new MatrixXd(even_rows, cols);
	MatrixXd *m_odd_x = new MatrixXd(odd_rows, cols);
	MatrixXd *m_even_y = new MatrixXd(even_rows, cols);
	MatrixXd *m_odd_y = new MatrixXd(odd_rows, cols);

	for (int i = 0; i < rows; ++i) {
		if (i % 2) {
			m_odd_field_index->row(i/2) = field->row(i);
			m_odd_x->row(i/2) = m_x->row(i);
			m_odd_y->row(i/2) = m_y->row(i);
		}
		else {
			m_even_field_index->row(i/2) = field->row(i);
			m_even_x->row(i/2) = m_x->row(i);
			m_even_y->row(i/2) = m_y->row(i);
		}
	}
	even_field_index.push_back(m_even_field_index);
	odd_field_index.push_back(m_odd_field_index);
	even_pos_x.push_back(m_even_x);
	even_pos_y.push_back(m_even_y);
	odd_pos_x.push_back(m_odd_x);
	odd_pos_y.push_back(m_odd_y);

	//fstream outFile("even_field_index.txt", ios_base::out);
	//for (int i = 0; i < even_field_index[0]->rows(); ++i) {
	//	for (int j = 0; j < even_field_index[0]->cols(); ++j)
	//		outFile << (*even_field_index[0])(i, j) << ' ';
	//	outFile << endl;
	//}
	//outFile.close();

	//outFile.open("odd_field_index.txt", ios_base::out);
	//for (int i = 0; i < odd_field_index[0]->rows(); ++i) {
	//	for (int j = 0; j < odd_field_index[0]->cols(); ++j)
	//		outFile << (*odd_field_index[0])(i, j) << ' ';
	//	outFile << endl;
	//}
	//outFile.close();
}


void FieldSegment::initFieldMatrix(MatrixXd* m_helio_index, vector<vector<vector<Vector2i>>>& segment_region_label) 
{
	// Step2. 分割矩阵，并记录每个分割区域的起始结束坐标
	int total_row = m_helio_index->rows();			// 镜场中所有定日镜行数
	int total_col = m_helio_index->cols();			// 镜场中所有定日镜列数
	int m_row = total_row / seg_row;						// 分割区域中定日镜行数
	int m_col = total_col / seg_col;						// 分割区域中定日镜列数
	int row_cnt = total_row%seg_row;
	int col_cnt = total_col%seg_col;
	int e_start_i, e_start_j, e_end_i, e_end_j;
	int tmp_row, tmp_col;
	int start_i = 0, start_j = 0;

	vector<vector<Vector2i>> seg_list;
	//fstream outFile("field_seg.txt", ios_base::out);
	for (int i = 0; i < seg_row; ++i) {
		if (row_cnt) tmp_row = m_row + 1, row_cnt--;
		else tmp_row = m_row;
		start_j = 0;
		col_cnt = total_col%seg_col;
		for (int j = 0; j < seg_col; ++j) {
			if (col_cnt) tmp_col = m_col + 1, col_cnt--;
			else tmp_col = m_col;
			Vector2i start(max(0, start_i - tmp_row), max(0, start_j - tmp_col));
			Vector2i end(min(total_row - 1, start_i + 2 * tmp_row - 1), min(total_col - 1, start_j + 2 * tmp_col - 1));
			Vector2i center_start(start_i, start_j);
			Vector2i center_end(start_i + tmp_row - 1, start_j + tmp_col - 1);
			seg_list.push_back(vector<Vector2i>{start, end, center_start, center_end});
			//outFile << i << ' ' << j << ", " << start_i << " " << start_j << ", " << center_end.x() << ' ' << center_end.y() << endl;
			start_j += tmp_col;
		}
		start_i += tmp_row;
	}
	//outFile.close();
	segment_region_label.push_back(seg_list);
}


void FieldSegment::initSample(const int sample_row, const int sample_col, MatrixXd* m_helio_index, vector<vector<Vector2i>>& segment_region_label,
	MatrixXd* sample_field_index, MatrixXd* sample_pos_x, MatrixXd* sample_pos_y,  vector<vector<vector<Vector2i>>>& sample_segment_region_label)
{
	// Step3. 采样矩阵，并更新采样分割的结果
	//for (auto&m : sample_index_field) delete m;
	//for (auto&m : sample_pos_x) delete m;
	//for (auto&m : sample_pos_y) delete m;
	//sample_index_field.clear();
	//sample_pos_x.clear();
	//sample_pos_y.clear();

	//MatrixXd *m_index, *m_x, *m_y;
	//vector<MatrixXd*> & helio_index_stores = solar_scene->layouts[0]->helio_index_store;
	vector<Heliostat*>& helios = solar_scene->helios;

	//MatrixXd* m_helio_index = helio_index_stores[k];
	int grid_sample_row = sample_row / seg_row;
	int grid_sample_col = sample_col / seg_col;
	int grid_row_cnt = sample_row%seg_row;
	int grid_col_cnt = sample_col%seg_col;
	int tmp_sample_row, tmp_sample_col;
	int start_i = 0, start_j = 0;
	//m_index = new MatrixXd(sample_row, sample_col);
	//m_x = new MatrixXd(sample_row, sample_col);
	//m_y = new MatrixXd(sample_row, sample_col);
		
	//fstream outFile("sample_seg.txt", ios_base::out);
	vector<vector<Vector2i>> seg_list;
	for (int i = 0; i < seg_row; ++i) {
		start_j = 0;
		grid_col_cnt = sample_col%seg_col;
		if (grid_row_cnt) tmp_sample_row = grid_sample_row + 1, --grid_row_cnt;
		else tmp_sample_row = grid_sample_row;
		for (int j = 0; j < seg_col; ++j) {
			if (grid_col_cnt) tmp_sample_col = grid_sample_col + 1, --grid_col_cnt;
			else tmp_sample_col = grid_sample_col;
			int s_i = segment_region_label[i*seg_col + j][2].x();
			int s_j = segment_region_label[i*seg_col + j][2].y();
			int e_i = segment_region_label[i*seg_col + j][3].x();
			int e_j = segment_region_label[i*seg_col + j][3].y();
			int f_i, f_j;

			for (int g_i = 0; g_i < tmp_sample_row; ++g_i) {
				if (g_i == tmp_sample_row - 1) f_i = e_i;
				else f_i = int((e_i - s_i + 1) *g_i / (tmp_sample_row -1)) + s_i;
				for (int g_j = 0; g_j < tmp_sample_col; ++g_j) {
					if (g_j == tmp_sample_col - 1) f_j = e_j;
					else f_j = int((e_j - s_j + 1) * g_j / (tmp_sample_col-1)) + s_j;
					//outFile << "(" << f_i << "," << f_j << ")\t";

					int index = (*m_helio_index)(f_i, f_j);
					(*sample_field_index)(start_i + g_i, start_j + g_j) = index;
					(*sample_pos_x)(start_i + g_i, start_j + g_j) = helios[index]->helio_pos.x();
					(*sample_pos_y)(start_i + g_i, start_i + g_j) = helios[index]->helio_pos.z();
				}
				//outFile << endl;
			}
				
			Vector2i sample_start(start_i, start_j);
			Vector2i sample_end(start_i + tmp_sample_row - 1, start_j + tmp_sample_col - 1);
			seg_list.push_back(vector<Vector2i>{sample_start, sample_end});
			//outFile << "start & end: " << i << ' ' << j << ", " << start_i << " " << start_j << ", " << sample_end.x() << ' ' << sample_end.y() << "\n" << endl;
			start_j += tmp_sample_col;
		}
		start_i += tmp_sample_row;
	}
	//sample_field_index = m_index;
	//sample_pos_x = m_x;
	//sample_pos_y = m_y;
	sample_segment_region_label.push_back(seg_list);

	//outFile.close();

}

void FieldSegment::clearParameters()
{
	for (int i = 0; i < even_field_index.size(); ++i) {
		delete even_field_index[i];
		delete odd_field_index[i];
		delete even_sample_field_index[i];
		delete odd_sample_field_index[i];
		delete even_sample_pos_x[i];
		delete odd_sample_pos_x[i];
		delete even_sample_pos_y[i];
		delete odd_sample_pos_y[i];
		delete even_pos_x[i];
		delete odd_pos_x[i];
		delete even_pos_y[i];
		delete odd_pos_y[i];
	}
	even_field_index.clear();
	odd_field_index.clear();
	even_sample_field_index.clear();
	odd_sample_field_index.clear();
	even_sample_pos_x.clear();
	odd_sample_pos_x.clear();
	even_sample_pos_y.clear();
	odd_sample_pos_y.clear();
	even_segment_region_label.clear();
	odd_segment_region_label.clear();
	even_sample_segment_region_label.clear();
	odd_sample_segment_region_label.clear();
	even_pos_x.clear();
	odd_pos_x.clear();
	even_pos_y.clear();
	odd_pos_y.clear();
}

void FieldSegment::_updateFieldSegment()
{
	vector<MatrixXd*>& field = solar_scene->layouts[0]->helio_index_store;
	vector<MatrixXd*>& m_x = solar_scene->layouts[0]->m_helio_x;
	vector<MatrixXd*>& m_y = solar_scene->layouts[0]->m_helio_y;

	for (int i = 0; i < field.size(); ++i) {
		_updatePos(m_x[i], m_y[i], i);
		_updateSamplePos(i);
	}

}

void FieldSegment::_updatePos(MatrixXd* m_x, MatrixXd*m_y, const int k)
{
	int rows = m_x->rows();
	int cols = m_x->cols();
	for (int i = 0; i < rows; ++i) {
		if (i % 2) {
			odd_pos_x[k]->row(i) = m_x->row(i);
			odd_pos_y[k]->row(i) = m_y->row(i);
		}
		else {
			even_pos_x[k]->row(i) = m_x->row(i);
			even_pos_y[k]->row(i) = m_y->row(i);
		}
	}
}

void FieldSegment::_updateSamplePos(const int k)
{
	vector<Heliostat*>& helios = solar_scene->helios;

	int even_row = even_sample_field_index[k]->rows();
	int even_col = even_sample_field_index[k]->cols();
	int odd_row = odd_sample_field_index[k]->rows();
	int odd_col = odd_sample_field_index[k]->cols();

	for(int i=0; i<even_row; ++i)
		for (int j = 0; j < even_col; ++j) {
			int index = (*even_sample_field_index[k])(i, j);
			(*even_sample_pos_x[k])(i, j) = helios[index]->helio_pos.x();
			(*even_sample_pos_y[k])(i, j) = helios[index]->helio_pos.z();
		}

	for(int i=0; i<odd_row;++i)
		for (int j = 0; j < even_col; ++j) {
			int index = (*odd_sample_field_index[k])(i, j);
			(*odd_sample_pos_x[k])(i, j) = helios[index]->helio_pos.x();
			(*odd_sample_pos_y[k])(i, j) = helios[index]->helio_pos.z();
		}
}