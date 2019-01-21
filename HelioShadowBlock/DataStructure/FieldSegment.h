#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"


class FieldSegment {

public:
	FieldSegment(SolarScene* _solar_scene) : solar_scene(_solar_scene){}
	void setSegmentParam(int _seg_row, int _seg_col, int _sample_row, int _sample_col) {
		seg_row = _seg_row;
		seg_col = _seg_col;
		sample_row = _sample_row;
		sample_col = _sample_col;
	}
	void initFieldSegment();
	void updateFieldSegment();
	vector<MatrixXd*> even_field_index, odd_field_index;				// ����������ż�������
	vector<MatrixXd*> even_pos_x, even_pos_y, odd_pos_x, odd_pos_y;
	vector<vector<MatrixXd*>> even_res, odd_res;						// ����ʱ�������о��ӵļ����� t x field_num x row x col
	vector<MatrixXd*> even_sample_field_index, odd_sample_field_index;	// field_num x sample_row x sample_col
	vector<MatrixXd*> even_sample_pos_x, even_sample_pos_y, odd_sample_pos_x, odd_sample_pos_y;		// field_num x sample_row x sample_col
	vector<vector<MatrixXd*>> even_sample_res, odd_sample_res;			// ����ʱ���²�����ļ����� t x field_num x sample_row x sample_col
	vector<vector<vector<Vector2i>>> even_segment_region_label, odd_segment_region_label;				// �洢�ָ���չ�������ʼ/��������, Ŀ���������ʼ/��������  seg_row x seg_col x 4(start_i, start_j, end_i, end_j)
	vector<vector<vector<Vector2i>>> even_sample_segment_region_label, odd_sample_segment_region_label;	// �洢����������ָ�����ʼ�ͽ�������  seg_row x seg_col x 4 (start_i, start_j, end_i, end_j)
	int seg_row, seg_col;						// the rows and cols of field segmentation
	int sample_row, sample_col;
	SolarScene *solar_scene;

	void matrixSeperate(MatrixXd* field, MatrixXd* m_x, MatrixXd* m_y);
	void initFieldMatrix(MatrixXd* m_helio_index, vector<vector<vector<Vector2i>>>& segment_region_label);
	void initSample(const int sample_row, const int sample_col, MatrixXd* m_helio_index, vector<vector<Vector2i>>& segment_region_label,
		MatrixXd* sample_field_index, MatrixXd* sample_pos_x, MatrixXd* sample_pos_y, vector<vector<vector<Vector2i>>>& sample_segment_region_label);
	void clearParameters();
	void _updateFieldSegment();
	void _updatePos(MatrixXd* m_x, MatrixXd*m_y, const int k);
	void _updateSamplePos(const int k);
};

