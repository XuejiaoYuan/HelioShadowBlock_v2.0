//
// Created by Amber on 2018/4/3.
//

#include "Layout.h"

void Layout::initLayout(fstream& inFile, InputMode& input_mode, int& helio_type) {
	stringstream line_stream;
	string line, word;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> layout_bound_pos.x() >> layout_bound_pos.y() >> layout_bound_pos.z();
		else if (word == "size")
			line_stream >> layout_size.x() >> layout_size.y() >> layout_size.z();
		else if (word == "inter")
			line_stream >> helio_interval.x() >> helio_interval.y() >> helio_interval.z();
		else if (word == "n")
			line_stream >> helio_num;
		else
			line_stream >> helio_type;
	}
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
}

void Layout::setHelioLayout(vector<Heliostat*> helios)
{
	// layout_first_helio_center.y() = helios[0]->helio_pos.y();
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));

	vector<int> row_col(2, 0);
	for (auto&helio : helios) {
		set<vector<int>> pos;
		for (auto&v : helio->vertex) {
			row_col[0] = (v.z() - layout_bound_pos.z()) / helio_interval.z();	// smaller z is, smaller row is
			row_col[1] = (v.x() - layout_bound_pos.x()) / helio_interval.x();	// smaller x is, smaller col is
			pos.find(row_col);
			if (pos.count(row_col) == 0) {
				pos.insert(row_col);
				helio_layout[row_col[0]][row_col[1]].push_back(helio);
			}
		}
	}
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
}

void Layout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	// ���ܾ��ξ����Ų�
	// 1. ����Layout�ĸ�������
	helio_interval = Vector3f((*field_args[0])[0], (*field_args[0])[1], (*field_args[0])[2]);			// ���վ����
	float z_start = (*field_args[1])[0];					// ���վ���һ�о����������֮�����
	int row = int((*field_args[2])[0]);						// ��������
	int col = int((*field_args[3])[0]);						// ��������
	layout_bound_pos = Vector3f(
		-(col / 2.0) * helio_interval.x(),
		1,
		z_start - (row + 0.5)*helio_interval.z()
	);
	layout_size = Vector3f(
		helio_interval.x()*col,
		helio_interval.y(),
		helio_interval.z()*row
	);

	// 2. ����helioλ�ò���
	int cnt = 0;
	for (int i = 0; i < row; i++) {
		z_start -= helio_interval.z();
		float x_start = -(col / 2.0 - 0.5)*helio_interval.x();
		for (int j = 0; j < col; j++) {
			helio = new Heliostat((HelioType)helio_type);
			helio->helio_gap = helio_gap;
			helio->helio_matrix = helio_matrix;
			helio->helio_size = helio_size;
			helio->helio_pos = Vector3f(x_start, helio_pos.y(), z_start);
			helio->helio_index = helios.size();
			calc_param(helio, recvs);
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}

	// 3. ����helio��layout�е��Ų�����
	setHelioLayout(helios);

}

void Layout::calc_param(Heliostat * helio, const vector<Receiver*>& recvs)
{
	float dis = helio->set_focus_center_index(recvs);
	if (dis <= 1000)
		helio->mAA = (float)(0.99321 - 0.0001176 * dis + 1.97 * pow(10.0, -8.0) * dis * dis);      //d<1000
	else
		helio->mAA = exp(-0.0001106 * dis);

	Vector3f reverse_sunray_dir = (helio->helio_pos - recvs[0]->focus_center[helio->focus_center_index]).normalized();
	helio->S = helio->helio_size.x() * helio->helio_size.y();

	// TODO: set helio sigma
	helio->sigma = 0;
}


void CrossRectLayout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	helio_interval = Vector3f((*field_args[0])[0], (*field_args[0])[1], (*field_args[0])[2]);		// ���վ����
	float z_start = (*field_args[1])[0];				// ���վ���һ�о����������֮�����
	int row = int((*field_args[2])[0]);						// ��������
	int col = int((*field_args[3])[0]);						// ��������
	layout_bound_pos = Vector3f(
		-(col / 2.0) * helio_interval.x(),
		1,
		z_start - (row + 0.5)*helio_interval.z()
	);
	layout_size = Vector3f(
		helio_interval.x()*col,
		helio_interval.y(),
		helio_interval.z()*row
	);
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;

	// 2. ����helioλ�ò���
	int cnt = 0;
	for (int i = 0; i < row; i++) {
		int tmp_col = col - i % 2;
		z_start -= helio_interval.z();
		float x_start = -(tmp_col / 2.0 - 0.5)*helio_interval.x();
		for (int j = 0; j < tmp_col; j++) {
			helio = new Heliostat((HelioType)helio_type);
			helio->helio_gap = helio_gap;
			helio->helio_matrix = helio_matrix;
			helio->helio_size = helio_size;
			helio->helio_pos = Vector3f(x_start, helio_pos.y(), z_start);
			helio->helio_index = helios.size();
			calc_param(helio, recvs);
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}

	// 3. ����helio��layout�е��Ų�����
	setHelioLayout(helios);
}

void FermatLayout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs)
{
	float dsep = (*field_args[0])[0];					// ���վ���Χ�а�ȫ����
	float dm = sqrt(pow(helio_size.x(), 2) + pow(helio_size.y(), 2))+ dsep;	// ���վ��Խ��߳���
	float helio_recv_dis1 = (*field_args[1])[0];		// ��һ��ͬ��Բ�������֮��ľ���
	float helio_recv_dis2 = 2 * helio_recv_dis1;		// �ڶ���ͬ�Ļ�����ȦԲ�������֮��ľ���
	float helio_recv_dis3 = 4 * helio_recv_dis1;		// ������ͬ�Ļ�����ȦԲ�������֮��ľ���
	float helio_recv_dis4 = 8 * helio_recv_dis1;		// ���ĸ�ͬ�Ļ�����ȦԲ�������֮��ľ���
	float helio_gap1 = (*field_args[2])[0] * dm;				// ��һ��ͬ�Ļ��ж��վ��ֲ����
	float helio_gap2 = (*field_args[3])[0] * dm;				// �ڶ���ͬ�Ļ��ж��վ��ֲ����
	float helio_gap3 = (*field_args[4])[0] * dm;				// ������ͬ�Ļ��ж��վ��ֲ����
	float angle_delta1 = dm / helio_recv_dis1;			// ��һ��ͬ��Բ�����վ��ڷŽǶȼ��
	float angle_delta2 = angle_delta1 / 2;				// �ڶ���ͬ�Ļ����վ��ڷŽǶȼ��
	float angle_delta3 = angle_delta1 / 4;				// ������ͬ�Ļ����վ��ڷŽǶȼ��
	int n_row1 = int(helio_recv_dis2 - helio_recv_dis1) / helio_gap1;
	int n_row2 = int(helio_recv_dis3 - helio_recv_dis2) / helio_gap2;
	int n_row3 = int(helio_recv_dis4 - helio_recv_dis3) / helio_gap3;
	helio_interval = Vector3f(dm, dm, dm);				// ���վ���Χ�м��

	// ����layout�Ĳ���
	layout_bound_pos = Vector3f(
		-helio_recv_dis4,
		1,
		-helio_recv_dis4
	);
	layout_size = Vector3f(
		2 * helio_recv_dis4,
		helio_interval.y(),
		2 * helio_recv_dis4
	);

	// ����helioλ�ò���
	setCircleHelios(helio_recv_dis1, helio_gap1, n_row1, angle_delta1, helios, recvs);
	setCircleHelios(helio_recv_dis2, helio_gap2, n_row2, angle_delta2, helios, recvs);
	setCircleHelios(helio_recv_dis3, helio_gap3, n_row3, angle_delta3, helios, recvs);

	// ����helio��layout�е��Ų�����
	setHelioLayout(helios);

#ifdef DEBUG
	fstream outFile("fermat_helio.scn", ios_base::out);
	for (auto&h : helios) {
		outFile << h->helio_pos.x() << ' ' << h->helio_pos.z() << endl;
	}
	outFile.close();
#endif // DEBUG

}

void FermatLayout::setCircleHelios(const float R, const float gap, const int rows, const float angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	int cnt = 0;
	int h_cnt = 2 * PI / angle_delta;

	for (int i = 0; i < rows; i++) {
		float angle_d = 2 * PI / h_cnt;
		float start_angle = (i + 1) % 2 * (angle_d / 2);
		float start_r = R + i * gap;

		for (int h = 0; h < h_cnt; h++) {
			helio = new Heliostat(helio_type);
			helio->helio_gap = helio_gap;
			helio->helio_matrix = helio_matrix;
			helio->helio_size = helio_size;
			helio->helio_pos = Vector3f(
				sin(start_angle + h*angle_d) * start_r,
				helio_pos.y(),
				cos(start_angle + h*angle_d) * start_r
			);
			helio->helio_poly_pos = Vector3f(
				start_angle + h*angle_d,
				helio_pos.y(),
				start_r
			);
			helio->helio_index = helios.size();
			calc_param(helio, recvs);
			helios.push_back(helio);
			helio = nullptr;
			cnt++;
		}
	}
	cout << "cnt: " << cnt << endl;
}