#include "Receiver.h"

void Receiver::init_recv(fstream& inFile, InputMode& input_mode)
{
	string line, word;
	stringstream line_stream;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> recv_pos.x() >> recv_pos.y() >> recv_pos.z();
		else if (word == "size")
			line_stream >> recv_size.x() >> recv_size.y() >> recv_size.z();
		else if (word == "norm")
			line_stream >> recv_normal.x() >> recv_normal.y() >> recv_normal.z();
		else
			line_stream >> recv_face;
	}
	focus_center.push_back(recv_pos + Vector3d(recv_normal.array() * recv_size.array()));
	recv_normal_list.push_back(recv_normal);
	double half_l = recv_size.x() / 2.0;
	double half_w = recv_size.z() / 2.0;
	Vector3d down_cor(0, -1, 0);
	Vector3d cor_dir = recv_normal.cross(down_cor).normalized();
	vector<Vector3d> vertex = {
		(focus_center[0] - down_cor* half_l - cor_dir*half_w),
		(focus_center[0] + down_cor* half_l - cor_dir*half_w),
		(focus_center[0] + down_cor* half_l + cor_dir*half_w),
		(focus_center[0] - down_cor* half_l + cor_dir*half_w),
	};
	recv_vertex.push_back(vertex);

	mask_rows = recv_size.x() / RECEIVER_SLICE;
	mask_cols = recv_size.z() / RECEIVER_SLICE;
}


void PolyhedronRecv::init_recv(fstream& inFile, InputMode& input_mode)
{
	string line, word;
	stringstream line_stream;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> recv_pos.x() >> recv_pos.y() >> recv_pos.z();
		else if (word == "size")
			line_stream >> recv_size.x() >> recv_size.y() >> recv_size.z();		 // 接收器边长，高度，0
		else if (word == "num")
			line_stream >> recv_face_num;
		else if (word == "norm")
			line_stream >> recv_normal.x() >> recv_normal.y() >> recv_normal.z();
		else
			line_stream >> recv_face;

	}
	Matrix3d m;
	double delta_angle = 2 * PI / recv_face_num;
	double pos = recv_size.z() / 2 / tan(delta_angle / 2);
	mask_rows = recv_size.x() / RECEIVER_SLICE;
	mask_cols = recv_size.z() / RECEIVER_SLICE;
	Vector3d down_cor(0, -1, 0);
	double half_l = recv_size.x() / 2.0;
	double half_w = recv_size.z() / 2.0;

	for (int i = 0; i < recv_face_num; i++) {
		m << cos(i*delta_angle), 0, sin(i*delta_angle),
			 0, 1, 0,
			-sin(i*delta_angle), 0, cos(i*delta_angle);
		recv_normal_list.push_back(m * recv_normal);
		focus_center.push_back(recv_pos + recv_normal_list[i]* pos);

		Vector3d cor_dir = recv_normal_list[i].cross(down_cor).normalized();
		vector<Vector3d> vertex = {
			(focus_center[i] - down_cor* half_l - cor_dir*half_w),
			(focus_center[i] + down_cor* half_l - cor_dir*half_w),
			(focus_center[i] + down_cor* half_l + cor_dir*half_w),
			(focus_center[i] - down_cor* half_l + cor_dir*half_w),
		};
		recv_vertex.push_back(vertex);
	}
}