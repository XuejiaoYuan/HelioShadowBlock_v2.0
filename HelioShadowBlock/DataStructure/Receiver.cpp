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
	focus_center.push_back(recv_pos + Vector3f(recv_normal.array() * recv_size.array()));
	recv_normal_list.push_back(recv_normal);
	vector<Vector3f> vertex;
	Matrix4f local2worldM, world2localM;
	GeometryFunc::setWorldVertex(recv_size, vertex, recv_normal_list[0], recv_pos, local2worldM, world2localM, true);
	recv_vertex.push_back(vertex);
	local2worldM_list.push_back(local2worldM);
	world2localM_list.push_back(world2localM);

	mask_rows = recv_size.x() / RECEIVER_SLICE;
	mask_cols = recv_size.y() / RECEIVER_SLICE;
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
			line_stream >> recv_size.x() >> recv_size.y() >> recv_size.z();		 // �������߳����߶ȣ�0
		else if (word == "num")
			line_stream >> recv_face_num;
		else if (word == "norm")
			line_stream >> recv_normal.x() >> recv_normal.y() >> recv_normal.z();
		else
			line_stream >> recv_face;

	}
	Matrix3f m;
	float delta_angle = 2 * PI / recv_face_num;
	float pos = recv_size.y() / 2 / tan(delta_angle / 2);
	mask_rows = recv_size.x() / RECEIVER_SLICE;
	mask_cols = recv_size.y() / RECEIVER_SLICE;
	for (int i = 0; i < recv_face_num; i++) {
		m << cos(i*delta_angle), 0, sin(i*delta_angle),
			 0, 1, 0,
			-sin(i*delta_angle), 0, cos(i*delta_angle);
		recv_normal_list.push_back(m * recv_normal);
		focus_center.push_back(recv_pos + recv_normal_list[i]* pos);
		//mask_x_list.push_back(MatrixXf(row, col));
		//mask_y_list.push_back(MatrixXf(row, col));

		vector<Vector3f> vertex;
		Matrix4f local2worldM, world2localM;
		GeometryFunc::setWorldVertex(recv_size, vertex, recv_normal_list[i], focus_center[i], local2worldM, world2localM, true);
		recv_vertex.push_back(vertex);
		local2worldM_list.push_back(local2worldM);
		world2localM_list.push_back(world2localM);
	}
}