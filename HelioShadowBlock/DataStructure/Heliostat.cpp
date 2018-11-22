//
// Created by Amber on 2018/4/3.
//
#include "Heliostat.h"

bool Heliostat::initSurfaceNormal(const vector<Vector3f> &focus_center, const Vector3f &sunray_dir) {
	float dis_min = INT_MAX;
	for (int i = 0; i < focus_center.size(); i++) {
		float dis = (focus_center[i] - helio_pos).norm();
		if (dis < dis_min) {
			dis_min = dis;
			focus_center_index = i;
		}
	}
	Vector3f reflectray_dir = focus_center[focus_center_index] - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();

	GeometryFunc::setWorldVertex(helio_size, vertex, helio_normal, helio_pos, local2worldM, world2localM, true);

	// getMatrixs(true);
	// calculate root heliostat's vertex
	// setLoacalVertex();
	// setWorldVertex();

	// initializeSubHelio(focus_center, sunray_dir);

	return true;
}

void Heliostat::changeSurfaceNormal(const vector<Vector3f>& focus_center, const Vector3f & sunray_dir)
{
	Vector3f reflectray_dir = focus_center[focus_center_index] - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();
	cos_w = (-sunray_dir).dot(helio_normal);

	GeometryFunc::setWorldVertex(helio_size, vertex, helio_normal, helio_pos, local2worldM, world2localM, true);

	calc_flux_param(focus_center[focus_center_index]);
	//getMatrixs(true);
	//setLoacalVertex();
	//setWorldVertex();
	
	// calc_l_w_ratio(reflectray_dir, focus_center[focus_center_index]);
	// changeSubHelio(focus_center, sunray_dir);
}

void Heliostat::changeSubHelio(const Vector3f & focus_center, const Vector3f & sunray_dir)
{
	vector<Vector3f> root_dir(2);
	root_dir[0] = (vertex[1] - vertex[0]).normalized();
	root_dir[1] = (vertex[3] - vertex[0]).normalized();
	int subHelio_num = helio_matrix.x()*helio_matrix.y();
	int i = 0;
	for (auto&subHelio : subhelios) {
		if (helio_type == RectangularHelioType)
			subHelio->helio_normal = helio_normal;
		subHelio->setVertex(this, root_dir, i, focus_center, sunray_dir, true);
		i++;
	}

}

float Heliostat::calcSunHelioAngle(const Vector3f & sunray_dir)
{
	Vector3f reverse_sunray_dir = -sunray_dir;
	return reverse_sunray_dir.dot(this->helio_normal);
}

void Heliostat::initHeliostat(stringstream& line_stream, fstream& inFile, LayoutType layout_type, const Vector2f& helio_gap,
	const Vector2i& helio_matrix, const vector<Vector3f>& focus_center, const Vector3f& sunray_dir)
{
	string line;
	if (layout_type == FermatLayoutType) {
		line_stream >> helio_poly_pos.x() >> helio_poly_pos.y() >> helio_poly_pos.z();
		helio_pos.x() = cos(helio_poly_pos.x())*helio_poly_pos.z();
		helio_pos.z() = sin(helio_poly_pos.x())*helio_poly_pos.z();
		helio_pos.y() = helio_poly_pos.y();
	}
	else
		line_stream >> helio_pos.x() >> helio_pos.y() >> helio_pos.z();
	getline(inFile, line);
	line_stream.clear();
	line_stream.str(line);
	line_stream >> helio_size.x() >> helio_size.y() >> helio_size.z();
	this->helio_gap = helio_gap;
	this->helio_matrix = helio_matrix;
	if(sunray_dir != Vector3f(0,0,0))
		bool flag = initSurfaceNormal(focus_center, sunray_dir);
}

void Heliostat::getSubHelioVertex(vector<Vector3f>& subhelio_vertex)
{
	if (helio_matrix.x() == 1 && helio_matrix.y() == 1)
		subhelio_vertex = vertex;

	for (auto&sub : subhelios) {
		for (auto&v : sub->vertex)
			subhelio_vertex.push_back(v);
	}
}

//void Heliostat::getMatrixs(bool init)
//{
//	// calcMatrix(helio_normal, helio_pos, local2worldM, world2localM);
//	Vector3f u[3];	// could be shared
//
//	u[1] = helio_normal;
//
//	if (abs(u[1].x()) > abs(u[1].z()))
//	{
//		u[2] = u[1].cross(Vector3f(0.0f, 1.0f, 0.0f)).normalized();
//		u[0] = u[1].cross(u[2]).normalized();
//	}
//	else
//	{
//		Vector3f tmp_u(0.0f, 1.0f, 0.0f);
//		u[0] = tmp_u.cross(u[1]).normalized();
//		u[2] = u[0].cross(u[1]).normalized();
//	}
//	for (int i = 0; i < 3; i++) {
//		local2worldM(i, 0) = u[i].x();
//		local2worldM(i, 1) = u[i].y();
//		local2worldM(i, 2) = u[i].z();
//		local2worldM(i, 3) = 0;
//	}
//	local2worldM(3, 0) = helio_pos.x();
//	local2worldM(3, 1) = helio_pos.y();
//	local2worldM(3, 2) = helio_pos.z();
//	local2worldM(3, 3) = 1;
//
//	world2localM = local2worldM.inverse();
//}
//
//void Heliostat::setLoacalVertex()
//{
//	float xlength = helio_size.x() / 2;
//	//float ylength = helio_size.y() / 2;
//	float zlength = helio_size.z() / 2;
//	vertex.clear();
//	vertex.push_back(Vector3f(-xlength, 0, -zlength));
//	vertex.push_back(Vector3f(-xlength, 0, +zlength));
//	vertex.push_back(Vector3f(+xlength, 0, +zlength));
//	vertex.push_back(Vector3f(+xlength, 0, -zlength));
//}
//
//void Heliostat::setWorldVertex()
//{
//	Vector3f tmp[4] = { vertex[0], vertex[1],vertex[2],vertex[3] };
//	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
//	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
//	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
//	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);
//
//}

void Heliostat::calc_flux_param(const Vector3f& focus_center)
{
	float t[3];
	vector<Vector3f> inter_v(3);
	Vector3f reverse_dir = (helio_pos - focus_center).normalized();
	for (int i = 0; i < 3; i++) {
		// ¼ÆËãimage plane¶¥µã
		inter_v[i] = GeometryFunc::calcIntersection(reverse_dir, focus_center, vertex[i], -reverse_dir);
	}
	float ip_w = (vertex[1] - vertex[0]).norm();
	float ip_l = (vertex[2] - vertex[1]).norm();
	l_w_ratio = ip_l / ip_w;

	flux_param = 0.5 * S * cos_w * rou * mAA * l_w_ratio / PI / sigma/ sigma;
}

float Heliostat::set_focus_center_index(const vector<Receiver*>& recvs)
{
	int recv_index = 0;
	int ret_index = 0;
	float min_d = INT_MAX;
	for (int i = 0; i < recvs.size(); i++) {
		for (int j = 0; j < recvs[i]->focus_center.size(); j++) {
			Vector3f fc = recvs[i]->focus_center[j];
			float dis = pow(fc.x() - helio_pos.x(), 2) 
				+ pow(fc.z() - helio_pos.z(), 2);
			if (dis < min_d) {
				recv_index = i;
				min_d = dis;
				ret_index = j;
			}
		}
	}
	focus_center_index = ret_index;

	Vector3f image_plane_normal = (helio_pos - recvs[recv_index]->focus_center[ret_index]).normalized();
	for (int i = 0; i < recvs.size(); i++) {
		for (int j = 0; j < recvs[i]->focus_center.size(); j++) {
			cos_phi.push_back(recvs[i]->recv_normal_list[j].dot(image_plane_normal));
		}
	}

	return min_d;
}

void Heliostat::initializeSubHelio(const Vector3f&focus_center, const Vector3f&sunray_dir)
{
	Vector2f subhelio_gap(0, 0);
	Vector2i subhelio_matrix(1, 1);

	int row = helio_matrix.x();
	int col = helio_matrix.y();
	Vector3f subhelio_size;
	subhelio_size.x() = (helio_size.x() - (col - 1)*helio_gap.x()) / col;
	subhelio_size.y() = helio_size.y();
	subhelio_size.z() = (helio_size.z() - (row - 1)*helio_gap.y()) / row;

	vector<Vector3f> root_dir(2);
	root_dir[0] = (vertex[1] - vertex[0]).normalized();
	root_dir[1] = (vertex[3] - vertex[0]).normalized();
	SubHelio* subHelio_ptr = nullptr;
	int subHelio_num = helio_matrix.x()*helio_matrix.y();
	for (int i = 0; i < subHelio_num; i++) {
		subHelio_ptr = new SubHelio();
		subHelio_ptr->helio_type = SubHelioType;
		subHelio_ptr->helio_gap = subhelio_gap;
		subHelio_ptr->helio_matrix = subhelio_matrix;
		subHelio_ptr->helio_size = subhelio_size;
		if (helio_type == RectangularHelioType)
			subHelio_ptr->helio_normal = helio_normal;
		subHelio_ptr->setVertex(this, root_dir, i, focus_center, sunray_dir);
		subhelios.push_back(subHelio_ptr);
		subHelio_ptr = nullptr;
	}
}

void SubHelio::setVertex(const Heliostat* root_helio, const vector<Vector3f>&root_dir, const int sub_index,
	const Vector3f&focus_center, const Vector3f&sunray_dir, const bool init)
{
	helio_index = sub_index;

	int row = root_helio->helio_matrix.x();
	int col = root_helio->helio_matrix.y();

	int current_row = sub_index / col;
	int current_col = sub_index % col;

	helio_pos = root_helio->vertex[0]
		+ current_row*(helio_size.z() + root_helio->helio_gap.y())*root_dir[0] + 0.5*helio_size.z()*root_dir[0]
		+ current_col*(helio_size.x() + root_helio->helio_gap.x())*root_dir[1] + 0.5*helio_size.x()*root_dir[1];

	if (root_helio->helio_type != RectangularHelioType) {
		Vector3f reflectray_dir = focus_center - helio_pos;
		reflectray_dir = reflectray_dir.normalized();
		helio_normal = (reflectray_dir - sunray_dir).normalized();
	}

	// TODO: consider sub helios
	//GeometryFunc::setWorldVertex();

}