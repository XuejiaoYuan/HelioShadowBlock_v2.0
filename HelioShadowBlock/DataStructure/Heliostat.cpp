//
// Created by Amber on 2018/4/3.
//
#include "Heliostat.h"

bool Heliostat::initSurfaceNormal(const Vector3f &focus_center, const Vector3f &sunray_dir) {
	if (focus_center.x() == 0 && focus_center.y() == 0 && focus_center.z() == 0) {
		cerr << "Please input receiver's data first!" << endl;
		return false;
	}
	Vector3f reflectray_dir = focus_center - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();

	getMatrixs();
	// calculate root heliostat's vertex
	setLoacalVertex();
	setWorldVertex();

	initializeSubHelio(focus_center, sunray_dir);

	return true;
}

void Heliostat::changeSurfaceNormal(const Vector3f & focus_center, const Vector3f & sunray_dir)
{
	Vector3f reflectray_dir = focus_center - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();

	getMatrixs(true);
	// calculate root heliostat's vertex
	setLoacalVertex();
	setWorldVertex();

	changeSubHelio(focus_center, sunray_dir);

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

void Heliostat::getSubHelioVertex(vector<Vector3f>& subhelio_vertex)
{
	if (helio_matrix.x() == 1 && helio_matrix.y() == 1)
		subhelio_vertex = vertex;

	for (auto&sub : subhelios) {
		for (auto&v : sub->vertex)
			subhelio_vertex.push_back(v);
	}
}

void Heliostat::getMatrixs(bool init)
{
	// calcMatrix(helio_normal, helio_pos, local2worldM, world2localM);
	Vector3f u[3];	// could be shared

	u[1] = helio_normal;

	if (abs(u[1].x()) > abs(u[1].z()))
	{
		u[2] = u[1].cross(Vector3f(0.0f, 1.0f, 0.0f)).normalized();
		u[0] = u[1].cross(u[2]).normalized();
	}
	else
	{
		Vector3f tmp_u(0.0f, 1.0f, 0.0f);
		u[0] = tmp_u.cross(u[1]).normalized();
		u[2] = u[0].cross(u[1]).normalized();
	}
	for (int i = 0; i < 3; i++) {
		local2worldM(i, 0) = u[i].x();
		local2worldM(i, 1) = u[i].y();
		local2worldM(i, 2) = u[i].z();
		local2worldM(i, 3) = 0;
	}
	local2worldM(3, 0) = helio_pos.x();
	local2worldM(3, 1) = helio_pos.y();
	local2worldM(3, 2) = helio_pos.z();
	local2worldM(3, 3) = 1;

	world2localM = local2worldM.inverse();
}

void Heliostat::setLoacalVertex()
{
	float xlength = helio_size.x() / 2;
	//float ylength = helio_size.y() / 2;
	float zlength = helio_size.z() / 2;
	vertex.clear();
	vertex.push_back(Vector3f(-xlength, 0, -zlength));
	vertex.push_back(Vector3f(-xlength, 0, +zlength));
	vertex.push_back(Vector3f(+xlength, 0, +zlength));
	vertex.push_back(Vector3f(+xlength, 0, -zlength));
}

void Heliostat::setWorldVertex()
{
	Vector3f tmp[4] = { vertex[0], vertex[1],vertex[2],vertex[3] };
	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);

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

	getMatrixs(init);
	// calculate root heliostat's vertex
	setLoacalVertex();
	setWorldVertex();

}