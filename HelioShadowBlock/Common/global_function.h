#pragma once
#ifndef GEOMETRY_FUNC_H
#define GEOMETRY_FUNC_H	
#include "CommonFunc.h"

namespace GeometryFunc
{
	inline Vector3d mulMatrix(const Vector3d&vertex, const Matrix4d& matrix, bool point = true) {
		RowVector4d v(vertex.x(), vertex.y(), vertex.z(), 1);		// vertex: 1; vector: 0 
		if (!point)
			v.w() = 0;

		RowVector4d res = v*matrix;
		// double res[4] = { 0 };
		//for (int i = 0; i<4; i++)
		//	for (int j = 0; j<4; j++)
		//		res(i) += v(j) * matrix[j][i];

		if (res(3) > Epsilon)
			res /= res(3);
		return Vector3d(res(0), res(1), res(2));

	}

	inline Vector3d calcIntersection(const Vector3d& normal, const Vector3d& origin_p, const Vector3d& v, const Vector3d& dir) {
		double div = dir.dot(normal);
		double t = (origin_p - v).dot(normal);
		Vector3d inter_v = v + dir*t / div;
		return inter_v;
	}

	inline void setLocalVertex(const double l, const double w, vector<Vector3d>& vertex) {
		double half_l = l / 2.0;
		double half_w = w / 2.0;
		vertex.clear();
		vertex.push_back(Vector3d(-half_l, 0, -half_w));
		vertex.push_back(Vector3d(-half_l, 0, +half_w));
		vertex.push_back(Vector3d(+half_l, 0, +half_w));
		vertex.push_back(Vector3d(+half_l, 0, -half_w));
	}

	inline void getHelioMatrix(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM) {
		Vector3d u[3];	// could be shared

		u[1] = normal;

		if (abs(u[1].x()) > abs(u[1].z()))
		{
			u[2] = u[1].cross(Vector3d(0.0f, 1.0f, 0.0f)).normalized();
			u[0] = u[1].cross(u[2]).normalized();
		}
		else
		{
			Vector3d tmp_u(0.0f, 1.0f, 0.0f);
			u[0] = tmp_u.cross(u[1]).normalized();
			u[2] = u[0].cross(u[1]).normalized();
		}

		for (int i = 0; i < 3; i++) {
			local2worldM(i, 0) = u[i].x();
			local2worldM(i, 1) = u[i].y();
			local2worldM(i, 2) = u[i].z();
			local2worldM(i, 3) = 0;
		}
		local2worldM(3, 0) = origin_p.x();
		local2worldM(3, 1) = origin_p.y();
		local2worldM(3, 2) = origin_p.z();
		local2worldM(3, 3) = 1;

		world2localM = local2worldM.inverse();
	}

	inline void getImgPlaneMatrixs(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM, unsigned int mode = 0) {
		Vector3d u[3];	// could be shared
		Vector3d tmp_u(0.0f, 1.0f, 0.0f);
		u[1] = normal;
			// 计算image plane坐标系变换矩阵
		u[0] = tmp_u.cross(u[1]).normalized();
		u[2] = u[0].cross(u[1]).normalized();

		for (int i = 0; i < 3; i++) {
			local2worldM(i, 0) = u[i].x();
			local2worldM(i, 1) = u[i].y();
			local2worldM(i, 2) = u[i].z();
			local2worldM(i, 3) = 0;
		}
		local2worldM(3, 0) = origin_p.x();
		local2worldM(3, 1) = origin_p.y();
		local2worldM(3, 2) = origin_p.z();
		local2worldM(3, 3) = 1;

		world2localM = local2worldM.inverse();
	}

	//inline void helio_getMatrixs(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM, bool init = false) {
	//	Vector3d u[3];	// could be shared

	//	u[1] = normal;
	//	u[0] = u[1].cross(Vector3d(0.0f, 1.0f, 0.0f)).normalized();
	//	u[2] = u[0].cross(u[1]).normalized();

	//	for (int i = 0; i < 3; i++) {
	//		local2worldM(i, 0) = u[i].x();
	//		local2worldM(i, 1) = u[i].y();
	//		local2worldM(i, 2) = u[i].z();
	//		local2worldM(i, 3) = 0;
	//	}
	//	local2worldM(3, 0) = origin_p.x();
	//	local2worldM(3, 1) = origin_p.y();
	//	local2worldM(3, 2) = origin_p.z();
	//	local2worldM(3, 3) = 1;

	//	world2localM = local2worldM.inverse();
	//}


	//inline void imgplane_getMatrixs(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM) {
	//	Vector3d u[3];	// could be shared

	//	u[1] = normal;

	//	Vector3d tmp_u(0.0f, 1.0f, 0.0f);
	//	u[0] = tmp_u.cross(u[1]).normalized();
	//	u[2] = u[0].cross(u[1]).normalized();

	//	for (int i = 0; i < 3; i++) {
	//		local2worldM(i, 0) = u[i].x();
	//		local2worldM(i, 1) = u[i].y();
	//		local2worldM(i, 2) = u[i].z();
	//		local2worldM(i, 3) = 0;
	//	}
	//	local2worldM(3, 0) = origin_p.x();
	//	local2worldM(3, 1) = origin_p.y();
	//	local2worldM(3, 2) = origin_p.z();
	//	local2worldM(3, 3) = 1;

	//	world2localM = local2worldM.inverse();
	//}

	//inline void setHelioWorldVertex(Vector3d& size, vector<Vector3d>& vertex, const Vector3d& normal, 
	//	const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM, bool init = false) {
	//	GeometryFunc::setLocalVertex(size.x(), size.z(), vertex);

	//	GeometryFunc::getMatrixs(normal, origin_p, local2worldM, world2localM, 0);

	//	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
	//	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
	//	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
	//	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);
	//}

	//inline void tmp_setWorldVertex(Vector3d& size, vector<Vector3d>& vertex, const Vector3d& normal,
	//	const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM, bool init = false) {
	//	GeometryFunc::setLocalVertex(size.x(), size.z(), vertex);

	//	GeometryFunc::matrix_test(normal, origin_p, local2worldM, world2localM);

	//	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
	//	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
	//	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
	//	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);
	//}
};

#endif // GEOMETRY_FUNC_H

