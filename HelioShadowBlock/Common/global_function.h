#pragma once
#ifndef GEOMETRY_FUNC_H
#define GEOMETRY_FUNC_H	
#include "CommonFunc.h"

namespace GeometryFunc
{

	inline void calcMatrix(const Vector3f&aligned_normal, const Vector3f& transform, Matrix4f& local2worldM, Matrix4f&world2loacalM) {
		Vector3f u[3];	// could be shared

		u[1] = aligned_normal;

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
		local2worldM(3, 0) = transform.x();
		local2worldM(3, 1) = transform.y();
		local2worldM(3, 2) = transform.z();
		local2worldM(3, 3) = 1;

		world2loacalM = local2worldM.inverse();
		// inverse(local2worldM, world2loacalM);
	}

	inline Vector3f mulMatrix(const Vector3f&vertex, Matrix4f& matrix, bool point = true) {
		RowVector4f v(vertex.x(), vertex.y(), vertex.z(), 1);		// vertex: 1; vector: 0 
		if (!point)
			v.w() = 0;

		RowVector4f res = v*matrix;
		// float res[4] = { 0 };
		//for (int i = 0; i<4; i++)
		//	for (int j = 0; j<4; j++)
		//		res(i) += v(j) * matrix[j][i];

		if (res(3) != 0)
			res /= res(3);
		return Vector3f(res(0), res(1), res(2));

	}


};

#endif // GEOMETRY_FUNC_H

