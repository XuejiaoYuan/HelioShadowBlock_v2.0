//
// Created by Amber on 2018/4/3.
//

#ifndef HELIOSHADOW_HELIOSTAT_H
#define HELIOSHADOW_HELIOSTAT_H
#pragma once

// #include "../Common/utils.h"
// #include "../Common/global_function.cuh"
#include "../Common/global_function.h"
// #include "../Common/vector_arithmetic.cuh"

typedef enum {
	RectangularHelioType, ParaboloidHelioType, SubHelioType
}HelioType;

class SubHelio;

class Heliostat {
public:
	Heliostat(const HelioType&_helio_type) {
		helio_type = helio_type;
		helio_matrix = Vector2i(1, 1);
		helio_gap = Vector2f(0, 0);
		helio_pos = Vector3f(0, 0, 0);
		helio_size = Vector3f(0, 0, 0);
		helio_normal = Vector3f(0, 0, 0);
	};
	void getSubHelioVertex(vector<Vector3f>& subhelio_vertex);
	void initializeSubHelio(const Vector3f&focus_center, const Vector3f&sunray_dir);
	bool initSurfaceNormal(const Vector3f& focus_center, const Vector3f& sunray_dir);   // Calculate the normal of heliostat surface
	void changeSurfaceNormal(const Vector3f& focus_center, const Vector3f& sunray_dir);
	void changeSubHelio(const Vector3f& focus_center, const Vector3f& sunray_dir);
	float calcSunHelioAngle(const Vector3f& sunray_dir);
	vector<Vector3f> vertex;		//Heliostat's vertex
	vector<SubHelio*> subhelios;
	HelioType helio_type;       //Heliostat's type
	Vector3f helio_pos;           //The position of the heliostat's center
	Vector3f helio_size;          //Heliostat's size:length, thickness, width 
	Vector2f helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col
	Vector3f helio_normal;          //Heliostat's surface normal
	unsigned int helio_index;			//Heliostat's index in the field
	Matrix4f local2worldM;		//Heliostat's transform matrixs
	Matrix4f world2localM;

protected:
	bool initialized;
	void getMatrixs(bool init = false);
	void setLoacalVertex();
	void setWorldVertex();
};

//class RectangularHelio :public Heliostat {
//public:
//	RectangularHelio() :Heliostat(RectangularHelioType) {};
//	bool initSurfaceNormal(const Vector3f& focus_center, const Vector3f& sunray_dir);
//
//private:
//
//};
//
//class ParaboloidHelio :public Heliostat {
//public:
//	ParaboloidHelio() :Heliostat(ParaboloidHelioType) {};
//	bool initSurfaceNormal(const Vector3f& focus_center, const Vector3f& sunray_dir) { return true; }
//
//};

class SubHelio :public Heliostat {
public:
	SubHelio() :Heliostat(SubHelioType) {}
	void setVertex(const Heliostat* root_helio, const vector<Vector3f>& root_dir,
		const int sub_index, const Vector3f&focus_center, const Vector3f&sunray_dir, const bool init = false);

private:
};

//class HeliostatCreator {
//public:
//	Heliostat* getHeliostat(const HelioType&helio_type) {
//		switch (helio_type) {
//		case RectangularHelioType:
//			return new RectangularHelio();
//		case ParaboloidHelioType:
//			return new ParaboloidHelio();
//		case SubHelioType:
//			return new SubHelio();
//		default:
//			return nullptr;
//		}
//	}
//};

#endif //HELIOSHADOW_HELIOSTAT_H
