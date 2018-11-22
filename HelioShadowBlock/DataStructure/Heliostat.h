//
// Created by Amber on 2018/4/3.
//

#ifndef HELIOSHADOW_HELIOSTAT_H
#define HELIOSHADOW_HELIOSTAT_H
#pragma once

// #include "../Common/utils.h"
// #include "../Common/global_function.cuh"
// #include "../Common/global_function.h"
// #include "../Common/vector_arithmetic.cuh"
#include "../DataStructure/Receiver.h"

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
		sd_bk = 0;
		rou = HELIOSTAT_REFLECTIVITY;
	};
	~Heliostat() {
		for (auto&sub : subhelios) {
			delete sub;
			sub = nullptr;
		}
		vertex.clear();
		subhelios.clear();
	}
	void initHeliostat(stringstream& line_stream, fstream& inFile, LayoutType layout_type, const Vector2f& helio_gap,
		const Vector2i& helio_matrix, const vector<Vector3f>& focus_center, const Vector3f& sunray_dir = Vector3f(0, 0, 0));
	void getSubHelioVertex(vector<Vector3f>& subhelio_vertex);
	void initializeSubHelio(const Vector3f&focus_center, const Vector3f&sunray_dir);
	bool initSurfaceNormal(const vector<Vector3f>& focus_center, const Vector3f& sunray_dir);   // Calculate the normal of heliostat surface
	void changeSurfaceNormal(const vector<Vector3f>& focus_center, const Vector3f& sunray_dir);
	void changeSubHelio(const Vector3f& focus_center, const Vector3f& sunray_dir);
	float calcSunHelioAngle(const Vector3f& sunray_dir);
	float set_focus_center_index(const vector<Receiver*>& recvs);

	vector<Vector3f> vertex;		//Heliostat's vertex
	vector<SubHelio*> subhelios;
	HelioType helio_type;			//Heliostat's type
	Vector3f helio_pos;           //The position of the heliostat's center
	Vector3f helio_poly_pos;		//The poly postion of the heliostat's center
	Vector3f helio_size;          //Heliostat's size:length, thickness, width 
	Vector2f helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col
	Vector3f helio_normal;          //Heliostat's surface normal
	unsigned int helio_index;			//Heliostat's index in the field
	unsigned int focus_center_index;	//Focus center index of reciver
	Matrix4f local2worldM;		//Heliostat's transform matrixs
	Matrix4f world2localM;
	float sd_bk;				// Heliostat's shadow and block ratio
	float mAA;					// Heliostat's atomospheric attenuation factor
	float cos_w;				// Heliostat's incidence cosine efficiency
	vector<float> cos_phi;		// Receiver cosine efficiencies
	float S;					// Heliostat's surface area
	float rou;					// Heliostat's reflectivity
	float l_w_ratio;			// Heliostat's projection length and width on image plane
	float sigma;				// Heliostat's sigma
	float flux_param;			// flux_param = 0.5 * S * cos_w * rou * mAA * l_w_ration / pi
	float flux_sum;				// flux_sum = 0.5 * S * cos_w * rou * mAA * l_w_ration / pi

protected:
	bool initialized;
	void calc_flux_param(const Vector3f& focus_center);
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
