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
	Heliostat(const HelioType& _type){}
	Heliostat(const HelioType&_type, const Vector2d& _gap, const Vector2i& _matrix, const Vector3d& _size, 
		const int _index, const Vector3d& _pos, const Vector3d& _poly_pos = Vector3d(0,0,0)) {
		helio_type = _type;
		helio_matrix = _matrix;
		helio_gap = _gap;
		helio_pos = _pos;
		helio_size = _size;
		helio_normal = Vector3d(0, 0, 0);
		sd_bk = 0;
		rou = HELIOSTAT_REFLECTIVITY;
		sigma = 1.31;	// TODO: update
		helio_index = _index;
	};
	~Heliostat() {
		for (auto&sub : subhelios) {
			delete sub;
			sub = nullptr;
		}
		vertex.clear();
		subhelios.clear();
	}
	void initHeliostat(stringstream& line_stream, fstream& inFile, LayoutType layout_type, const Vector2d& helio_gap,
		const Vector2i& helio_matrix, const vector<Vector3d>& focus_center, const Vector3d& sunray_dir = Vector3d(0, 0, 0));
	
	void getSubHelioVertex(vector<Vector3d>& subhelio_vertex);
	void initializeSubHelio(const Vector3d&focus_center, const Vector3d&sunray_dir);
	bool initSurfaceNormal(const vector<Vector3d>& focus_center, const Vector3d& sunray_dir);   // Calculate the normal of heliostat surface
	void changeSurfaceNormal(const vector<Vector3d>& focus_center, const Vector3d& sunray_dir);
	void changeSubHelio(const Vector3d& focus_center, const Vector3d& sunray_dir);
	double calcSunHelioAngle(const Vector3d& sunray_dir);
	double set_focus_center_index(const vector<Receiver*>& recvs);
	void calcFluxParam(const vector<Receiver*>& recvs);

	vector<Vector3d> vertex;		//Heliostat's vertex
	vector<SubHelio*> subhelios;
	HelioType helio_type;			//Heliostat's type
	Vector3d helio_pos;           //The position of the heliostat's center
	Vector3d helio_poly_pos;		//The poly postion of the heliostat's center
	Vector3d helio_size;          //Heliostat's size:length, thickness, width 
	Vector2d helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col
	Vector3d helio_normal;          //Heliostat's surface normal
	unsigned int helio_index;			//Heliostat's index in the field
	unsigned int focus_center_index;	//Focus center index of reciver
	Matrix4d local2worldM;		//Heliostat's transform matrixs
	Matrix4d world2localM;
	double sd_bk;				// Heliostat's shadow and block ratio
	double mAA;					// Heliostat's atomospheric attenuation factor
	double cos_w;				// Heliostat's incidence cosine efficiency
	vector<double> cos_phi;		// Receiver cosine efficiencies
	double S;					// Heliostat's surface area
	double rou;					// Heliostat's reflectivity
	double l_w_ratio;			// Heliostat's projection length and width on image plane
	double sigma;				// Heliostat's sigma
	double flux_param;			// flux_param = 0.5 * S * cos_w * rou * mAA * l_w_ration / pi
	double flux_sum;				// flux_sum = 0.5 * S * cos_w * rou * mAA * l_w_ration / pi
	double max_rela_dis;			// 计算阴影与遮挡时最大无关距离
	double min_rela_dis;			// 计算阴影与遮挡时最小相关距离
	double approx_rela_dis;		// 由公式计算得到的阴影与遮挡时最大无关距离

protected:
	bool initialized;
	void calc_flux_param(const Vector3d& focus_center);
	void setHelioVertex();
};

//class RectangularHelio :public Heliostat {
//public:
//	RectangularHelio() :Heliostat(RectangularHelioType) {};
//	bool initSurfaceNormal(const Vector3d& focus_center, const Vector3d& sunray_dir);
//
//private:
//
//};
//
//class ParaboloidHelio :public Heliostat {
//public:
//	ParaboloidHelio() :Heliostat(ParaboloidHelioType) {};
//	bool initSurfaceNormal(const Vector3d& focus_center, const Vector3d& sunray_dir) { return true; }
//
//};

class SubHelio :public Heliostat {
public:
	SubHelio() :Heliostat(SubHelioType) {}
	void setVertex(const Heliostat* root_helio, const vector<Vector3d>& root_dir,
		const int sub_index, const Vector3d&focus_center, const Vector3d&sunray_dir, const bool init = false);

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
