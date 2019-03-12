#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../DataStructure/Clip/Clipper/clipper.hpp"
#include "../GaussLegendre/GaussLegendre.h"
#include "../DataStructure/FieldSegment.h"
using namespace ClipperLib;


typedef enum {
	RectFieldType, CrossRectFieldType, FermatFieldType, RadialFieldType
}FieldType;


class SdBkCalc
{
public:
	SdBkCalc(const FieldType& _field_type, SolarScene* _solar_scene, GaussLegendre* _gl) {
		this->field_type = _field_type;
		this->solar_scene = _solar_scene;
		this->gl = _gl;
	}
	double calcSingleShadowBlock(int helio_index);
	double calcSingleFluxSum(int helio_index, const double DNI);
	double calcShadowBlock(const double DNI);
	vector<MatrixXd*> calcSampleShadowBlock(vector<MatrixXd*>& sample_index, const double DNI);
	void calcExcludeShadowBlock(const double DNI);
	void calcSampleShadowBlock(int sample_row, int sample_col, const double DNI);
	void saveCalcRes(const string s);

	FieldType field_type;
	SolarScene* solar_scene;
	GaussLegendre* gl;


protected:

	double helioClipper(Heliostat*helio, const Vector3d&dir, const set<vector<int>>& estimate_grid);
	double helioClipper(Heliostat*helio, const vector<Vector3d>&dir, const vector<set<vector<int>>>& estimate_grid);
	double calcAccurateIntersection(Heliostat* helio, const Vector3d&dir, set<vector<int>>&relative_grid_label);		// used ray tracing calculate accurate relative grids
	double calcAccurateIntersection(Heliostat* helio, const vector<Vector3d>& dir, vector<set<vector<int>>>& relative_helio_label);
	void calcIntersection3DDDA(Heliostat* helio, const Vector3d&dir, set<vector<int>> & relative_grid_label);			// using 3DDDA for relative grid's prediction
	double checkForRelativeHelio(const set<vector<int>>& accurate_grid, const set<vector<int>>& estimate_grid);
	double calcIntersectionPoint(const Vector3d&orig, const Vector3d&dir, const Vector3d&A, const Vector3d&B, const Vector3d&C);
	double calcFluxMap(Heliostat*helio, const double DNI);
	double _calc_flux_sum(vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat* helio, const double cos_phi, const double DNI);
	double _calc_flux_sum(vector<Vector2d>& proj_v, Heliostat* helio, const double cos_phi, const double DNI);
	double _multi_inte_flux_sum(vector<Vector2d>& proj_v, int n, Heliostat* helio, const double cos_phi, const double DNI);
	double ray_tracing_flux_sum(vector<Vector3d>& recv_v, Vector3d& recv_pos, Vector3d& recv_normal, Heliostat* helio, const Vector3d& dir, const double DNI);
	double inte_infinite_flux_sum(Heliostat* helio, const Vector3d& recv_pos, const double cos_phi, const double DNI);
	double _helio_calc(int index, int DNI);

	void flux_sum_matrix_grid(vector<Vector3d>& _recv_v, vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat* helio, const double cos_phi, const double DNI);
	void flux_sum_matrix_inte(Vector3d& recv_normal, Vector3d& fc, vector<Vector3d>& _recv_v, Matrix4d& local2world, vector<Vector2d>& proj_v, Heliostat * helio, const double cos_phi, const double DNI);

};

class RectSdBkCalc :public SdBkCalc {
public:
	RectSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl) : SdBkCalc(RectFieldType, _solar_scene, _gl) {}
	void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};
private:
	void get_row_col(const int index, int& r, int& c);
	MatrixXd* field_data_pre() { return nullptr; }
	//MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num);

};

class CrossRectSdBkCalc :public SdBkCalc {
public:
	CrossRectSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl) :SdBkCalc(CrossRectFieldType, _solar_scene, _gl) {}
	//void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};
	void save_clipper_res(const string save_path, int month, int day, int hour, int minute);

private:
	void get_row_col(const int index, int& r, int& c);
	//MatrixXd* field_data_pre();
	//MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num);
};

class FermatSdBkCalc :public SdBkCalc {
public:
	FermatSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl):SdBkCalc(FermatFieldType, _solar_scene, _gl){}
	void save_clipper_res(const string save_path, int month, int day, int hour, int minute);
	//void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};

private:
	void get_row_col(const int index, int& r, int& c) {}
	MatrixXd* field_data_pre() { return nullptr; }
	MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num) { return nullptr; }

};

class RadialFieldCalc :public SdBkCalc {
public:
	RadialFieldCalc(SolarScene* _solar_scene, GaussLegendre* _gl):SdBkCalc(RadialFieldType, _solar_scene, _gl){}
	//void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};

private:
	void get_row_col(const int index, int& r, int& c) {}
	MatrixXd* field_data_pre() { return nullptr; }
	MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num) { return nullptr; }
};

class SdBkCalcCreator {
public:
	SdBkCalc* getSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl=NULL) {
		switch (_solar_scene->layouts[0]->layout_type)
		{
		case RectLayoutType:
			return new RectSdBkCalc(_solar_scene, _gl);
		case CrossRectLayoutType:
			return new CrossRectSdBkCalc(_solar_scene, _gl);
		case FermatLayoutType:
			return new FermatSdBkCalc(_solar_scene, _gl);
		case RadialLayoutType:
			return new RadialFieldCalc(_solar_scene, _gl);
		default:
			return nullptr;
		}
	}
};