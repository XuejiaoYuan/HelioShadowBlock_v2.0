#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../DataStructure/Clip/Clipper/clipper.hpp"
#include "../GaussLegendre/GaussLegendre.h"
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
		//this->clipper_res_store.resize(2);
		//this->sample_clipper_res_store.resize(2);
	}
	double calcSingleShadowBlock(int helio_index);
	void calcShadowBlock(const double DNI);
	void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false);
	MatrixXd* calcSampleShadowBlock();
	virtual void save_clipper_res(const string save_path, int month, int day, int hour, int minute) {};

	FieldType field_type;
	SolarScene* solar_scene;
	GaussLegendre* gl;

	MatrixXd* field_index;					// 定日镜场所有定日镜序号
	MatrixXd* sample_field_index;
	vector<MatrixXd*> field_data;			// 定日镜场所有定日镜x与z坐标
	vector<MatrixXd*> sample_field_data;		// 定日镜场采样定日镜x与z坐标
	MatrixXd* sample_sd_bk_res;

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

	void flux_sum_matrix_grid(vector<Vector3d>& _recv_v, vector<Vector2d>& proj_v, const int rows, const int cols, Heliostat* helio, const double cos_phi, const double DNI);
	void flux_sum_matrix_inte(Vector3d& recv_normal, Vector3d& fc, vector<Vector3d>& _recv_v, Matrix4d& local2world, vector<Vector2d>& proj_v, Heliostat * helio, const double cos_phi, const double DNI);
	virtual void get_row_col(const int index, int& r, int& c) = 0;
	// 采样计算预处理操作
	virtual MatrixXd* field_data_pre() = 0;
	virtual MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num) = 0;

};

class RectSdBkCalc :public SdBkCalc {
public:
	RectSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl) : SdBkCalc(RectFieldType, _solar_scene, _gl) {}
	void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};
private:
	void get_row_col(const int index, int& r, int& c);
	MatrixXd* field_data_pre() { return nullptr; }
	MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num);

};

class CrossRectSdBkCalc :public SdBkCalc {
public:
	CrossRectSdBkCalc(SolarScene* _solar_scene, GaussLegendre* _gl) :SdBkCalc(CrossRectFieldType, _solar_scene, _gl) {}
	//void sample_calc_preprocess(const int sample_row_num, const int sample_col_num, bool calc_s = false, bool calc_f = false) {};
	void save_clipper_res(const string save_path, int month, int day, int hour, int minute);

private:
	void get_row_col(const int index, int& r, int& c);
	MatrixXd* field_data_pre();
	MatrixXd* sample_field_data_pre(const int sample_row_num, const int sample_col_num);
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