#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"


class FieldSegment {

public:
	FieldSegment(const LayoutType& _layout_type) : layout_type(_layout_type){}
	void setSegmentParam(int _seg_row, int _seg_col, int _sample_row, int _sampe_col) {
		seg_row = _seg_row;
		seg_col = _seg_col;
		sample_row = _sample_row;
		sample_col = _sampe_col;
	}
	virtual void initFieldMatrix(SolarScene* solar_scene);

	vector<MatrixXd*> field;					// heliostats' positions and index store in matrix
	vector<vector<MatrixXd*>> field_seg;		// field segmentation under one field parameters
												// m x 3: m the number of segmentation; 3: x, z, helio index
	vector<vector<MatrixXd*>> sample_field_seg;	// sample of field segmentation
												// m x 3: m the number of segmentation; 3: x, z, sample helio index
	int seg_row, seg_col;						// the rows and cols of field segmentation
	int sample_row, sample_col;					// the rows and cols in sample field
	SolarScene *solar_scene;
	LayoutType layout_type;
};

class RectFieldSegment : public FieldSegment {
public:
	RectFieldSegment() : FieldSegment(RectLayoutType) {};
};

class CrossFieldSegment : public FieldSegment {
public:
	CrossFieldSegment() :FieldSegment(CrossRectLayoutType) {}
	void initFieldMatrix(SolarScene* solar_scene);
};

class RadialFieldSegment : public FieldSegment {
public:
	RadialFieldSegment() : FieldSegment(RadialLayoutType) {}
};

class FermatFieldSegment : public FieldSegment {
public:
	FermatFieldSegment() : FieldSegment(FermatLayoutType) {}
};

class FieldSegmentCreator {
public:
	FieldSegment* getFieldSegment(const LayoutType& layout_type) {
		switch (layout_type)
		{
		case RectLayoutType:
			return new RectFieldSegment();
		case CrossRectLayoutType:
			return new CrossFieldSegment();
		case FermatLayoutType:
			return new FermatFieldSegment();
		default:
			return nullptr;
		}
	}
};