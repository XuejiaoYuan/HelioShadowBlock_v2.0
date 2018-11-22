//
// Created by Amber on 2018/4/3.
//
// Layout
// Define the layout of the heliostats field
//

#ifndef HELIOSHADOW_LAYOUT_H
#define HELIOSHADOW_LAYOUT_H
#pragma once

#include "../DataStructure/Heliostat.h"
#include "../DataStructure/Receiver.h"


class Layout {
public:
    Layout(const LayoutType&_layout_type){
        layout_type = _layout_type;
		helio_interval = Vector3f(0, 0, 0);
        helio_num = 0;
        grid_num = 0;
		layout_bound_pos = Vector3f(0, 0, 0);
        layout_size = Vector3f(0, 0, 0);
		layout_row_col = Vector2i(0, 0);
    }
	void initLayout(fstream& inFile, InputMode& input_mode, int& helio_type);				// 输入文件时，初始化镜场分布
	virtual void setHelioLayout(vector<Heliostat*> helios);									// 设置定日镜在网格中的排布
	virtual void adjustHelioLayout(vector<Heliostat*>& helios,								// 镜场优化时调整镜场网格参数
		const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs);	
	void calc_param(Heliostat* helio, const vector<Receiver*>& recvs);


	Vector3f helio_interval;					//Interval between heliostat's center
    int helio_num;								//The number of heliostat in the field
    LayoutType layout_type;						//Heliostat field's layout type
    int grid_num;								//The number of grid in the heliostat field's layout
	Vector3f layout_bound_pos;					// The bounding box of layout
	Vector3f layout_first_helio_center;			// The first heliostat center's position in the field
	Vector3f layout_size;						//Size of the layout, length/width/thickness
	Vector2i layout_row_col;					//The rows and cols of the layout
	vector<vector<vector<Heliostat*>>> helio_layout;				//List the index of heliostats in the field
	HelioType helio_type;			//Heliostat's type
	Vector3f helio_pos;           //The position of the heliostat's center
	Vector3f helio_size;          //Heliostat's size:length, thickness, width 
	Vector2f helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col

};

//
// 180417 Only consider rectangular field in one side
//
class RectLayout:public Layout{
public:
    RectLayout():Layout(RectLayoutType){};
};

class CrossRectLayout :public Layout {
public:
	CrossRectLayout() :Layout(CrossRectLayoutType){};
	void adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs);
};

class FermatLayout:public Layout{
public:
    FermatLayout():Layout(FermatLayoutType){}
	void adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<float>*>& field_args, const vector<Receiver*>& recvs);

private:
	void setCircleHelios(const float R, const float gap, const int rows, const float angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs);
};

class RadialLayout:public Layout{
public:
    RadialLayout():Layout(RadialLayoutType){}
};

class LayoutCreator{
public:
    Layout* getLayout(const LayoutType& layout_type){
        switch(layout_type){
            case RectLayoutType:
                return new RectLayout();
			case CrossRectLayoutType:
				return new CrossRectLayout();
            case FermatLayoutType:
                return new FermatLayout();
            case RadialLayoutType:
                return new RadialLayout();
            default:
                return nullptr;
        }
    }
};

#endif //HELIOSHADOW_LAYOUT_H
