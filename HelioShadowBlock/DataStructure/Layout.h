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

typedef enum{
    RectLayoutType, CrossRectLayoutType, FermatLayoutType, RadialLayoutType
}LayoutType;

typedef enum {
	Initial, GroundMode, ReceiverMode, LayoutMode, HeliostatMode
}InputMode;


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
	void initLayout(fstream& inFile, InputMode& input_mode, int& helio_type);
	virtual void setHelioLayout(vector<Heliostat*> helios);
	Vector3f helio_interval;                  //Interval between heliostat's center
    int helio_num;                        //The number of heliostat in the field
    LayoutType layout_type;               //Heliostat field's layout type
    int grid_num;                         //The number of grid in the heliostat field's layout
	Vector3f layout_bound_pos;                // The bounding box of layout
	Vector3f layout_first_helio_center;		// The first heliostat center's position in the field
	Vector3f layout_size;                     //Size of the layout, length/width/thickness
	Vector2i layout_row_col;					//The rows and cols of the layout
	vector<vector<vector<Heliostat*>>> helio_layout;				//List the index of heliostats in the field
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
};

class FermatLayout:public Layout{
public:
    FermatLayout():Layout(FermatLayoutType){}
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
