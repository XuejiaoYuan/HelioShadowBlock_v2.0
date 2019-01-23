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
		helio_interval = Vector3d(0, 0, 0);
        helio_num = 0;
        grid_num = 0;
		layout_bound_pos = Vector3d(0, 0, 0);
        layout_size = Vector3d(0, 0, 0);
		layout_row_col = Vector2i(0, 0);
    }
	void initLayout(fstream& inFile, InputMode& input_mode, int& helio_type);				// �����ļ�ʱ����ʼ�������ֲ�
	inline void setHelioLayout(Heliostat* helio);											// ���õ������վ����������������е��Ų�
	void setHelioLayout(vector<Heliostat*>& helios);										// �������������ж��վ�������
	virtual void adjustHelioLayout(vector<Heliostat*>& helios,								// �����Ż�ʱ���������������
		const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs);	


	Vector3d helio_interval;					//Interval between heliostat's center
    int helio_num;								//The number of heliostat in the field
    LayoutType layout_type;						//Heliostat field's layout type
    int grid_num;								//The number of grid in the heliostat field's layout
	Vector3d layout_bound_pos;					// The bounding box of layout
	Vector3d layout_first_helio_center;			// The first heliostat center's position in the field
	Vector3d layout_size;						//Size of the layout, length/width/thickness
	Vector2i layout_row_col;					//The rows and cols of the layout
	vector<vector<vector<Heliostat*>>> helio_layout;				//List the index of heliostats in the field
	HelioType helio_type;			//Heliostat's type
	Vector3d helio_pos;           //The position of the heliostat's center
	Vector3d helio_size;          //Heliostat's size:length, thickness, width 
	Vector2d helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col
	vector<MatrixXd*> helio_index_store;				// Store helio index into matrix, exclude the last helio in every even row
	vector<vector<int>> exclude_helio_index;		    // Stote helio index which not included in the helio index matrix
	vector<vector<vector<double>>> exclude_helio_res;	// t x field_num x helio_num
	vector<MatrixXd*> m_helio_x, m_helio_y;
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
	void adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs);
};

class FermatLayout:public Layout{
public:
    FermatLayout():Layout(FermatLayoutType){}
	void adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs);

private:
	void setCircleHelios(const int filed_index, const double R, const double gap, const int rows, const double angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs);
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
