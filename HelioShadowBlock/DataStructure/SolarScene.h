//
// Created by Amber on 2018/4/3.
//
#pragma once

#ifndef HELIOSHADOW_SOLARSCENE_H
#define HELIOSHADOW_SOLARSCENE_H


#include "Layout.h"
#include "Receiver.h"



class SolarScene {
public:
    SolarScene(){
        scene_length = 0;
        scene_width = 0;
        layouts.clear();
        helios.clear();
        recvs.clear();
    }
	bool initSolarScene(const string&scene_filepath, const Vector3d&sunray_dir);
	bool changeSolarScene(const Vector3d&sunray_dir);		//change the surface normal of heiliostats

	// function for parameter optimisim
	bool initFieldParam(const string&file_name);
	bool adjustFieldParam(const vector<vector<double>*>& field_args);


    double scene_length;           //heliostat ground's length and width
    double scene_width;

    vector<Layout*>  layouts;
    vector<Heliostat*> helios;
	// vector<vector<Heliostat*>> helios;
	vector<Receiver*> recvs;
	Vector3d sunray_dir;
};

#endif //HELIOSHADOW_SOLARSCENE_H
