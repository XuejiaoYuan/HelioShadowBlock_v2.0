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

    bool initSolarScene(const string&scene_filepath, const Vector3f&sunray_dir);
	bool changeSolarScene(const Vector3f&sunray_dir);

    float scene_length;           //heliostat ground's length and width
    float scene_width;

    vector<Layout*>  layouts;
    vector<Heliostat*> helios;
    vector<Receiver*> recvs;
	Vector3f sunray_dir;
};

#endif //HELIOSHADOW_SOLARSCENE_H
