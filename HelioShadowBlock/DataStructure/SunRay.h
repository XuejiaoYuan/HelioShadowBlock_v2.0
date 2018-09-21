//
// Created by Amber on 2018/4/17.
//

#ifndef HELIOSHADOW_SUNRAY_H
#define HELIOSHADOW_SUNRAY_H

// #include "../Common/utils.h"
//#include "../Common/global_function.h"
#include "../Common/CommonFunc.h"
#include "SPA.h"
//#include <string>
//#include <fstream>
//#include <iostream>
//#include <sstream>
//#include <vector>
//using namespace std;

class SunRay {
public:
    SunRay()= default;
    Vector3f calcSunRay(const string& spa_data_file);
	Vector3f changeSunRay(const vector<int>& time_param);
	Vector3f changeSunRay(const float&altitude, const float&azimuth);
	Vector3i getSunSet();
	float current_altitude, current_azimuth;

private:
	spa_data spa;
    Vector3f sunray_dir;			// From sun to ground
};


#endif //HELIOSHADOW_SUNRAY_H
