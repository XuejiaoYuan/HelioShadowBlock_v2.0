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
    Vector3d calcSunRay(const string& spa_data_file);
	Vector3d changeSunRay(const vector<int>& time_param);
	Vector3d changeSunRay(const double&altitude, const double&azimuth);
	Vector3i getSunSet();
	double calcDNI(const vector<int>& time_param);
	double current_altitude, current_azimuth;
	double current_DNI;
private:
	int calcDay(const vector<int>& time_param);
	spa_data spa;
    Vector3d sunray_dir;			// From sun to ground
};


#endif //HELIOSHADOW_SUNRAY_H
