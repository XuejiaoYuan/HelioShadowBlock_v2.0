//
// Created by Amber on 2018/4/17.
//

#include "SunRay.h"
#include "SPA.h"

Vector3f SunRay::calcSunRay(const string &spa_data_file) {
    fstream inFile(spa_data_file);
    if(inFile.fail()){
        cout << "Can't open the file!"<< endl;
        return Vector3f(0, 0, 0);
    }

    string str, word;
    stringstream ss;
    while(getline(inFile, str)){
        ss.clear();
        ss.str(str);
        ss >> word;
        if(word == "year")
            ss >> spa.year;
        else if(word == "month")
            ss >> spa.month;
        else if(word == "day")
            ss >> spa.day;
        else if(word == "hour")
            ss >> spa.hour;
        else if(word == "minute")
            ss >> spa.minute;
        else if(word == "second")
            ss >> spa.second;
        else if(word == "timezone")
            ss >> spa.timezone;
        else if(word == "delta_ut1")
            ss >> spa.delta_ut1;
        else if(word == "delta_t")
            ss >> spa.delta_t;
        else if(word == "longitude")
            ss >> spa.longitude;
        else if(word == "latitude")
            ss >> spa.latitude;
        else if(word == "elevation")
            ss >> spa.elevation;
        else if(word == "pressure")
            ss >> spa.pressure;
        else if(word == "temperature")
            ss >> spa.temperature;
        else if(word == "slope")
            ss >> spa.slope;
        else if(word == "azm_rotation")
            ss >> spa.azm_rotation;
        else if(word == "atmos_refract")
            ss >> spa.atmos_refract;
        else if(word == "function")
            ss >> spa.function;
    }

	spa.function = SPA_ALL;
    int res = spa_calculate(&spa);
    if(res == 0){
		double altitude = deg2rad(90 - spa.zenith);
		double azimuth = deg2rad(spa.azimuth-90);
		sunray_dir.x() = cos(altitude) * sin(azimuth);
		sunray_dir.y() = -sin(altitude);
		sunray_dir.z() =  -cos(altitude) * cos(azimuth);
    }
	sunray_dir = sunray_dir.normalized();
    inFile.close();

	float min = 60.0*(spa.sunset - (int)(spa.sunset));
	float sec = 60.0*(min - (int)min);
	printf("Sunset:        %02d:%02d:%02d Local Time\n", (int)(spa.sunset), (int)min, (int)sec);

    return sunray_dir;
}

Vector3f SunRay::changeSunRay(const vector<int>& time_param)
{
	spa.month = time_param[0];
	spa.day = time_param[1];
	spa.hour = time_param[2];
	spa.minute = time_param[3];
	int res = spa_calculate(&spa);
	current_altitude = 90 - spa.zenith;
	current_azimuth = spa.azimuth;
	//if (spa.hour > (int)spa.sunset ||
	//	(spa.hour==(int)spa.sunset && 
	//		spa.minute>=60.0*(spa.sunset - (int)(spa.sunset))))
	//	return make_float3(-1, -1, -1);
	if (res == 0) {
		double altitude = deg2rad(90 - spa.zenith);
		double azimuth = deg2rad(spa.azimuth - 90);
		sunray_dir.x() = cos(altitude) * sin(azimuth);
		sunray_dir.y() = -sin(altitude);
		sunray_dir.z() = -cos(altitude) * cos(azimuth);
	}
	sunray_dir = sunray_dir.normalized();
	//cout << "Month: " << spa.month << " Day: " << spa.day << " Hour: " << spa.hour << endl;
	//cout << "altitude: " << 90 - spa.zenith << " azimuth: " << spa.azimuth << endl;

	return sunray_dir;
}

Vector3f SunRay::changeSunRay(const float & altitude, const float & azimuth)
{

	float local_altitude = deg2rad(altitude);
	float local_azimuth = deg2rad(azimuth - 90);
	sunray_dir.x() = cos(local_altitude) * sin(local_azimuth);
	sunray_dir.y() = -sin(local_altitude);
	sunray_dir.z() =  -cos(local_altitude) * cos(local_azimuth);
	sunray_dir = sunray_dir.normalized();
	return sunray_dir;
}
