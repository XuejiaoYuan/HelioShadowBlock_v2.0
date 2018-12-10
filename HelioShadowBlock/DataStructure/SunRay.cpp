//
// Created by Amber on 2018/4/17.
//

#include "SunRay.h"
#include "SPA.h"

Vector3d SunRay::calcSunRay(const string &spa_data_file) {
    fstream inFile(spa_data_file);
    if(inFile.fail()){
        cout << "Can't open the file!"<< endl;
        return Vector3d(0, 0, 0);
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
		sunray_dir.x() = cos(altitude) * cos(azimuth);
		sunray_dir.y() = -sin(altitude);
		sunray_dir.z() = -cos(altitude) * sin(azimuth);
    }
	sunray_dir = sunray_dir.normalized();
    inFile.close();

	double min = 60.0*(spa.sunset - (int)(spa.sunset));
	double sec = 60.0*(min - (int)min);
	printf("Sunset:        %02d:%02d:%02d Local Time\n", (int)(spa.sunset), (int)min, (int)sec);

    return sunray_dir;
}

Vector3d SunRay::changeSunRay(const vector<int>& time_param)
{
	spa.month = time_param[0];
	spa.day = time_param[1];
	spa.hour = time_param[2];
	spa.minute = time_param[3];
	int res = spa_calculate(&spa);
	current_altitude = 90 - spa.zenith;
	current_azimuth = spa.azimuth;
	calcDNI(time_param);

	if (res == 0) {
		double altitude = deg2rad(90 - spa.zenith);
		double azimuth = deg2rad(spa.azimuth - 90);
		sunray_dir.x() = cos(altitude) * cos(azimuth);
		sunray_dir.y() = -sin(altitude);
		sunray_dir.z() = -cos(altitude) * sin(azimuth);
	}
	sunray_dir = sunray_dir.normalized();

	return sunray_dir;
}

Vector3d SunRay::changeSunRay(const double & altitude, const double & azimuth)
{
	double local_altitude = deg2rad(altitude);
	double local_azimuth = deg2rad(azimuth - 90);
	sunray_dir.x() = cos(local_altitude) * sin(local_azimuth);
	sunray_dir.y() = -sin(local_altitude);
	sunray_dir.z() =  -cos(local_altitude) * cos(local_azimuth);
	sunray_dir = sunray_dir.normalized();
	return sunray_dir;
}

Vector3i SunRay::getSunSet()
{
	double min = 60.0*(spa.sunset - (int)(spa.sunset));
	double sec = 60.0*(min - (int)min);
	printf("Sunset:        %02d:%02d:%02d Local Time\n", (int)(spa.sunset), (int)min, (int)sec);
	Vector3i sunset(int(spa.sunset), int(min), int(sec));
	return sunset;
}

double SunRay::calcDNI(const vector<int>& time_param)
{
	if (this->current_altitude < 0)
		return 0;

	double ALT = 0.0;  // 当地海拔, 默认为0, 单位km
	int day_count = calcDay(time_param);
	double E0 = 1.353 + 0.045 * cos(2 * PI* (day_count + 10) / 365); // 太阳辐射进入地球大气之后单位面积的能量，单位kW / m2

	double a = 0.4237 - 0.00821 * (6 - ALT) * (6 - ALT);
	double b = 0.5055 + 0.00595 * (6.5 - ALT) * (6.5 - ALT);
	double c = 0.2711 + 0.01858 * (2.5 - ALT) * (2.5 - ALT);

	double fair = this->current_altitude * PI / 180.0;
	double DNI = 1000 * E0 * (a + b * exp(-c / sin(fair)));
	// double DNI = E0 * (a + b * exp(-c / sin(fair)));  // 为减小数量级修改了计算式
	this->current_DNI = DNI;
	return DNI;
}

int SunRay::calcDay(const vector<int>& time_param)
{
	int year = 2018;
	int month = time_param[0];
	int day = time_param[1];
	vector<int> Month = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	bool flag = false;
	if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))
		flag = true;

	int cnt_day = 0;
	for (int i = 0; i < month - 1; i++)
		cnt_day += Month[i];
	if (flag && month > 2)
		cnt_day += 1;
	cnt_day += day;
	return cnt_day;
}
