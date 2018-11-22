#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <chrono>
#include <omp.h>
#include <cstdarg>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

#define Epsilon		1e-6
#define VERTEXSCALE 1000000
#define PI acos(float(-1))
#define HELIOSTAT_REFLECTIVITY 0.88 
#define RECEIVER_SLICE 0.05		// *ATTENTION*: w与l应被RECEIVER_SLICE整除

//#define DEBUG
//#define OUTPUTRES
//#define READFILE
#define CLIPPER


typedef enum {
	Initial, GroundMode, ReceiverMode, LayoutMode, HeliostatMode
}InputMode;

typedef enum {
	RectLayoutType, CrossRectLayoutType, FermatLayoutType, RadialLayoutType
}LayoutType;

