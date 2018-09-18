#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <chrono>
#include <omp.h>

using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

#define Epsilon		1e-6
#define VERTEXSCALE 100000


//#define DEBUG
//#define OUTPUTRES
//#define READFILE
#define CLIPPER