#pragma once
#include "../Common/CommonFunc.h"

class Timer {
public:
	static chrono::time_point<chrono::steady_clock> start;
	static void resetStart();
	static void printDuration(const string s);
	static double getDuration();
};
