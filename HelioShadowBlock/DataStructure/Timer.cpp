#include "Timer.h"


void Timer::resetStart()
{
	start = chrono::high_resolution_clock::now();
}

void Timer::printDuration(const string s)
{
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	std::cout << s << ": " << time << "s." << endl;
}


chrono::time_point<chrono::steady_clock> Timer::start = chrono::high_resolution_clock::now();