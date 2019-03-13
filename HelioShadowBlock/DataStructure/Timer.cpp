#include "Timer.h"


void Timer::resetStart()
{
	start = chrono::high_resolution_clock::now();
}

void Timer::printDuration(const string s)
{
	std::cout << s << ": " << getDuration() << "s." << endl;
}


chrono::time_point<chrono::steady_clock> Timer::start = chrono::high_resolution_clock::now();

double Timer::getDuration() {
	auto elapsed = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
	auto time = double(elapsed.count())*chrono::microseconds::period::num / chrono::microseconds::period::den;
	return time;
}