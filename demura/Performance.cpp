#include "Performance.h"

#include <iostream>
using std::cout;
using std::endl;

Performance::Performance():time(0)
{
	QueryPerformanceCounter(&start_time);
}


Performance::~Performance()
{
}


void Performance::start() {
	QueryPerformanceCounter(&start_time);
}

void Performance::restart() {
	QueryPerformanceCounter(&start_time);
}

double Performance::end() {
	QueryPerformanceCounter(&end_time);
	if (!QueryPerformanceFrequency(&frequency))
		return -1;
	return (double)(end_time.QuadPart - start_time.QuadPart) / (double)frequency.QuadPart;
}


double Performance::endAndPrint() {
	QueryPerformanceCounter(&end_time);
	if (!QueryPerformanceFrequency(&frequency))
		return -1;
	time = (double)(end_time.QuadPart - start_time.QuadPart) / (double)frequency.QuadPart;
	cout << "Running time: " << time * 1000 << "ms" << endl;
	return time;
}