#pragma once
#include <windows.h>
class Performance
{
private:
	LARGE_INTEGER start_time;
	LARGE_INTEGER end_time;
	LARGE_INTEGER frequency;
	double time;
public:
	Performance();
	~Performance();
	void start();
	void restart();
	double end();
	double endAndPrint();
};

