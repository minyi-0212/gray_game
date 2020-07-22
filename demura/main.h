#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Performance.h"
#include "exr.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <map>
#include <unordered_set>
#include <Eigen/core>
#include <set>

enum RGB { BLUE, GREEN, RED };
enum LED_state
{
	VALID, INVALID, CROSS, LINE_ONE, FIRST
};

struct LED_info {
	cv::Point pixel, locate;
	LED_state state;
};

#define SINGLE BLUE
//#define ONLY_LOCATION
//#define SIGMA_COMPUTE
#define INPUT_CROSS_COOR
#define SQRT
