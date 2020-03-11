#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Performance.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <map>
#include <unordered_set>
#include <Eigen/core>
#include <set>

struct info
{
	Eigen::VectorXd p;
	Eigen::VectorXd distribution;
};