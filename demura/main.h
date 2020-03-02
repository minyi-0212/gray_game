#pragma once
#include "Performance.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/core>

struct info
{
	Eigen::VectorXd p;
	Eigen::VectorXd distribution;
};