#pragma once
#include "main.h"

void compute_dumura(const std::unordered_map<int, Eigen::VectorXd>& centers,
	std::vector<cv::Point>& centers_error);
