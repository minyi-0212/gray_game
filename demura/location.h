#pragma once
#include "main.h"

void find_OLED_location(std::unordered_map<int, Eigen::VectorXd>& centers,
	std::vector<cv::Point>& centers_error);
