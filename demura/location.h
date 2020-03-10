#pragma once
#include "main.h"

//void find_OLED_location(std::map<int, Eigen::VectorXd>& centers,
//	std::vector<cv::Point>& centers_error);

void find_OLED_location(std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error);
void find_OLED_location_with_mask(std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error, bool green = false);