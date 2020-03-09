#pragma once
#include "main.h"

//void compute_dumura(const std::map<int, Eigen::VectorXd>& centers,
//	std::vector<cv::Point>& centers_error);

void compute_dumura(std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error);
