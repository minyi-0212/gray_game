#pragma once
#include "main.h"

//void compute_dumura(const std::map<int, Eigen::VectorXd>& centers,
//	std::vector<cv::Point>& centers_error);

void compute_dumura(std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error, int xy, // x:0, y:1, no sigma compute:3
	double from, double to, double another, double add, const char* prefix, RGB select_rgb,
	const char* outputfile);
