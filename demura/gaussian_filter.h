#pragma once
#include "main.h"

//void compute_dumura(const std::map<int, Eigen::VectorXd>& centers,
//	std::vector<cv::Point>& centers_error);
//void compute_dumura(std::vector<std::vector<cv::Point>>& centers_vec,
//	std::vector<std::vector<Eigen::VectorXd>>& data,
//	std::vector<cv::Point>& centers_error, int xy, // x:0, y:1, no sigma compute:3
//	double from, double to, double another, double add, const char* prefix, RGB select_rgb,
//	const char* outputfile);
//
//void compute_dumura_by_combile_single_rgb(
//	std::vector<std::vector<std::vector<std::pair<cv::Point, cv::Point>>>>& centers_vec,
//	std::vector<std::vector<std::vector<Eigen::VectorXd>>>& data,
//	std::vector<std::vector<cv::Point>>& centers_error,
//	std::vector<std::vector<std::pair<cv::Point, cv::Point>>>& cross_points, 
//	std::vector<cv::Mat>& pic,
//	const char* output_prefix, int width, int height);

void compute_dumura(const std::vector<std::vector<std::vector<LED_info>>>& centers_vec,
	const std::vector<int>& capture_pentile_g_value,
	const std::vector<int>& expo,
	const std::vector<cv::Mat>& pic, const int primary_pic,
	const char* output_prefix, int width, int height);


void compute_dumura_single_pic(std::vector<std::vector<std::vector<LED_info>>>& relationship,
	const std::vector<int>& capture_pentile_g_value, const cv::Mat& img, 
	const char* input_prefix, const char* output_prefix, int width, int height);

void compute_intensity_to_exr(std::vector<std::vector<std::vector<LED_info>>>& relationship,
	const cv::Mat& img, const char* input_prefix, const char* output_prefix, int width, int height);