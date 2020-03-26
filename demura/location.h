#pragma once
#include "main.h"

//void find_OLED_location(std::map<int, Eigen::VectorXd>& centers,
//	std::vector<cv::Point>& centers_error);

void find_OLED_location(std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error);
void find_OLED_location_with_mask(
	const char *infile_name, const char *mask_file_name,
	const char *outfile_selected_point_name, const char *outfile_set_3x3_region_name,
	std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error, bool green = false);
void find_OLED_cross(const char *img_file, const char *cross_file,
	const char *output_select_point_file, 
	std::vector<std::vector<cv::Point>>& centers_vec,
	std::vector<std::vector<Eigen::VectorXd>>& data,
	std::vector<cv::Point>& centers_error);
void find_OLED_location_with_rgb_combination(
	std::vector<cv::Mat>& rgb, cv::Mat& mask,
	const char *cross_file,
	const char *output_selected_points_prefix,
	int pentile_height,
	std::vector<std::vector<std::vector<std::pair<cv::Point, cv::Point>>>>& centers_vec,
	std::vector<std::vector<std::vector<Eigen::VectorXd>>>& data,
	std::vector<std::vector<cv::Point>>& centers_error,
	std::vector<std::vector<std::pair<cv::Point, cv::Point>>> cross_points);