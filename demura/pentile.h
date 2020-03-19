#pragma once
#include "main.h"
void rgb2pentile(const cv::Mat& rgb, cv::Mat& pentile);
void pentile2rgb(const cv::Mat& pentile, cv::Mat& rgb);
void find_pentile_rgb_relationship();