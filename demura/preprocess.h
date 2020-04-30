#pragma once
#include "main.h"
void preprocess(std::vector<cv::Mat>& rgb, int low_limit, const char *outfile, cv::Mat& mask);
void preprocess(cv::Mat& img, int low_limit, const char *outfile, cv::Mat& mask);