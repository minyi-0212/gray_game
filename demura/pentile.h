#pragma once
#include "main.h"
void rgb2pentile(const cv::Mat& rgb, cv::Mat& pentile);
void pentile2rgb(const cv::Mat& pentile, cv::Mat& rgb);
void draw_pattern(); 
void draw_pattern2(const char* prefix, int w, int h);