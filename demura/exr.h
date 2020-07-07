#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "ImfRgbaFile.h"
#include "ImfTiledRgbaFile.h"
#include "ImfArray.h"
#include <vector>

bool save_exr_with_rgb(const char fileName[], int width, int height, 
	const Imf::Array<Imf::Rgba> &pixels, const int compression = 3);
bool save_exr_with_float(const char fileName[], int width, int height, 
	const std::vector<float> &pixel, const int compression = 3);
bool save_exr_with_float_gray(const char fileName[], int width, int height, const std::vector<float>& pixels);
void write_tiled_exr(const char fileName[], int width, int height, const std::vector<float>& pixels,
	int tileWidth = 32, int tileHeight = 32);

void readPixel(const char *fileName, Imf::Rgba& data, int x, int y);
bool read_exr(const char fileName[], std::vector<float> &pixels);
bool read_exr(const char fileName[], int& width, int& height, std::vector<float> &pixels, bool bgr = false);
bool read_exr_in_row(const char fileName[], const int& start_row, const int& read_number_of_lines, int& width, int& height, std::vector<float> &pixels);
bool read_exr_in_row(const char fileName[], const int& start_row, const int& read_number_of_lines, std::vector<float> &pixels);
bool read_exrs_in_row(const std::vector<std::string>& files_name,
	const int& start_row, const int& read_number_of_lines, std::vector<std::vector<float>> &pixels);
bool read_exr_in_block(const char fileName[],
	const int& x_1, const int& x_2,
	const int& y_1, const int& y_2, std::vector<float> &pixels);
void read_tiled_exr(const char fileName[], int &width, int &height, std::vector<float> &pixels);
void read_tiled_exr(const char fileName[], int &width, int &height, std::vector<float> &pixels, const int x0, const int x1, const int y0, const int y1);
void read_tiled_exr(const char fileName[], const int tile_x0, const int tile_x1, const int tile_y0, const int tile_y1, std::vector<float> &pixels);


bool is_file_exist(const char fileName[]);

// multi_channel
bool write_multi_channel_exr(const char fileName[], const int channel_size, const float* gPixels,
	const int width, const int height);

// test deep 
void writeDeepScanlineFile();
