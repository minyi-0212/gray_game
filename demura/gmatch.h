#pragma once
#include <vector>
#include <Eigen/Eigen>
class gmatch
{
public:
	gmatch(float sigma, int ksize = 3, int sample_of_center = 65, int sample_of_guass_point = 64);
	float get_intensity(std::vector<float>& region);
	void reset_sigma(float sigma);
	float get_single_loss(std::vector<float>& region);
	float get_total_loss(std::vector< std::vector<float>>& regions);
private:
	int kernel_size, kernel_size_1d;
	float sigmax, sigmay;
	int sample_of_center, sample_of_guass_point, base;
	float center_interval, point_interval;
	std::vector<Eigen::VectorXf> kernels;
	std::vector<Eigen::VectorXf> kernels_normalized;
	Eigen::VectorXf input, input_normalized;
	Eigen::MatrixXf A;

	void get_kernels();
	void compute_kernel(const int id, const float cx, const float cy);
	float gauss_value(const float x, const float y, const float cx, const float cy);
	void output(const Eigen::VectorXf& a);
};