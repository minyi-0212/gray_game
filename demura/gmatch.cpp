#include "gmatch.h"
#include <iostream>

using namespace std;
const float PI = 4.0*atan(1.0);

//#define DEBUG
gmatch::gmatch(float sigma, int ksize, int sample_of_center, int sample_of_guass_point)
	:sigmax(sigma), sigmay(sigma), kernel_size(ksize),
	sample_of_center(sample_of_center), center_interval(1.0 / sample_of_center),
	sample_of_guass_point(sample_of_guass_point), point_interval(1.0 / sample_of_guass_point)
{
	kernel_size_1d = ksize * ksize;
	base = sample_of_center * sample_of_center;
	input = Eigen::VectorXf(kernel_size_1d);
	input_normalized = Eigen::VectorXf(kernel_size_1d);
	kernels.resize(base, Eigen::VectorXf(kernel_size_1d));
	kernels_normalized.resize(base, Eigen::VectorXf(kernel_size_1d));
	A = input;
	get_kernels();
}

float gmatch::get_intensity(vector<float>& region)
{
	try
	{
		if (region.size() != kernel_size_1d)
		{
			//debug info, can delete.
			cout << "[GET_INTENSITY]:input kernel size error: expect " << kernel_size_1d << " not " << region.size() << endl;
			return -1;
		}
		for (int i = 0; i < kernel_size_1d; i++)
		{
			input[i] = region[i];
		}
		input_normalized = input.normalized();

		float tmp_dot = 0, t;
		int index = -1;
		for (int k = 0; k < kernels_normalized.size(); k++)
		{
			t = input_normalized.dot(kernels_normalized[k]);
			if (t > tmp_dot)
			{
				tmp_dot = t;
				index = k;
			}
		}
		if (index == -1)
		{
			//debug info, can delete.
			cout << "[GET_INTENSITY]:can't find result." << endl;
			return -1;
		}
		A = kernels[index];
#ifdef DEBUG
		output(input);
		output(kernels[index]);
#endif
		return A.colPivHouseholderQr().solve(input)[0];
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}
}

void gmatch::reset_sigma(float sigma)
{
	sigmax = sigmay = sigma;
	get_kernels();
}

float gmatch::get_single_loss(vector<float>& region)
{
	try
	{
		if (region.size() != kernel_size_1d)
		{
			//debug info, can delete.
			cout << "[GET_SINGLE_LOSS]:input kernel size error: expect "
				<< kernel_size_1d << " not " << region.size() << endl;
			return -1;
		}
		for (int i = 0; i < kernel_size_1d; i++)
		{
			input[i] = region[i];
		}
		input_normalized = input.normalized();

		float tmp_dot = 0, t;
		int index = -1;
		for (int k = 0; k < kernels_normalized.size(); k++)
		{
			t = input_normalized.dot(kernels_normalized[k]);
			if (t > tmp_dot)
			{
				tmp_dot = t;
				index = k;
			}
		}
		if (index == -1)
		{
			//debug info, can delete.
			cout << "[GET_SINGLE_LOSS]:can't match kernels :" << endl;
			output(input);
			return -1;
		}
		// L2 loss, can be change
		Eigen::VectorXf vk(input_normalized - kernels[index]);
		return vk.dot(vk);
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}
}

float gmatch::get_total_loss(vector<vector<float>>& regions)
{
	try
	{
		if (regions.size() == 0 || regions[0].size() != kernel_size_1d)
		{
			//debug info, can delete.
			cout << "[GET_TOTAL_LOSS]: size error!" << endl;
			return -1;
		}
		float tmp_dot, t, loss = 0;
		int index;
		for (auto region : regions)
		{
			for (int i = 0; i < kernel_size_1d; i++)
			{
				input[i] = region[i];
			}
			input_normalized = input.normalized();
			tmp_dot = 0, index = -1;
			for (int k = 0; k < kernels_normalized.size(); k++)
			{
				t = input_normalized.dot(kernels_normalized[k]);
				if (t > tmp_dot)
				{
					tmp_dot = t;
					index = k;
				}
			}
			if (index == -1)
			{
				//debug info, can delete.
				cout << "[GET_TOTAL_LOSS]:can't match kernels :" << endl;
				output(input);
				//return -1;
			}
			else
			{
				// L2 loss, can be change
				Eigen::VectorXf vk(input_normalized - kernels[index]);
				loss += vk.dot(vk);
			}
		}
		return loss;
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}
}

float gmatch::gauss_value(const float x, const float y, const float cx, const float cy)
{
	return (1 / (2 * PI*sigmax*sigmay)) * exp(-
		((x - cx)*(x - cx) / (2 * sigmax * sigmax) + (y - cy)*(y - cy) / (2 * sigmay * sigmay)));
}

void gmatch::compute_kernel(const int id, const float cx, const float cy)
{
	for (int ki = 0; ki < kernel_size*kernel_size; ki++)
	{
		float v = 0;
		float x = ki % kernel_size, y = ki / kernel_size;
		for (float yi = point_interval / 2; yi < 1; yi += point_interval)
		{
			for (float xi = point_interval / 2; xi < 1; xi += point_interval)
			{
				v += gauss_value(x + xi, y + yi, cx, cy);
			}
		}
		kernels[id][ki] = v / sample_of_guass_point / sample_of_guass_point;
	}
}

void gmatch::get_kernels()
{
	float cx_begin = 1, cy_begin = 1;
	int index = 0;
	for (float yi = center_interval / 2; yi <= 1; yi += center_interval)
	{
		for (float xi = center_interval / 2; xi <= 1; xi += center_interval)
		{
			compute_kernel(index, cx_begin + xi, cy_begin + yi);
			kernels_normalized[index] = kernels[index].normalized();
			index++;
		}
	}
}

void gmatch::output(const Eigen::VectorXf& a)
{
	cout << "[";
	for (int i = 0; i < kernel_size_1d - 1; i++)
	{
		cout << a[i] << ",";
		if ((i + 1) % kernel_size == 0)
			cout << endl;
	}
	cout << a[kernel_size_1d - 1] << "]" << endl;
}