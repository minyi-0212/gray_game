#include "gaussian_filter.h"
#include <Eigen/Dense>
#include <set>
#include <direct.h>

const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
using namespace std;
using namespace cv;
using namespace Eigen;

extern Mat img;
extern int IMG_WIDTH;
Mat img_output;

void output(const VectorXd& a)
{
	cout << "[";
	for (int i = 0; i < 8; i++)
		cout << a[i] << ",";
	cout << a[8] << "]" << endl;
}

const int sample_of_center = 65;
const int sample_of_guass_point = 65;
const double center_interval = 1.0 / (sample_of_center - 1);
const double point_interval = 1.0 / (sample_of_guass_point - 1);

double gauss_value(const double x, const double y,
	const double cx, const double cy, const double sigma)
{
	return (1 / (2 * PI*sigma*sigma)) * exp(-((x - cx)*(x - cx)
		+ (y - cy)*(y - cy)) / (2 * sigma*sigma));
}

void compute_kernel(VectorXd& k,
	const double cx, const double cy, const double sigma)
{
	for (int ki = 0; ki < 9; ki++)
	{
		double v = 0;
		double x = ki % 3, y = ki / 3;
		for (double yi = 0; yi <= 1; yi += point_interval)
		{
			for (double xi = 0; xi <= 1; xi += point_interval)
			{
				v += gauss_value(x + xi, y + yi, cx, cy, sigma);
			}
		}
		k[ki] = v / sample_of_guass_point / sample_of_guass_point;
	}
}

void get_kernels(vector<VectorXd>& kernels, const double sigma)
{
	double cx_begin = 1, cy_begin = 1;
	int index = 0;
	for (double yi = 0; yi <= 1; yi += center_interval)
	{
		for (double xi = 0; xi <= 1; xi += center_interval)
		{
			compute_kernel(kernels[index++], cx_begin + xi, cy_begin + yi, sigma);
		}
	}
}


void combination_kernels(const vector<VectorXd>& kernels)
{
	int size = sample_of_center, kernel_size = 3, scale = 1000;
	cout << "kernels count: " << size << "*" << size << endl;
	Mat k_img(Size(size * kernel_size, size * kernel_size), CV_8UC3);
	for (int y = 0; y < size; y++)
	{
		for (int x = 0; x < size; x++)
		{
			for (int yy = 0; yy < kernel_size; yy++)
			{
				for (int xx = 0; xx < kernel_size; xx++)
				{
					int t = __min((int)(kernels[y*size + x][yy*kernel_size + xx] * scale),255);
					k_img.at<Vec3b>(y * kernel_size + yy, x * kernel_size + xx) = Vec3b(t, t, t);
				}
			}
		}
	}
	imwrite("./output/kernels_new.png", k_img);
}

int find_most_similar(const VectorXd& point, const vector<VectorXd>& kernels)
{
	double tmp_dot = 0, t;
	int index = 0;
	for (int k = 0; k < kernels.size(); k++)
	{
		t = point.dot(kernels[k]);
		if (t > tmp_dot)
		{
			tmp_dot = t;
			index = k;
		}
	}
	return index;
}

const vector<int> region({ -1,0,1 });
void draw_pic_with_kernel(const int y, const int x,
	const double val, const VectorXd& kernel, Mat& img)
{
	for (int i = 0; i < region.size(); i++)
	{
		for (int j = 0; j < region.size(); j++)
		{
			img.at<Vec3b>(y + region[i], x + region[j]) =
				Vec3b(val * kernel[i * 3 + j], 
					val * kernel[i * 3 + j], 
					val * kernel[i * 3 + j]);
		}
	}
}

void match(const unordered_map<int, VectorXd>& centers,
	const vector<VectorXd>& kernels)
{
	vector<VectorXd> kernels_normalized;
	kernels_normalized.resize(kernels.size());
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels_normalized[i] = kernels[i].normalized();
	}
	
	Mat img_480 = Mat(Size(img.cols, img.rows), CV_8UC3, Scalar(0, 0, 0));
	set<int> index_set;
	for (auto v : centers)
	{
		int index = find_most_similar(v.second.normalized(), kernels_normalized);
		index_set.insert(index);
		MatrixXd A(kernels[index]);
		draw_pic_with_kernel(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			A.colPivHouseholderQr().solve(v.second)[0], kernels[index], img_output);
		draw_pic_with_kernel(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			480, kernels[index], img_480);
	}
	cout <<"used kernels count: "<< index_set.size() << endl;
	imwrite("./output/result_of_kernels.png", img_output);
	imwrite("./output/result_of_480.png", img_480);
	cvtColor(img_output, img_output, CV_BGR2GRAY);
	cvtColor(img_480, img_480, CV_BGR2GRAY);
	imwrite("./output/result_of_kernels.bmp", img_output);
	imwrite("./output/result_of_480.bmp", img_480);
}

double compute_sigma(const unordered_map<int, VectorXd>& centers,
	vector<VectorXd>& kernels, const double sigma)
{
	Performance p;
	get_kernels(kernels, sigma);
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels[i].normalize();
	}
	double loss = 0, tmp_dot;
	for (auto v : centers)
	{
		int tmp_k_index = find_most_similar(v.second.normalized(), kernels);
		VectorXd vk(v.second.normalized() - kernels[tmp_k_index]);
		loss += vk.dot(vk);
	}
	cout << "compute loss: " << p.end() << "s" << endl;
	cout << "sigma: " << sigma << ", loss: " << loss << endl << endl;
	return loss;
}

//#define SIGMA_COMPUTE
void compute_dumura(const unordered_map<int, VectorXd>& centers)
{
	Performance p;
	_mkdir("./output");
	img_output = Mat(Size(img.cols, img.rows), CV_8UC3, Scalar(0, 0, 0));
	vector<VectorXd> kernels(sample_of_center*sample_of_center, VectorXd(9));
#ifdef SIGMA_COMPUTE
	double sigma = 0.5, loss = 100000, loss_old, flag = -1;
	for (; sigma < 2; sigma += 0.01)
	{
		loss_old = loss;
		loss = compute_sigma(centers, kernels, sigma);
		if (loss - loss_old > 0 && flag < 0)
		{
			sigma -= 0.01;
			break;
		}
		flag = loss - loss_old;
	}
	cout << "sigma:" << sigma << endl;
#else
	double sigma = 0.6;
#endif
	get_kernels(kernels, sigma);
	combination_kernels(kernels);
	match(centers, kernels);
	p.endAndPrint();
}