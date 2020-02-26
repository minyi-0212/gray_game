#include "gaussian_filter.h"
#include <Eigen/core>

const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
using namespace std;
using namespace cv;
using namespace Eigen;

extern Mat img;

void get_gaussian_kernel(vector<vector<double>>& gaus, const int size, const double sigma)
{
	if (size==0 || gaus.size() != size || gaus[0].size() != size)
	{
		cout << "get_gaussian_kernel(): the para not correct!" << endl;
		return;
	}
	double center = (size-1) / 2.0;
	double sum = 0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			gaus[i][j] = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			sum += gaus[i][j];
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			gaus[i][j] /= sum;
		}
	}
	return;
}

const int half_kernel_size = 1;  // 3*3, half is 1
void get_region_kernels(vector<VectorXd>& kernels, const int size, const int sigma)
{
	vector<vector<double>> gaus(size, vector<double>(size));
	get_gaussian_kernel(gaus, size, sigma);
	/*for (auto g : gaus)
	{
		for (auto gg : g)
		{
			cout << gg << " ";
		}
		cout << endl;
	}*/
	int total_kernels_cnt = size - 2 * half_kernel_size;
	kernels.resize(total_kernels_cnt*total_kernels_cnt, VectorXd(9));
	for (int i = half_kernel_size; i < size - half_kernel_size; i++)
	{
		for (int j = half_kernel_size; j < size - half_kernel_size; j++)
		{
			kernels[(i - half_kernel_size)*total_kernels_cnt + j - half_kernel_size] <<
				gaus[i - 1][j - 1], gaus[i - 1][j], gaus[i - 1][j + 1],
				gaus[i][j - 1], gaus[i][j], gaus[i][j + 1],
				gaus[i + 1][j - 1], gaus[i + 1][j], gaus[i + 1][j + 1];
		}
	}
	//cout << kernels[0];
}

void compute_sigma(const vector<VectorXd>& centers_vec,const int size, const int sigma)
{
	//int size = 16, sigma = 10;
	vector<VectorXd> kernels;
	Performance p;
	get_region_kernels(kernels, size, sigma);
	cout << "get_kernels: " << p.end() << "s" << endl;
	double loss = 0, tmp_loss;
	cout << centers_vec.size() << endl;
	p.start();
#pragma parallel omp for
	for (int v = 0; v < centers_vec.size(); v++)
	{
		tmp_loss = 100;
		for (int k = 0; k < kernels.size(); k++)
		{
			VectorXd vk(centers_vec[v] - kernels[k]);
			cout << vk.dot(vk) << endl;
			if (vk.dot(vk) < tmp_loss)
				tmp_loss = vk.dot(vk);
		}
		if(tmp_loss<0.1)
			cout << tmp_loss << endl;
		loss += tmp_loss;
	}
	cout << "compute loss: " << p.end() <<"s"<< endl;
	cout <<"sigma: "<< sigma << ", loss: " << loss << endl << endl;
}

void gaussian(const vector<Point>& centers)
{
	vector<VectorXd> centers_vec(centers.size(), VectorXd(9));  //to do->hash°¢unique
	
	for (int i=0;i<centers.size();i++)
	{
		centers_vec[i] << 
			img.at<Vec3b>(centers[i].y-1,centers[i].x-1)[0],
			img.at<Vec3b>(centers[i].y-1,centers[i].x)[0],
			img.at<Vec3b>(centers[i].y-1,centers[i].x+1)[0],
					   
			img.at<Vec3b>(centers[i].y,centers[i].x-1)[0],
			img.at<Vec3b>(centers[i].y,centers[i].x)[0],
			img.at<Vec3b>(centers[i].y,centers[i].x+1)[0],
					   
			img.at<Vec3b>(centers[i].y+1,centers[i].x-1)[0],
			img.at<Vec3b>(centers[i].y+1,centers[i].x)[0],
			img.at<Vec3b>(centers[i].y+1, centers[i].x+1)[0];
		centers_vec[i].normalize();
	}
	Performance p;
	compute_sigma(centers_vec, 16, 1);
	compute_sigma(centers_vec, 16, 2);
	compute_sigma(centers_vec, 16, 3);
	compute_sigma(centers_vec, 16, 4);
	compute_sigma(centers_vec, 16, 5);
	p.endAndPrint();
}