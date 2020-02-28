#include "gaussian_filter.h"

const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
using namespace std;
using namespace cv;
using namespace Eigen;

extern Mat img;
extern int IMG_WIDTH;
Mat img_output;

double gauss_sample(const int size, const double sigma,
	const int x, const int y, const int center)
{
	double v=0, range= 0.5, interval = 0.1;
	for (double yi = y - range; yi < y + range; yi += interval)
	{
		for (double xi = x - range; xi < x + range; xi += interval)
		{
			v += (1 / (2 * PI*sigma*sigma))*exp(-((xi - center)*(xi - center)
				+ (yi - center)*(yi - center)) / (2 * sigma*sigma));
		}
	}
	return v / (interval * 2 / interval);
}

void get_gaussian_kernel(vector<vector<double>>& gaus, const int size, const double sigma)
{
	if (size==0 || gaus.size() != size+2 || gaus[0].size() != size+2)
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
			gaus[i+1][j+1] = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			//gaus[i + 1][j + 1] = gauss_sample(size, sigma, i, j, center);
			sum += gaus[i+1][j+1];
		}
	}
	for (int i = 1; i < size+1; i++)
	{
		for (int j = 1; j < size+1; j++)
		{
			gaus[i][j] /= sum;
		}
	}
	return;
}

void get_region_kernels(vector<VectorXd>& kernels, const int size, const double sigma)
{
	vector<vector<double>> gaus(size+2, vector<double>(size+2));
	get_gaussian_kernel(gaus, size, sigma);
	/*for (auto g : gaus)
	{
		for (auto gg : g)
		{
			cout << gg << " ";
		}
		cout << endl;
	}*/
	kernels.resize(size*size, VectorXd(9));
	for (int i = 1; i < size + 1; i++)
	{
		for (int j = 1; j < size + 1; j++)
		{
			kernels[(i-1)*size + j-1] <<
				gaus[i - 1][j - 1], gaus[i - 1][j], gaus[i - 1][j + 1],
				gaus[i][j - 1], gaus[i][j], gaus[i][j + 1],
				gaus[i + 1][j - 1], gaus[i + 1][j], gaus[i + 1][j + 1];
			kernels[(i - 1)*size + j - 1].normalize();
		}
	}
	//cout << kernels[0];
}

void combination_kernels(const vector<VectorXd>& kernels)
{
	/*cout << omp_get_num_procs() << endl;
	omp_set_num_threads(12);*/
	int size = sqrt(kernels.size()), kernel_size = 3;
	cout << "kernels count: " << size << endl;
	Mat k(Size(size * kernel_size, size * kernel_size), CV_8UC1);
//#pragma parallel omp for
	for (int y = 0; y < size; y++)
	{
		for (int x = 0; x < size; x++)
		{
			for (int yy = 0; yy < kernel_size; yy++)
			{
				for (int xx = 0; xx < kernel_size; xx++)
				{
					/*cout << x << " " << " " << yy << " " << x * kernel_size + xx << ","
						<< xx << " " << yy << " " << x << " " << yy * kernel_size + xx << endl;*/
					k.at<byte>(y * kernel_size + yy, x * kernel_size + xx)
						= (byte)(kernels[y*size + x][yy*kernel_size + xx] *255);
				}
			}
		}
	}
	imwrite("./output/kernels.png", k);
}

double compute_sigma(const unordered_map<int, VectorXd>& centers,const int size, const double sigma)
{
	//int size = 16, sigma = 10;
	vector<VectorXd> kernels;
	Performance p;
	get_region_kernels(kernels, size, sigma);
	//combination_kernels(kernels);
	cout << "get_kernels: " << p.end() << "s" << endl;
	double loss = 0, tmp_loss;
	p.start();
#pragma parallel omp for
	for (auto v: centers)
	{
		tmp_loss = 100;
		for (int k = 0; k < kernels.size(); k++)
		{
			VectorXd vk(v.second.normalized() - kernels[k]);
			/*cout << "[";
			for (int i = 0; i < 9; i++)
				cout << vk[i] << ", ";
			cout << "]" << endl;
			cout << vk.dot(vk) << endl;*/
			//cout << v.second.normalized().norm() << endl;
			if (vk.dot(vk) < tmp_loss)
				tmp_loss = vk.dot(vk);
		}
		/*if(tmp_loss<0.1)
			cout << tmp_loss << endl;*/
		loss += tmp_loss;
	}
	cout << "compute loss: " << p.end() <<"s"<< endl;
	cout << "sigma: "<< sigma << ", loss: " << loss << endl << endl;
	return loss;
}

vector<int> region({-1,0,1});
void draw(const int y, const int x, const byte val)
{
	for (int i = 0; i < region.size(); i++)
	{
		for (int j = 0; j < region.size(); j++)
		{
			img_output.at<byte>(y+region[i], x+region[j]) = val;
		}
	}
}

void get_result(const unordered_map<int, VectorXd>& centers, const int size, const double sigma)
{
	vector<VectorXd> kernels;
	get_region_kernels(kernels, size, sigma);
	double tmp_loss;
	int index;
	VectorXd s1, s2;
	for (auto v : centers)
	{
		tmp_loss = 100;
		index = v.first;
		for (int k = 0; k < kernels.size(); k++)
		{
			VectorXd vk(v.second.normalized() - kernels[k]);
			//cout << vk.dot(vk) << endl;
			if (vk.dot(vk) < tmp_loss)
			{
				tmp_loss = vk.dot(vk);
				s1 = v.second;
				s2 = vk;
			}
		}
		int result = s1.dot(s2) / s1.normalized().dot(s2);
		cout << result << "," << s1.norm() << endl;
		draw(index / IMG_WIDTH, index % IMG_WIDTH, result);
	}
	imwrite("./output/result.png", img_output);
}

void gaussian(const unordered_map<int, VectorXd>& centers)
{
	img_output= Mat(Size(img.cols, img.rows), CV_8UC1);
	/*vector<VectorXd> centers_vec(centers.size(), VectorXd(9));  //to do->hash°¢unique
	
	for (int i=0;i<centers.size();i++)
	{
		centers_vec[i] << 
			img.at<byte>(centers[i].y-1,centers[i].x-1),
			img.at<byte>(centers[i].y-1,centers[i].x),
			img.at<byte>(centers[i].y-1,centers[i].x+1),
					   
			img.at<byte>(centers[i].y,centers[i].x-1),
			img.at<byte>(centers[i].y,centers[i].x),
			img.at<byte>(centers[i].y,centers[i].x+1),
					   
			img.at<byte>(centers[i].y+1,centers[i].x-1),
			img.at<byte>(centers[i].y+1,centers[i].x),
			img.at<byte>(centers[i].y+1, centers[i].x+1);
		centers_vec[i].normalize();
	}*/
	Performance p;
	// 1.8 0.7
	/*for (double sigma = 1.6; sigma < 2.5; sigma+=0.1) 
	{
		compute_sigma(centers, 16, sigma);
	}*/
	//compute_sigma(centers_vec, 16, 1.8);
	get_result(centers, 16, 0.7);
	p.endAndPrint();
}