#include "gaussian_filter.h"
#include <Eigen/Dense>
#include <set>

const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
using namespace std;
using namespace cv;
using namespace Eigen;

extern Mat img;
extern int IMG_WIDTH;
Mat img_output;

double gauss_sample(const int size, const double sigma,
	const int x, const int y, const double center)
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
			//gaus[i+1][j+1] = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			gaus[i + 1][j + 1] = gauss_sample(size, sigma, i, j, center);
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
	vector<vector<double>> gaus(size+2, vector<double>(size+2)),
		gaus_add1(size + 3, vector<double>(size + 3));
	get_gaussian_kernel(gaus, size, sigma);
	get_gaussian_kernel(gaus_add1, size+1, sigma);
	/*for (auto g : gaus)
	{
		for (auto gg : g)
		{
			cout << gg << " ";
		}
		cout << endl;
	}*/
	//kernels.resize(size*size, VectorXd(9));
	kernels.clear();
	VectorXd tmp(9);
	//for (int i = 1; i < size + 1; i++)
	for (int i = 1; i < size + 1; i++)
	{
		for (int j = 1; j < size + 1; j++)
		{
			tmp /*kernels[(i-1)*size + j-1]*/ <<
				gaus[i - 1][j - 1], gaus[i - 1][j], gaus[i - 1][j + 1],
				gaus[i][j - 1], gaus[i][j], gaus[i][j + 1],
				gaus[i + 1][j - 1], gaus[i + 1][j], gaus[i + 1][j + 1];
			kernels.push_back(tmp);
			tmp /*kernels[(i-1)*size + j-1]*/ <<
				gaus_add1[i - 1][j - 1], gaus_add1[i - 1][j], gaus_add1[i - 1][j + 1],
				gaus_add1[i][j - 1],	 gaus_add1[i][j],	  gaus_add1[i][j + 1],
				gaus_add1[i + 1][j - 1], gaus_add1[i + 1][j], gaus_add1[i + 1][j + 1];
			kernels.push_back(tmp);
		}
	}
	//cout << kernels[0];
}

void combination_kernels(const vector<VectorXd>& kernels)
{
	int size = sqrt(kernels.size()/2), kernel_size = 3;
	cout << "kernels count: " << size << endl;
	Mat k_img(Size(size * kernel_size, size * kernel_size), CV_8UC3);
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
					int t = (int)(kernels[(y*size + x)*2][yy*kernel_size + xx]  * 255 );
					k_img.at<Vec3b>(y * kernel_size + yy, x * kernel_size + xx)
						= Vec3b(t, t, t);
				}
			}
		}
	}
	imwrite("./output/kernels.png", k_img);
}

void output(const VectorXd& a)
{
	cout << "[";
	for (int i = 0; i < 9; i++)
		cout << a[i] << ",";
	cout << "]" << endl;
}

int find_most_similar(const VectorXd& point, const vector<VectorXd>& kernels)
{
	double tmp_dot = 0;
	int index = 0;
	for (int k = 0; k < kernels.size(); k++)
	{
		/*VectorXd vk(v.second.normalized() - kernels[k].normalized());
		if (vk.dot(vk) < tmp_loss)
			tmp_loss = vk.dot(vk);
		output(v.second.normalized());
		output(kernels[k].normalized());
		cout << tmp_loss << " " << v.second.normalized().dot(kernels[k].normalized()) << endl;*/
		double t = point.dot(kernels[k]);
		if (t > tmp_dot)
		{
			tmp_dot = t;
			index = k;
		}
	}
	return index;
}

double compute_sigma(const unordered_map<int, VectorXd>& centers,
	const int size, const double sigma)
{
	//int size = 16, sigma = 10;
	vector<VectorXd> kernels;
	Performance p;
	get_region_kernels(kernels, size, sigma);
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels[i].normalize();
	}
	//combination_kernels(kernels);
	cout << "get_kernels: " << p.end() << "s" << endl;
	double loss = 0, tmp_dot;
	p.start();
//#pragma parallel omp for
	for (auto v: centers)
	{
		int tmp_k_index = find_most_similar(v.second.normalized(), kernels);
		VectorXd vk(v.second.normalized() - kernels[tmp_k_index]);
		loss += vk.dot(vk);
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
			img_output.at<byte>(y + region[i], x + region[j]) = val;
		}
	}
}

void draw_pic_with_kernel(const int y, const int x, 
	const double val, const VectorXd& kernel)
{
	for (int i = 0; i < region.size(); i++)
	{
		for (int j = 0; j < region.size(); j++)
		{
			img_output.at<byte>(y + region[i], x + region[j]) = val * kernel[i * 3 + j];
		}
	}
}

void get_result(const unordered_map<int, VectorXd>& centers, 
	const int size, const double sigma)
{	
	//omp_set_num_threads(12);
	vector<VectorXd> kernels, kernels_normalized;
	get_region_kernels(kernels, size, sigma);
	kernels_normalized.resize(kernels.size());
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels_normalized[i] = kernels[i].normalized();
	}
	combination_kernels(kernels_normalized);
	//double max_val;
	set<int> index_set;
//#pragma parallel omp for
	for (auto v : centers)
	{
		int index = find_most_similar(v.second.normalized(), kernels_normalized);
		index_set.insert(index);
		MatrixXd A(kernels[index]);
		/*cout << v.first / IMG_WIDTH << " " << v.first % IMG_WIDTH << " "
			<< A.colPivHouseholderQr().solve(v.second)[0] << endl;*/
		/*double val = A.colPivHouseholderQr().solve(v.second)[0];
		draw(v.first / IMG_WIDTH, v.first % IMG_WIDTH, 
			__min(val, 255));
		if (val > max_val)
			max_val = val;*/
		/*draw(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			__min(A.colPivHouseholderQr().solve(v.second)[0], 255));*/

		draw_pic_with_kernel(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			A.colPivHouseholderQr().solve(v.second)[0], kernels[index]);
		/*draw_pic_with_kernel(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			480, kernels[index]);*/

		/*draw(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			s1.dot(s2) / s1.normalized().dot(s2));*/
		/*draw(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			__min(255, s1.maxCoeff() / s2.maxCoeff()));*/

	}
	//cout << "max: " << max_val << endl;
	/*for (auto i : index_set)
		cout << i % 16 << "," << i / 16 << endl;*/
	cout << index_set.size() << endl;
	imwrite("./output/result_of_kernels.png", img_output);
	imwrite("./output/result_of_kernels.bmp", img_output);

	for (auto v : centers)
	{
		int index = find_most_similar(v.second.normalized(), kernels_normalized);
		draw_pic_with_kernel(v.first / IMG_WIDTH, v.first % IMG_WIDTH,
			480, kernels[index]);
	}
	imwrite("./output/result_of_480.png", img_output);
	imwrite("./output/result_of_480.bmp", img_output);
	//imwrite("./output/result.bmp", img_output);
}

void gaussian(const unordered_map<int, VectorXd>& centers)
{
	img_output = Mat(Size(img.cols, img.rows), CV_8UC1);
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
	/*for (double sigma = 0.6; sigma < 0.8; sigma+=0.01)
	{
		compute_sigma(centers, 16, sigma);
	}*/
	//compute_sigma(centers, 16, 0.66);
	get_result(centers, 16, 0.66);
	p.endAndPrint();
}