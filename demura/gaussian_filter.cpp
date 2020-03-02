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

double compute_gauss_value(const double sigma,
	const double x, const double y, const double center)
{
	return (1 / (2 * PI*sigma*sigma)) * exp(-((x - center)*(x - center)
		+ (y - center)*(y - center)) / (2 * sigma*sigma));
}

vector<double> sample({-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5});
double gauss_sample(const int size, const double sigma,
	const int x, const int y, const double center)
{
	double v=0;
	for (int j=0;j< sample.size();j++)
	{
		for (int i = 0; i < sample.size(); i++)
		{
			v += compute_gauss_value(sigma, x + sample[i], y + sample[j], center);
		}
	}
	return v / sample.size() / sample.size();
}

void get_gaussian_kernel(vector<vector<double>>& gaus, const int size, const double sigma)
{
	if (size==0 || gaus.size() != size+2 || gaus[0].size() != size+2)
	{
		cout << "get_gaussian_kernel(): the para not correct!" << endl;
		return;
	}
	double center = (size-1) / 2.0;
	//cout << "center: " << center << endl;
	double sum = 0;
	for (int i = 0; i <= size+1; i++)
	{
		for (int j = 0; j <= size+1; j++)
		{
			//gaus[i][j] = compute_gauss_value(sigma, i - 1, j - 1, center);
			gaus[i][j] = gauss_sample(size, sigma, i-1, j-1, center);
			//sum += gaus[i][j];
		}
	}
	
	Mat tmp= Mat::zeros(Size(size + 2, size + 2), CV_8UC3);
	for (int i = 0; i < size+2; i++)
	{
		for (int j = 0; j < size+2; j++)
		{
			//gaus[i][j] /= sum;
			cout << gaus[i][j] << " ";
			byte t = gaus[i][j] * 1000;
			tmp.at<Vec3b>(j, i) = Vec3b(t,t,t);
		}
		cout << endl;
	}
	cout << endl;
	imwrite("./output/kernels_look.png", tmp);

	return;
}

void get_region_kernels(vector<VectorXd>& kernels, const int size, const double sigma)
{
	vector<vector<double>> gaus(size + 2, vector<double>(size + 2, 0)),
		gaus_add1(size + 3, vector<double>(size + 3, 0));
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
					int t = (int)(kernels[(y*size + x)*2][yy*kernel_size + xx]  * 2550 );
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
			img_output.at<Vec3b>(y + region[i], x + region[j])
				= Vec3b(val, val, val);
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
			img_output.at<Vec3b>(y + region[i], x + region[j]) = 
				Vec3b(val * kernel[i * 3 + j], val * kernel[i * 3 + j], val * kernel[i * 3 + j]);
			//img_output.at<byte>(y + region[i], x + region[j]) = val * kernel[i * 3 + j];
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
	combination_kernels(kernels);
	
	//double max_val;
	set<int> index_set;
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
	img_output = Mat(Size(img.cols, img.rows), CV_8UC3, Scalar(255, 0, 0));
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
	/*for (double sigma = 0.6; sigma < 0.7; sigma+=0.01)
	{
		compute_sigma(centers, 16, sigma);
	}*/
	//compute_sigma(centers, 16, 0.66);
	get_result(centers, 16, 0.67);
	p.endAndPrint();
}