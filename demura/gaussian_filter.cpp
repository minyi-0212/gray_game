#include "gaussian_filter.h"
#include <Eigen/Dense>
#include <fstream>

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

const int sample_of_center = 64;
const double center_interval = 1.0 / (sample_of_center-1);
const int sample_of_guass_point = 64;
const double point_interval = 1.0 / sample_of_guass_point;

double gauss_value(const double x, const double y,
	const double cx, const double cy, const double sigmax, const double sigmay)
{
	/*return (1 / (2 * PI*sigma*sigma)) * exp(-((x - cx)*(x - cx)
		+ (y - cy)*(y - cy)) / (2 * sigma*sigma));*/
	return (1 / (2 * PI*sigmax*sigmay)) * exp(-
		((x - cx)*(x - cx) / (2 * sigmax * sigmax) + (y - cy)*(y - cy) / (2 * sigmay * sigmay)));
}

void compute_kernel(VectorXd& k,
	const double cx, const double cy, const double sigmax, const double sigmay)
{
	for (int ki = 0; ki < 9; ki++)
	{
		double v = 0;
		double x = ki % 3, y = ki / 3;
		for (double yi = point_interval / 2; yi < 1; yi += point_interval)
		{
			for (double xi = point_interval / 2; xi < 1; xi += point_interval)
			{
				//cout << x + xi << "," << y + yi << endl;
				v += gauss_value(x + xi, y + yi, cx, cy, sigmax, sigmay);
			}
		}
		k[ki] = v / sample_of_guass_point / sample_of_guass_point;
	}
}

void get_kernels(vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	double cx_begin = 1, cy_begin = 1;
	int index = 0;
	for (double yi = 0; yi <= 1; yi += center_interval)
	{
		for (double xi = 0; xi <= 1; xi += center_interval)
		{
			compute_kernel(kernels[index++],
				cx_begin + xi, cy_begin + yi, 
				sigmax, sigmay);
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
	imwrite("./output_b/kernels_new.png", k_img);
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
	byte t;
	for (int i = 0; i < region.size(); i++)
	{
		for (int j = 0; j < region.size(); j++)
		{
			t = __min(val * kernel[i * 3 + j], 255);
			img.at<Vec3b>(y + region[i], x + region[j]) = Vec3b(t, t, t);
		}
	}
}

void draw_pic_with_scalar(const int y, const int x, const Vec3b val, Mat& img)
{
	for (int i = 0; i < region.size(); i++)
	{
		for (int j = 0; j < region.size(); j++)
		{
			img.at<Vec3b>(y + region[i], x + region[j]) = val;
		}
	}
}

void match(const vector<vector<Point>>& centers_vec,
	const vector<vector<VectorXd>>& data,
	const vector<VectorXd>& kernels)
{
	vector<VectorXd> kernels_normalized;
	kernels_normalized.resize(kernels.size());
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels_normalized[i] = kernels[i].normalized();
	}
	
	Mat img_480 = Mat(Size(img.cols, img.rows), CV_8UC3, Scalar(0, 0, 0)),
		img_ones = Mat(Size(5000, 1500), CV_8UC3, Scalar(0, 0, 0));
	set<int> index_set;
	ofstream out("./output_b/result.csv");
	bool need_indent = true;
	VectorXd ones(9);
	ones << 1, 1, 1,
		1, 1, 1,
		1, 1, 1;
	for (int vj = 0; vj < centers_vec.size(); vj++)
	{
		if (need_indent)
		{
			out << ",";
		}
		need_indent = !need_indent;
		for (int vi = 0; vi < centers_vec[vj].size(); vi++)
		{
			if (data[vj][vi][0] < 0)
			{
				out << -1 << ",,";
				continue;
			}
			int index = find_most_similar(data[vj][vi].normalized(), kernels_normalized);
			index_set.insert(index);
			MatrixXd A(kernels[index]);
			double value = A.colPivHouseholderQr().solve(data[vj][vi])[0];
			draw_pic_with_kernel(centers_vec[vj][vi].y, centers_vec[vj][vi].x, 
				value, kernels[index], img_output);
			draw_pic_with_kernel(centers_vec[vj][vi].y, centers_vec[vj][vi].x, 
				480, kernels[index], img_480);
			img_ones.at<Vec3b>(vj, 2 * vi + vj % 2) = Vec3b(value / 12, 0, 0);
			//out << v.first / IMG_WIDTH << "," << v.first % IMG_WIDTH << "," << value << endl;
			out << value << ",,";
		}
		out << endl;
	}
	out.close();
	cout <<"used kernels count: "<< index_set.size() << endl;
	imwrite("./output_b/result_of_kernels.png", img_output);
	imwrite("./output_b/result_of_480_mul_kernels.png", img_480);
	imwrite("./output_b/B16_result.png", img_ones);
	//cvtColor(img_output, img_output, CV_BGR2GRAY);
	//cvtColor(img_480, img_480, CV_BGR2GRAY);
	//imwrite("./output_b/result_of_kernels.bmp", img_output);
	//imwrite("./output_b/result_of_480.bmp", img_480);
}

double compute_sigma(const vector<vector<VectorXd>>& data,
	vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	Performance p;
	get_kernels(kernels, sigmax, sigmay);
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels[i].normalize();
	}
	double loss = 0, tmp_dot;
	for (int vj = 0; vj < data.size(); vj++)
	{
		for (int vi = 0; vi < data[vj].size(); vi++)
		{
			if (data[vj][vi][0] < 0)
				continue;
			int tmp_k_index = find_most_similar(data[vj][vi].normalized(), kernels);
			VectorXd vk(data[vj][vi].normalized() - kernels[tmp_k_index]);
			loss += vk.dot(vk);
		}
	}
	cout << "compute loss: " << p.end() << "s" << endl;
	cout << "sigma: " << sigmax << "," << sigmay << ", loss: " << loss << endl << endl;
	return loss;
}

//#define SIGMA_COMPUTE
void compute_dumura(vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error)
{
	Performance p;
	img_output = Mat(Size(img.cols, img.rows), CV_8UC3, Scalar(0, 0, 0));
	for (auto ce : centers_error)
	{
		//draw_pic_with_scalar(ce.y, ce.x, Vec3b(0, 0, 255), img_output);
		img_output.at<Vec3b>(ce.y, ce.x) = Vec3b(0, 0, 255);
	}
	vector<VectorXd> kernels(sample_of_center*sample_of_center, VectorXd(9));
#ifdef SIGMA_COMPUTE
	/*ofstream out("./output_b/B16_loss.csv");
	double sigma = 0.5, loss = 100000, loss_old, flag = -1;
	for (; sigma < 2; sigma += 0.01)
	{
		loss_old = loss;
		loss = compute_sigma(data, kernels, sigma);
		if (loss - loss_old > 0 && flag < 0)
		{
			sigma -= 0.01;
			break;
		}
		flag = loss - loss_old;
		out << sigma << "," << loss << endl;
	}
	cout << "sigma:" << sigma << endl;
	out.close();*/
#if 0
	char outfile[MAX_PATH];
	double sigmax = 1.06, sigmay = 0.1, to = 2,
		loss = 1000000, loss_old, flag = -1;
	sprintf_s(outfile, "./output_b/%.2f__%.2f,%.2f_c%d_g%d.csv",
		sigmax, sigmay, to,
		sample_of_center, sample_of_guass_point);
	ofstream out(outfile);
	for (; sigmay < to; sigmay += 0.1)
	{
		loss_old = loss;
		loss = compute_sigma(data, kernels, sigmax, sigmay);
		out << sigmax << "," << sigmay << "," << loss << endl;
	}
	cout << "sigma:" << sigmax <<"," << sigmay << endl;
	out.close();
#else
	char outfile[MAX_PATH];
	double sigmax = 0.1, sigmay = 0.91, to = 2,
		loss = 1000000, loss_old, flag = -1;
	sprintf_s(outfile, "./output_b/%.2f,%.2f__%.2f_c%d_g%d.csv",
		sigmax, to, sigmay,
		sample_of_center, sample_of_guass_point);
	ofstream out(outfile);
	for (; sigmax < to; sigmax += 0.1)
	{
		loss_old = loss;
		loss = compute_sigma(data, kernels, sigmax, sigmay);
		out << sigmax << "," << sigmay << "," << loss << endl;
	}
	cout << "sigma:" << sigmax << "," << sigmay << endl;
	out.close();
#endif

#else
	double sigmax = 1.06, sigmay = 0.92;
	cout << "sigma: " << sigmax << "," << sigmay << endl;
	get_kernels(kernels, sigmax, sigmay);
	combination_kernels(kernels);
	match(centers_vec, data, kernels);
#endif
	p.endAndPrint();
}