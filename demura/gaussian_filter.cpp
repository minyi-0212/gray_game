#include "gaussian_filter.h"
#include <Eigen/Dense>
#include <fstream>
#include <flann/flann.hpp>
#include <flann/util/matrix.h>
#include <unordered_map>

const double PI = 4.0*atan(1.0); //‘≤÷‹¬ ¶–∏≥÷µ
using namespace std;
using cv::Mat;
using cv::Vec3b;
using cv::Point;
using cv::Size;
using cv::Scalar;
using namespace Eigen;

//extern Mat img;
extern int IMG_WIDTH, IMG_HEIGHT;
Mat img_output;

void output(const VectorXd& a)
{
	cout << "[";
	for (int i = 0; i < 8; i++)
		cout << a[i] << ",";
	cout << a[8] << "]" << endl;
}

const int sample_of_center = 64, base = sample_of_center* sample_of_center;
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

void rotate_xy(double& x, double& y, const double rx, const double ry, const double theta)
{
	double c = cos(theta*PI / 180), s = sin(theta*PI / 180),
		xx = x - rx, yy = y - ry;
	x = c * xx - s * yy + rx;
	y = s * xx + c * yy + ry;
}

void compute_kernel(VectorXd& k,
	const double cx, const double cy,
	const double sigmax, const double sigmay)
{
	for (int ki = 0; ki < 9; ki++)
	{
		double v = 0;
		double x = ki % 3, y = ki / 3;
		for (double yi = point_interval / 2; yi < 1; yi += point_interval)
		{
			for (double xi = point_interval / 2; xi < 1; xi += point_interval)
			{
				v += gauss_value(x + xi, y + yi, cx, cy, sigmax, sigmay);
			}
		}
		k[ki] = v / sample_of_guass_point / sample_of_guass_point;
	}
}

void compute_kernel(VectorXd& k,
	const double cx, const double cy, 
	const double sigmax, const double sigmay,
	const double theta)
{
	for (int ki = 0; ki < 9; ki++)
	{
		double v = 0;
		double x = ki % 3, y = ki / 3, rotate_x, rotate_y;
		for (double yi = point_interval / 2; yi < 1; yi += point_interval)
		{
			for (double xi = point_interval / 2; xi < 1; xi += point_interval)
			{
				rotate_x = x + xi;
				rotate_y = y + yi;
				/*Vector2d tmp1, tmp2;
				tmp1 << rotate_x - cx, rotate_y - cy;
				cout << rotate_x << "," << rotate_y << "  "
					<< (rotate_x - cx)*(rotate_x - cx) + (rotate_y - cy)* (rotate_y - cy) << " -> ";*/
				rotate_xy(rotate_x, rotate_y, cx, cy, theta);
				/*tmp2 << rotate_x - cx, rotate_y - cy;
				cout << rotate_x << "," << rotate_y << "  "
					<< (rotate_x - cx)*(rotate_x - cx) + (rotate_y - cy)* (rotate_y - cy) << "  "
					<< endl;
				tmp1.normalize();
				tmp2.normalize();
				cout << tmp1.dot(tmp2) <<"," << cos(45*PI/180)<< endl;*/
				//v += gauss_value(x + xi, y + yi, cx, cy, sigmax, sigmay);
				v += gauss_value(rotate_x, rotate_y, cx, cy, sigmax, sigmay);
			}
		}
		k[ki] = v / sample_of_guass_point / sample_of_guass_point;
	}
}

void rotate_kernel_degree_90(const VectorXd& k, VectorXd& to)
{
	to << k(6), k(3), k(0), k(7), k(4), k(1), k(8), k(5), k(2);
}

void get_kernels(vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	double cx_begin = 1, cy_begin = 1;
	int index = 0;
	if (sigmax == sigmay)
	{
		for (double yi = 0; yi <= 1; yi += center_interval)
		{
			for (double xi = 0; xi <= 1; xi += center_interval)
			{
				compute_kernel(kernels[index++], cx_begin + xi, cy_begin + yi, sigmax, sigmay);
			}
		}
	}
	else
	{
		for (double yi = 0; yi <= 1; yi += center_interval)
		{
			for (double xi = 0; xi <= 1; xi += center_interval)
			{
				for (int ri = 0; ri < 4; ri++)
				{
					compute_kernel(kernels[index++],
						cx_begin + xi, cy_begin + yi,
						sigmax, sigmay, ri * 45);
				}
			}
		}
	}
	//int base = index;
	/*cout << "need rotate " << kernels.size() / base - 1 << " times." << endl;
	while (index < kernels.size())
	{
		rotate_kernel_degree_90(kernels[index - base], kernels[index++]);
		for (int ii =0;ii< index-1;ii++)
			if (kernels[index - 1] == kernels[ii])
			{
				cout << "equal " << (index-1) / base <<","<< (index-1) % base << endl;
				break;
			}
	}
	cout << index << "," << sample_of_center * sample_of_center << endl;*/
}

void combination_kernels(const vector<VectorXd>& kernels)
{
	int size = sample_of_center, kernel_size = 3, scale = 1000;
	cout << "combination kernels count: " << size << "*" << size << endl;
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

void query_by_kdtree
(const vector<vector<VectorXd>>& data, const vector<VectorXd>& kernels,
	flann::Matrix<int>& indices, flann::Matrix<double>& dists)
{
	vector<double> kd_tree_kernel(kernels.size() * 9, 0),
		kd_tree_data;
	for (int i = 0; i < kernels.size(); i++)
	{
		VectorXd tmp = kernels[i];
		tmp.normalize();
		for (int mi = 0; mi < 9; mi++)
			kd_tree_kernel[i * 9 + mi] = tmp[mi];
	}
	for (auto d : data)
		for (auto dd : d)
		{
			dd.normalize();
			for (int mi = 0; mi < 9; mi++)
				kd_tree_data.push_back(dd[mi]);
		}
	flann::Matrix<double> points_mat = flann::Matrix<double>(&kd_tree_kernel[0], kernels.size(), 9);
	flann::Matrix<double> query = flann::Matrix<double>(&kd_tree_data[0], kd_tree_data.size() / 9, 9);
	int nns_number = 1;
	//flann::Matrix<int> indices(new int[query.rows*nns_number], query.rows, nns_number);
	//flann::Matrix<double> dists(new double[query.rows*nns_number], query.rows, nns_number);
	indices = flann::Matrix<int>(new int[query.rows*nns_number], query.rows, nns_number);
	dists = flann::Matrix<double>(new double[query.rows*nns_number], query.rows, nns_number);
	flann::Index<flann::L2<double>> index(points_mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	index.knnSearch(query, indices, dists, nns_number, flann::SearchParams(128));
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
	
	Mat img_480 = Mat(Size(IMG_WIDTH, IMG_HEIGHT), CV_8UC3, Scalar(0, 0, 0)),
		img_ones = Mat(Size(2500, 1250), CV_8UC3, Scalar(0, 0, 0));
	set<int> index_set;
	ofstream out("./output/result.csv");
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
			/*Vec3b color(0, 0, 0);
			if (index % 4 == 0)
				color = Vec3b(value / 12, 0, 0);
			else if (index % 4 == 1)
				color = Vec3b(0, value / 12, 0);
			else if (index % 4 == 2)
				color = Vec3b(0, 0, value / 12);
			else if (index % 4 == 3)
				color = Vec3b(0, value / 12, value / 12);*/
			Vec3b color(value / 12, value / 12, value / 12);
			img_ones.at<Vec3b>(vj, 2 * vi + vj % 2 /*vi*/) = color;
			//out << v.first / IMG_WIDTH << "," << v.first % IMG_WIDTH << "," << value << endl;
			out << value << ",,";
		}
		out << endl;
	}
	out.close();
	cout <<"used kernels count: "<< index_set.size() << endl;
	imwrite("./output/result_of_kernels.png", img_output);
	imwrite("./output/result_of_480_mul_kernels.png", img_480);
	imwrite("./output/B16_result.png", img_ones);
	//cvtColor(img_output, img_output, CV_BGR2GRAY);
	//cvtColor(img_480, img_480, CV_BGR2GRAY);
	//imwrite("./output/result_of_kernels.bmp", img_output);
	//imwrite("./output/result_of_480.bmp", img_480);
}

void match_new(const vector<vector<Point>>& centers_vec,
	const vector<vector<VectorXd>>& data,
	const vector<VectorXd>& kernels, RGB select_rgb)
{
	vector<int> first_point_indent(centers_vec.size(), 0);
	{
		int to_alignment = 0;
		for (int i=1;i< centers_vec.size(); i++)
		{
			if (centers_vec[i][0].x < centers_vec[to_alignment][0].x)
				to_alignment = i;
		}
		cout << "the alignment line " << to_alignment <<"[0].x = "<< centers_vec[to_alignment][0].x << endl;
		for (int vj = 0; vj < centers_vec.size(); vj+=2)
		{
			for (int vi = 0; vi < centers_vec[to_alignment].size(); vi++)
			{
				if (centers_vec[vj][0].x > centers_vec[to_alignment][vi].x+6)
				{
					/*cout << (centers_vec[vj][0].x - centers_vec[to_alignment][vi].x) * 1.0
						/ (centers_vec[vj][0].y - centers_vec[to_alignment][vi].y) << endl;*/
					first_point_indent[vj]++;
				}
				else
					break;
			}
			for (int vi = 0; vi < centers_vec[to_alignment].size(); vi++)
			{
				if (centers_vec[vj+1][0].x > centers_vec[to_alignment][vi].x+6+5)
				{
					/*cout << (centers_vec[vj][0].x - centers_vec[to_alignment][vi].x) * 1.0
						/ (centers_vec[vj][0].y - centers_vec[to_alignment][vi].y) << endl;*/
					first_point_indent[vj+1]++;
				}
				else
					break;
			}
		}
		for (int vj = 0; vj < centers_vec.size(); vj ++)
		{
			if (first_point_indent[vj] <= 2)
				first_point_indent[vj] = 0;
		}
		//for (int i = 0; i < 400; i++)
			//cout << first_point_indent[i] << endl;
		/*for (int vj = 1; vj < centers_vec.size(); vj += 2)
		{
			for (int vi = 0; vi < centers_vec[to_alignment].size(); vi++)
			{
				if (centers_vec[vj][0].x < centers_vec[to_alignment][vi].x + 14)
				{
					first_point_indent[vj]++;
				}
				else
					break;
			}
		}*/
	}

	flann::Matrix<int> indices;
	flann::Matrix<double> dists;
	query_by_kdtree(data, kernels, indices, dists);

	Mat img_480 = Mat(Size(IMG_WIDTH, IMG_HEIGHT), CV_8UC3, Scalar(0, 0, 0)),
		img_ones = Mat(Size(2436, 752*3), CV_8UC3, Scalar(0, 0, 0));
	set<int> index_set;
	VectorXd ones(9);
	ones << 1, 1, 1,
		1, 1, 1,
		1, 1, 1;
	int centers_id = 0;
	for (int vj = 0; vj < centers_vec.size(); vj++)
	{
		for (int vi = 0; vi < centers_vec[vj].size(); vi++)
		{
			if (data[vj][vi][0] < 0)
			{
				continue;
			}
			//int index = find_most_similar(data[vj][vi].normalized(), kernels_normalized);
			int index = indices[centers_id++][0];
			index_set.insert(index);
			MatrixXd A(kernels[index]);
			double value = A.colPivHouseholderQr().solve(data[vj][vi])[0];
			draw_pic_with_kernel(centers_vec[vj][vi].y, centers_vec[vj][vi].x,
				value, kernels[index], img_output);
			draw_pic_with_kernel(centers_vec[vj][vi].y, centers_vec[vj][vi].x,
				480, kernels[index], img_480);
			Vec3b color(0, 0, 0);
			color[select_rgb] = 5 * sqrt(value);
			if (select_rgb == BLUE)
			{
				img_ones.at<Vec3b>(vj * 2, 2 * (vi+ first_point_indent[vj]) + (vj + 1) % 2) = color;
			}
			else if (select_rgb == GREEN)
			{
				img_ones.at<Vec3b>(vj * 2 + 1, vi + first_point_indent[vj]) = color;
			}
			else
			{
				img_ones.at<Vec3b>(vj * 2, 2 * (vi + first_point_indent[vj]) + vj % 2) = color;
			}
		}
	}
	cout << "used kernels count: " << index_set.size() << endl;
	imwrite("./output/result_of_kernels.png", img_output);
	imwrite("./output/result_of_480_mul_kernels.png", img_480);
	imwrite("./output/B16_result.png", img_ones);
}

void pentile_rgb_relationship
(const vector<vector<Point>>& centers_vec,
	const vector<vector<VectorXd>>& data,
	const vector<VectorXd>& kernels,
	const char* outputfile)
{
	flann::Matrix<int> indices;
	flann::Matrix<double> dists;
	query_by_kdtree(data, kernels, indices, dists);
	set<int> index_set;
	//ofstream out("./output/pentile_rgb_relationship.csv");
	ofstream out(outputfile);
	//int centers_id = 0, base = 255 - centers_vec.size()/4;
	//double brightness = 0, gray = 0;
	//for (int vj = 0; vj < centers_vec.size(); vj++)
	//{
	//	double sum_brightness = 0, sum_gray = 0, valid_cnt = 0;
	//	for (int vi = 0; vi < centers_vec[vj].size(); vi++)
	//	{
	//		if (data[vj][vi][0] < 0)
	//		{
	//			centers_id++;
	//			continue;
	//		}
	//		valid_cnt++;
	//		int index = indices[centers_id++][0];
	//		index_set.insert(index);
	//		MatrixXd A(kernels[index]);
	//		double value = A.colPivHouseholderQr().solve(data[vj][vi])[0];
	//		sum_brightness += value;
	//		sum_gray += data[vj][vi][4];
	//		//out << base << "," << value << endl;
	//	}
	//	sum_brightness /= valid_cnt;
	//	sum_gray /= valid_cnt;
	//	//out << sum_gray << "," << sum_brightness << ",";
	//	if ((vj + 4 - centers_vec.size() % 4) % 4 == 3)
	//	{
	//		out << base++ << "," << (brightness + sum_brightness) / (vj < 4 ? vj + 1 : 4)
	//			<<"," << (gray + sum_gray) / (vj < 4 ? vj + 1 : 4) << endl;
	//		brightness = 0;
	//		gray = 0;
	//	}
	//	else
	//	{
	//		brightness += sum_brightness;
	//		gray += sum_gray;
	//	}
	//}

	{
		int centers_id = 0;
		double brightness = 0, gray = 0;
		for (int vj = 0; vj < centers_vec.size(); vj++)
		{
			double sum_brightness = 0, sum_gray = 0, valid_cnt = 0;
			for (int vi = 0; vi < centers_vec[vj].size(); vi++)
			{
				if (data[vj][vi][0] < 0)
				{
					centers_id++;
					continue;
				}
				valid_cnt++;
				int index = indices[centers_id++][0];
				index_set.insert(index);
				MatrixXd A(kernels[index]);
				double value = A.colPivHouseholderQr().solve(data[vj][vi])[0];
				sum_brightness += value;
				sum_gray += data[vj][vi][4];
			}
			sum_brightness /= valid_cnt;
			sum_gray /= valid_cnt;
			out << sum_brightness << "," << sum_gray << endl;
		}
	}

	out.close();
	cout << "total lines: " << centers_vec.size() << endl
		<< "used kernels count: " << index_set.size() << endl;
}

double compute_sigma(const vector<vector<VectorXd>>& data,
	vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	Performance p;
	get_kernels(kernels, sigmax, sigmay);
	//cout << kernels.size() << endl;
	vector<double> kd_tree_kernel(kernels.size() * 9, 0),
		kd_tree_data;
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels[i].normalize();
		for (int mi = 0; mi < 9; mi++)
			kd_tree_kernel[i * 9 + mi] = kernels[i][mi];
	}
	for(auto d:data)
		for (auto dd : d)
		{
			dd.normalize();
			for (int mi = 0; mi < 9; mi++)
				kd_tree_data.push_back(dd[mi]);
		}
	flann::Matrix<double> points_mat = flann::Matrix<double>(&kd_tree_kernel[0], kernels.size(), 9);
	flann::Matrix<double> query = flann::Matrix<double>(&kd_tree_data[0], kd_tree_data.size() / 9, 9);
	int nns_number = 1;
	flann::Matrix<int> indices(new int[query.rows*nns_number], query.rows, nns_number);
	flann::Matrix<double> dists(new double[query.rows*nns_number], query.rows, nns_number);
	flann::Index<flann::L2<double>> index(points_mat, flann::KDTreeIndexParams(4));
	index.buildIndex();
	index.knnSearch(query, indices, dists, nns_number, flann::SearchParams(128));
	double loss = 0;
	//cout << dists.rows << endl;
	for (int i = 0; i < dists.rows; i++) {
		loss += dists[i][0];
	}
	
	/*for (int i = 0; i < kernels.size(); i++)
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
	}*/
	cout << "compute loss: " << p.end() << "s" << endl;
	cout << "sigma: " << sigmax << "," << sigmay << ", loss: " << loss << endl << endl;
	return loss;
}

void compute_dumura(vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error, int xy, 
	// change sigmax:0, change sigmax y:1, the same sigma:2, no sigma compute:3
	// compute relationship between rgb & pentile-4
	double from, double to, double another, double add, const char* prefix, RGB select_rgb,
	const char* outputfile)
{
	Performance p;
	img_output = Mat(Size(IMG_WIDTH, IMG_HEIGHT), CV_8UC3, Scalar(0, 0, 0));
	for (auto ce : centers_error)
	{
		//draw_pic_with_scalar(ce.y, ce.x, Vec3b(0, 0, 255), img_output);
		img_output.at<Vec3b>(ce.y, ce.x) = Vec3b(0, 0, 255);
	}
	vector<VectorXd> kernels;
	if (xy == 0 || xy == 1)
		kernels.resize(sample_of_center*sample_of_center * 4, VectorXd(9));
	else
		kernels.resize(sample_of_center*sample_of_center, VectorXd(9));

	/*ofstream out("./output/B16_loss.csv");
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
	
	char outfile[MAX_PATH];
	double sigmax, sigmay,
		loss = 1000000, loss_old, flag = -1;
	if (xy == 0)
	{
		cout << " change sigmax" << endl;
		sigmax = from;
		sigmay = another;
		sprintf_s(outfile, "./output/%s_%.2f,%.2f__%.2f_c%d_g%d.csv",
			prefix, sigmax, to, sigmay,
			sample_of_center, sample_of_guass_point);
		ofstream out(outfile);
		for (; sigmax < to; sigmax += add)
		{
			loss_old = loss;
			loss = compute_sigma(data, kernels, sigmax, sigmay);
			out << sigmax << "," << sigmay << "," << loss << endl;
		}
		cout << "sigma:" << sigmax << "," << sigmay << endl;
		out.close();
	}
	else if(xy==1)
	{
		cout << " change sigmay" << endl;
		sigmax = another;
		sigmay = from;
		loss = 1000000;
		flag = -1;
		sprintf_s(outfile, "./output/%s_%.2f__%.2f,%.2f_c%d_g%d.csv",
			prefix, sigmax, sigmay, to,
			sample_of_center, sample_of_guass_point);
		ofstream out(outfile);
		for (; sigmay < to; sigmay += add)
		{
			loss_old = loss;
			loss = compute_sigma(data, kernels, sigmax, sigmay);
			out << sigmax << "," << sigmay << "," << loss << endl;
		}
		cout << "sigma:" << sigmax << "," << sigmay << endl;
		out.close();
	}
	else if (xy == 2)
	{
		cout << " the same sigmax and sigmay" << endl;
		sigmax = from;
		loss = 1000000;
		flag = -1;
		sprintf_s(outfile, "./output/%s_iso_c%d_g%d.csv",
			prefix, sample_of_center, sample_of_guass_point);
		cout << endl << "write to " << outfile << endl;
		ofstream out(outfile);
		for (; sigmax < to; sigmax += add)
		{
			loss_old = loss;
			loss = compute_sigma(data, kernels, sigmax, sigmax);
			out << sigmax << "," << sigmax << "," << loss << endl;
		}
		cout << "sigma:" << sigmax << "," << sigmax << endl;
		out.close();
	}
	else if(xy == 3)
	{
		cout << " match guassian" << endl;
		// blue 1.06 0.92
		// green 0.71, 0.85
		double sigmax = 0.77, sigmay = 0.77;
		cout << "sigma: " << sigmax << "," << sigmay << endl;
		get_kernels(kernels, sigmax, sigmay);
		cout <<"total kernels count:"<< kernels.size() << endl;
		combination_kernels(kernels);
		match_new(centers_vec, data, kernels, select_rgb);
	}
	else if (xy == 4)
	{
		cout << " compute relationship between rgb & pentile" << endl;
		double sigmax = 0.59, sigmay = 0.59;
		cout << "sigma: " << sigmax << "," << sigmay << endl;
		get_kernels(kernels, sigmax, sigmay);
		cout << "total kernels count:" << kernels.size() << endl;
		cout<<endl<<"write to "<< outputfile << endl;
		pentile_rgb_relationship(centers_vec, data, kernels, outputfile);
	}
	else
	{
		cout << "select: \
			0-change sigmax, \
			1-change sigmax y, \
			2-the same sigma, 3-no sigma compute\
			4-compute relationship between rgb & pentile" << endl;
	}
	p.endAndPrint();
}