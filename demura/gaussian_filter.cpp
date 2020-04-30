#include "gaussian_filter.h"
#include "pentile.h"
#include <Eigen/Dense>
#include <fstream>
#include <flann/flann.hpp>
#include <flann/util/matrix.h>
#include <unordered_map>
#include <numeric>

const double PI = 4.0*atan(1.0);
using namespace std;
using cv::Mat;
using cv::Vec3b;
using cv::Point;
using cv::Size;
using cv::Scalar;
using namespace Eigen;

const int sample_of_center = 64, base = sample_of_center * sample_of_center;
const double center_interval = 1.0 / sample_of_center;
const int sample_of_guass_point = 64;
const double point_interval = 1.0 / sample_of_guass_point;

void output(const VectorXd& a)
{
	cout << "[";
	for (int i = 0; i < 8; i++)
		cout << a[i] << ",";
	cout << a[8] << "]" << endl;
}

double gauss_value(const double x, const double y,
	const double cx, const double cy, const double sigmax, const double sigmay)
{
	/*return (1 / (2 * PI*sigma*sigma)) * exp(-((x - cx)*(x - cx)
		+ (y - cy)*(y - cy)) / (2 * sigma*sigma));*/
	return (1 / (2 * PI*sigmax*sigmay)) * exp(-
		((x - cx)*(x - cx) / (2 * sigmax * sigmax) + (y - cy)*(y - cy) / (2 * sigmay * sigmay)));
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

void rotate_xy(double& x, double& y, const double rx, const double ry, const double theta)
{
	double c = cos(theta*PI / 180), s = sin(theta*PI / 180),
		xx = x - rx, yy = y - ry;
	x = c * xx - s * yy + rx;
	y = s * xx + c * yy + ry;
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

void get_kernels(vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	double cx_begin = 1, cy_begin = 1;
	int index = 0;
	if (sigmax == sigmay)
	{
		for (double yi = center_interval/2; yi <= 1; yi += center_interval)
		{
			for (double xi = center_interval/2; xi <= 1; xi += center_interval)
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

void query_by_kdtree(vector<double>& kd_tree_data, const vector<VectorXd>& kernels,
	flann::Matrix<int>& indices, flann::Matrix<double>& dists)
{
	vector<double> kd_tree_kernel(kernels.size() * 9, 0);
		//kd_tree_data;
	for (int i = 0; i < kernels.size(); i++)
	{
		VectorXd tmp = kernels[i];
		tmp.normalize();
		for (int mi = 0; mi < 9; mi++)
			kd_tree_kernel[i * 9 + mi] = tmp[mi];
	}
	/*for (auto d : data)
		for (auto dd : d)
		{
			dd.normalize();
			for (int mi = 0; mi < 9; mi++)
				kd_tree_data.push_back(dd[mi]);
		}*/
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

//#define SQRT
void match_with_location(Mat& result, 
	const char *outfile, const char *result_txt, const char *result_csv,
	const vector<vector<LED_info>>& relationship,
	const vector<int>& capture_pentile_g_value,
	const vector<int>& expo,
	const vector<Mat>& pic, const int primary_pic,
	const vector<VectorXd>& kernels, RGB select_rgb)
{
	try {
		flann::Matrix<int> indices;
		flann::Matrix<double> dists;
		vector<double> kd_tree_data;
		VectorXd center_region(9);
		for (auto d : relationship)
			for (auto dd : d)
			{
				center_region << pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x - 1),
					pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x),
					pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x + 1),
					pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x - 1),
					pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x),
					pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x + 1),
					pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x - 1),
					pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x),
					pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x + 1);
				center_region.normalize();
				for (int mi = 0; mi < 9; mi++)
					kd_tree_data.push_back(center_region[mi]);
			}
		query_by_kdtree(kd_tree_data, kernels, indices, dists);

		vector<vector<double>> guassian_mul_value(pic.size(), vector<double>(indices.rows, 0));
		set<int> index_set;
		VectorXd ones(9);
		ones << 1, 1, 1, 1, 1, 1, 1, 1, 1;
		int centers_id = 0;
		Point cur;
		for (int vj = 0; vj < relationship.size(); vj++)
		{
			for (int vi = 0; vi < relationship[vj].size(); vi++)
			{
				int index = indices[centers_id][0];
				MatrixXd A(kernels[index]);
				if (relationship[vj][vi].state == VALID)
				{
					for (int pici = 0; pici < pic.size(); pici++)
					{
						cur = relationship[vj][vi].pixel;
						center_region << pic[pici].at<byte>(cur.y - 1, cur.x - 1),
							pic[pici].at<byte>(cur.y - 1, cur.x),
							pic[pici].at<byte>(cur.y - 1, cur.x + 1),
							pic[pici].at<byte>(cur.y, cur.x - 1),
							pic[pici].at<byte>(cur.y, cur.x),
							pic[pici].at<byte>(cur.y, cur.x + 1),
							pic[pici].at<byte>(cur.y + 1, cur.x - 1),
							pic[pici].at<byte>(cur.y + 1, cur.x),
							pic[pici].at<byte>(cur.y + 1, cur.x + 1);
						if (center_region.sum() < 0.00001)
							guassian_mul_value[pici][centers_id] = 0;
						else
						{
#ifdef SQRT
							guassian_mul_value[pici][centers_id] = sqrt(A.colPivHouseholderQr().solve(center_region)[0]) / expo[pici];
#else
							guassian_mul_value[pici][centers_id] = A.colPivHouseholderQr().solve(center_region)[0] / expo[pici];
#endif
						}
					}
				}
				index_set.insert(index);
				centers_id++;
			}
		}
		//int primary_pic = 4;
		// 计算当前gaussian val的mean
		double mean;
		{
			int cnt = 0;
			double sum = 0.0, accum = 0.0,
				mmax = guassian_mul_value[primary_pic][0],
				mmin = INT_MAX;
			for (auto d : guassian_mul_value[primary_pic])
			{
				if (d > 0.001)
				{
					sum += d;
					cnt++;
					if (mmin > d) mmin = d;
					if (mmax < d) mmax = d;
				}
			}
			mean = sum / cnt;
			for (auto d : guassian_mul_value[primary_pic])
			{
				if (d > 0.001)
				{
					accum += (d - mean)*(d - mean);
				}
			}
			double stdev = sqrt(accum / cnt); //方差

			cout << "used kernels count: " << index_set.size() << endl;
			cout << "sum = " << sum << ", size = " << guassian_mul_value[primary_pic].size() << endl
				<< "min = " << mmin << ", max = " << mmax << endl
				<< "mean = " << mean << ", stdev = " << stdev << endl;
		}
		//vector<int> capture_pentile_g_value({ 8,12,16,20,24,28,32 });
		centers_id = 0;
		vector<double> mura_value;
		long long cnt_small = 0, cnt_large = 0;
		cout << "compute gaussian" << endl;
		for (int vj = 0; vj < relationship.size(); vj++)
		{
			for (int vi = 0; vi < relationship[vj].size(); vi++)
				//for(auto re: relationship[vj])
			{
				//cout << guassian_mul_value[centers_id] << endl;
				if (relationship[vj][vi].locate.y >= 0 && relationship[vj][vi].locate.y < result.rows
					&& relationship[vj][vi].locate.x >= 0 && relationship[vj][vi].locate.x < result.cols)
				{
					if (relationship[vj][vi].state == INVALID)
						result.at<cv::Vec3b>(relationship[vj][vi].locate.y,
							relationship[vj][vi].locate.x)[select_rgb] = 16;
					else if (guassian_mul_value[0][centers_id] >= mean)
					{
						result.at<cv::Vec3b>(relationship[vj][vi].locate.y,
							relationship[vj][vi].locate.x)[select_rgb] = capture_pentile_g_value[0];
						cnt_small++;
					}
					else if (guassian_mul_value[guassian_mul_value.size() - 1][centers_id] <= mean)
					{
						cnt_large++;
						result.at<cv::Vec3b>(relationship[vj][vi].locate.y,
							relationship[vj][vi].locate.x)[select_rgb] = *capture_pentile_g_value.rbegin();
					}
					else
					{
						int l = 0, r = pic.size() - 1, mid;
						while (l < r)
						{
							mid = (l + r) >> 1;
							if (guassian_mul_value[mid][centers_id] == mean)
							{
								r = mid;
								break;
							}
							else if (guassian_mul_value[mid][centers_id] > mean)
								r = mid;
							else
								l = mid + 1;
						}
						/*cout << "binary search end ..." << endl;
						cout <<"l="<<l<<", r="<<r<<","
							<< capture_pentile_g_value[r - 1] << ", "
							<< guassian_mul_value[r-1][centers_id] << ", "
							<< capture_pentile_g_value[r] << ","
							<< guassian_mul_value[r][centers_id] << ","
							<< (mean*(capture_pentile_g_value[r] - capture_pentile_g_value[r - 1]) +
								capture_pentile_g_value[r - 1] * guassian_mul_value[r][centers_id] -
								capture_pentile_g_value[r] * guassian_mul_value[r - 1][centers_id]) /
								(guassian_mul_value[r][centers_id] - guassian_mul_value[r - 1][centers_id]) << ","
							<< mean
							<< endl;*/
						/*double tmp = (mean*(capture_pentile_g_value[r] - capture_pentile_g_value[r - 1]) +
							capture_pentile_g_value[r - 1] * guassian_mul_value[r][centers_id] -
							capture_pentile_g_value[r] * guassian_mul_value[r - 1][centers_id]) /
							(guassian_mul_value[r][centers_id] - guassian_mul_value[r - 1][centers_id]);*/
						double tmp;
						if (capture_pentile_g_value[r] - mean < mean - capture_pentile_g_value[r - 1])
							tmp = capture_pentile_g_value[r];
						else
							tmp = capture_pentile_g_value[r-1];
						result.at<cv::Vec3b>(relationship[vj][vi].locate.y,
							relationship[vj][vi].locate.x)[select_rgb] = tmp;
						mura_value.push_back(tmp);
					}
				}
				centers_id++;
			}
		}

		cout << "compute invalid" << endl;
		vector<int> dx({ -1,0,1, -1,1, -1,0,1 });
		vector<int> dy({ -1,-1,-1, 0,0, 1,1,1 });
		int cnt;
		for (int vj = 0; vj < relationship.size(); vj++)
		{
			for (auto r : relationship[vj])
			{
				if (r.locate.x >= 0 && r.locate.y >= 0 && r.state == INVALID)
				{
					cnt = 0;
					result.at<cv::Vec3b>(r.locate.y, r.locate.x)[select_rgb] = 0;
					for (int i = 0; i < dx.size(); i++)
					{
						if (r.locate.y + dy[i] >= 0 && r.locate.y + dy[i] < result.rows
							&&  r.locate.x + dx[i] >= 0 && r.locate.x + dx[i] < result.cols)
						{
							result.at<cv::Vec3b>(r.locate.y, r.locate.x)[select_rgb] +=
								result.at<cv::Vec3b>(r.locate.y + dy[i], r.locate.x + dx[i])[select_rgb];
							cnt++;
						}
					}
					if (cnt) result.at<cv::Vec3b>(r.locate.y, r.locate.x)[select_rgb] /= cnt;
				}
			}
		}

		cout << "output gaussian, mean, mura" << endl;
		ofstream out;
		out.open(result_txt);
		{
			centers_id = 0;
			for (int vj = 0; vj < relationship.size(); vj++)
			{
				for (int vi = 0; vi < relationship[vj].size(); vi++)
				{
					out << relationship[vj][vi].locate.y << " " << relationship[vj][vi].locate.x << " ";
					if (relationship[vj][vi].locate.x >= 0 && relationship[vj][vi].locate.y >= 0 
						&& relationship[vj][vi].state == VALID)
					{
						for (int r = 0; r < guassian_mul_value.size(); r++)
						{
							out << guassian_mul_value[r][centers_id] << " ";
						}
						out << mean << " " << (int)result.at<cv::Vec3b>(
							relationship[vj][vi].locate.y, relationship[vj][vi].locate.x)[select_rgb];
						centers_id++;
					}
					out << endl;
				}
			}
		}
		out.close();

		cout << "output mura value distribution" << endl;
		out.open(result_csv);
		{
			double sum = 0.0, accum = 0.0,
				mmax = mura_value[0],
				mmin = INT_MAX;
			int cnt = 0;
			for (auto d : mura_value)
			{
				if (d > 0.001)
				{
					sum += d;
					cnt++;
					if (mmin > d) mmin = d;
					if (mmax < d) mmax = d;
				}
			}
			mean = sum / cnt;
			for (auto d : mura_value)
			{
				if (d > 0.001)
				{
					accum += (d - mean)*(d - mean);
				}
			}
			double stdev = sqrt(accum / cnt); //方差

			cout << endl
				<< "----------------------------" << endl
				<< "sum = " << sum << ", size = " << centers_id << endl
				<< "min = " << mmin << ", max = " << mmax << endl
				<< "mean = " << mean << ", stdev = " << stdev << endl
				<< "small cnt: " << cnt_small << endl
				<< "large cnt: " << cnt_large << endl;
			out << "sum,mean,min,max,stdev" << endl
				<< sum << "," << mean << "," << mmin << "," << mmax << "," << stdev << endl;
			for (int i = 0; i < mura_value.size(); i++)
			{
				out << mura_value[i] << ",";
				if ((i + 1) % 100 == 0)
					out << endl;
			}
		}
		out.close();
		//imwrite(outfile, result);
		cout << "end compute." << endl;
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}
}

double compute_sigma(vector<double>& kd_tree_data,
	vector<VectorXd>& kernels, const double sigmax, const double sigmay)
{
	Performance p;
	get_kernels(kernels, sigmax, sigmay);
	//cout << kernels.size() << endl;
	vector<double> kd_tree_kernel(kernels.size() * 9, 0);
	for (int i = 0; i < kernels.size(); i++)
	{
		kernels[i].normalize();
		for (int mi = 0; mi < 9; mi++)
			kd_tree_kernel[i * 9 + mi] = kernels[i][mi];
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

extern void draw_cross(Mat& img, int x, int y, Vec3b color);
void compute_dumura(const vector<vector<vector<LED_info>>>& relationship,
	const vector<int>& capture_pentile_g_value,
	const vector<int>& expo,
	const vector<cv::Mat>& pic, const int primary_pic,
	const char* output_prefix, int width, int height)
{
	Performance p;
	//img_output = Mat(Size(IMG_WIDTH, IMG_HEIGHT), CV_8UC3, Scalar(0, 0, 0));
	vector<VectorXd> kernels;
	kernels.resize(sample_of_center*sample_of_center, VectorXd(9));
	cout << "[match guassian] matching..." << endl;
	// blue 1.06 0.92
	// green 0.71, 0.85
	//double sigmax = 0.77, sigmay = 0.77; 
	Mat img_result = Mat(Size(3000, 2000), CV_8UC3, Scalar(0, 0, 0));
	vector<double> sigma_bgr({ 0.66, 0.62, 0.64 });
	//for (int i = 0; i < relationship.size(); i++)
	int select_rgb = 1;
	{
		/*{
			vector<double> kd_tree_data;
			VectorXd center_region(9);
			for (auto d : relationship[select_rgb])
				for (auto dd : d)
				{
					center_region << pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x - 1),
						pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x),
						pic[primary_pic].at<byte>(dd.pixel.y - 1, dd.pixel.x + 1),
						pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x - 1),
						pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x),
						pic[primary_pic].at<byte>(dd.pixel.y, dd.pixel.x + 1),
						pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x - 1),
						pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x),
						pic[primary_pic].at<byte>(dd.pixel.y + 1, dd.pixel.x + 1);
					center_region.normalize();
					for (int mi = 0; mi < 9; mi++)
						kd_tree_data.push_back(center_region[mi]);
				}
			double sigma_init, sigma, loss_old, loss = 1000000;
			bool flag = false;
			for (sigma_init = 0.5; sigma_init < 1.5; sigma_init += 0.1)
			{
				loss_old = loss;
				loss = compute_sigma(kd_tree_data, kernels, sigma_init, sigma_init);
				if (!flag && loss_old - loss > 0)
				{
					flag = true;
				}
				else if (flag && loss_old - loss < 0)
				{
					break;
				}
			}
			sigma_init -= 0.1;
			cout << "init sigma: " << sigma_init << " ";
			flag = false;
			for (sigma = sigma_init - 0.1; sigma < sigma_init + 0.1; sigma += 0.01)
			{
				loss_old = loss;
				loss = compute_sigma(kd_tree_data, kernels, sigma, sigma);
				if (!flag && loss_old - loss > 0)
				{
					flag = true;
				}
				else if (flag && loss_old - loss < 0)
				{
					break;
				}
			}
			sigma -= 0.01;
			sigma_bgr[select_rgb] = sigma;
		}*/

		cout << "use sigma: " << sigma_bgr[select_rgb] << endl;
		get_kernels(kernels, sigma_bgr[select_rgb], sigma_bgr[select_rgb]);
		char result_file[MAX_PATH], result_txt[MAX_PATH], result_csv[MAX_PATH];
		sprintf(result_file, "%s/result_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		sprintf(result_txt, "%s/gaussian_val.txt", output_prefix);
		sprintf(result_csv, "%s/distribution.csv", output_prefix);
		match_with_location(img_result, result_file, result_txt, result_csv,
			relationship[select_rgb], capture_pentile_g_value, expo,
			pic, primary_pic, kernels, (RGB)select_rgb);
	}
	p.endAndPrint();

	cout << "total kernels count:" << kernels.size() << endl;
	char result_file[MAX_PATH];
	sprintf(result_file, "%s/result.png", output_prefix);
	img_result = img_result(cv::Rect(0, 0, width, height));
	imwrite(result_file, img_result);
	Scalar scal(0, 0, 0);
	Mat rgb(Size(height, width), CV_8UC3, scal);
	
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			rgb.at<Vec3b>(width - 1 - x, y)[1] =
				img_result.at<Vec3b>(y, x)[1] == 0 ? 32 : (int)img_result.at<Vec3b>(y, x)[1];
		}
	}
	vector<int> need_x({ 60,1060 }), need_y({ 130, width - 130 });
	for (auto y : need_y)
	{
		for (auto x : need_x)
		{
			draw_cross(rgb, x, y, Vec3b(0, 32, 0));
		}
	}
	sprintf(result_file, "%s/result_rotate_90.bmp", output_prefix);
	imwrite(result_file, rgb);
	Mat pentile;
	rgb2pentile(rgb, pentile);
	sprintf(result_file, "%s/result_pentile.bmp", output_prefix);
	imwrite(result_file, pentile);
}

void match_with_location2(Mat& result, const char *outfile, 
	const char *input_txt_file, const char *output_csv_file,
	const vector<vector<LED_info>>& relationship,
	const vector<int>& capture_pentile_g_value,
	const Mat& img, const vector<VectorXd>& kernels)//, RGB select_rgb)
{
	flann::Matrix<int> indices;
	flann::Matrix<double> dists;
	vector<double> kd_tree_data;
	VectorXd center_region(9);
	Mat pic = img.clone();
	cvtColor(pic, pic, cv::COLOR_BGR2GRAY);
	for (auto d : relationship)
		for (auto dd : d)
		{
			center_region << pic.at<byte>(dd.pixel.y - 1, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y - 1, dd.pixel.x),
				pic.at<byte>(dd.pixel.y - 1, dd.pixel.x + 1),
				pic.at<byte>(dd.pixel.y, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y, dd.pixel.x),
				pic.at<byte>(dd.pixel.y, dd.pixel.x + 1),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x + 1);
			center_region.normalize();
			for (int mi = 0; mi < 9; mi++)
				kd_tree_data.push_back(center_region[mi]);
		}
	query_by_kdtree(kd_tree_data, kernels, indices, dists);

	vector<double> guassian_mul_value(indices.rows, -1);
	set<int> index_set;
	VectorXd ones(9);
	ones << 1, 1, 1, 1, 1, 1, 1, 1, 1;
	int centers_id = 0;
	Point cur;
	for (int vj = 0; vj < relationship.size(); vj++)
	{
		for (int vi = 0; vi < relationship[vj].size(); vi++)
		{
			int index = indices[centers_id][0];
			MatrixXd A(kernels[index]);
			if (relationship[vj][vi].state == VALID)
			{
				cur = relationship[vj][vi].pixel;
				center_region << pic.at<byte>(cur.y - 1, cur.x - 1),
					pic.at<byte>(cur.y - 1, cur.x),
					pic.at<byte>(cur.y - 1, cur.x + 1),
					pic.at<byte>(cur.y, cur.x - 1),
					pic.at<byte>(cur.y, cur.x),
					pic.at<byte>(cur.y, cur.x + 1),
					pic.at<byte>(cur.y + 1, cur.x - 1),
					pic.at<byte>(cur.y + 1, cur.x),
					pic.at<byte>(cur.y + 1, cur.x + 1);
				guassian_mul_value[centers_id] = A.colPivHouseholderQr().solve(center_region)[0];
			}
			index_set.insert(index);
			centers_id++;
		}
	}

	double mean;
	{
		int cnt = 0;
		double sum = 0.0, accum = 0.0,
			mmax = guassian_mul_value[0],
			mmin = INT_MAX;
		for (auto d : guassian_mul_value)
		{
			if (d > 0.001)
			{
				sum += d;
				cnt++;
				if (mmin > d) mmin = d;
				if (mmax < d) mmax = d;
			}
		}
		mean = sum / cnt;
		for (auto d : guassian_mul_value)
		{
			if (d > 0.001)
			{
				accum += (d - mean)*(d - mean);
			}
		}
		double stdev = sqrt(accum / cnt); //方差

		cout << "used kernels count: " << index_set.size() << endl;
		cout << "sum = " << sum << ", size = " << guassian_mul_value.size() << endl
			<< "min = " << mmin << ", max = " << mmax << endl
			<< "mean = " << mean << ", stdev = " << stdev << endl;
	}
	struct pair_hash
	{
		inline size_t operator()(const pair<int, int>& a) const {
			return a.first * 12'000 + a.second;
		}
	};
	unordered_map<pair<int, int>, vector<double>, pair_hash> mmap(300'000);
	{
		int size = capture_pentile_g_value.size();
		vector<double> val(size);
		pair<int, int> tmp;
		ifstream in(input_txt_file);
		if (!in) 
		{
			cout << "open fail: " << input_txt_file << endl;
			return ;
		}
		while (!in.eof())
		{
			in >> tmp.first >> tmp.second;
			//cout << d1 <<","<< d2 << endl;
			for (int i = 0; i < size; i++)
			{
				in >> val[i];
			}
			mmap[tmp] = val;
		}
		cout << "mmap size: " << mmap.size() << endl;
		cout << "mmap[0, 129] " << mmap[{0, 129}][0] << " " << mmap[{0, 129}][1] << " " << mmap[{0, 129}][2] << " " << endl;
	}

	ofstream out(output_csv_file); //"./output/1.csv"
	centers_id = 0;
	vector<double> mura_value;
	pair<int, int> tmp;
	out << "locate_x, locate_y, pixel_x, pixel_y, gaussian_val, mura_val" << endl;
	for (int vj = 0; vj < relationship.size(); vj++)
	{
		for (auto re : relationship[vj])
		{
			/*tmp.first = re.locate.y;
			tmp.second = re.locate.x;
			if (re.locate.y >= 0 && re.locate.x >= 0)
				if (re.state == INVALID)
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(0, 0, 0);
				else if (mmap.find(tmp) == mmap.end())
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(255, 255, 255);
				else if (guassian_mul_value[centers_id] < mmap[tmp][0])
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(255, 0, 0);
				else if (guassian_mul_value[centers_id] < mmap[tmp][1])
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(0, 255, 0);
				else if (guassian_mul_value[centers_id] < mmap[tmp][2])
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(0, 0, 255);
				else
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(0, 255, 255);*/
			tmp.first = re.locate.y;
			tmp.second = re.locate.x;
			if (re.state == VALID && re.locate.y >= 0 && re.locate.x >= 0)
			{
				if (mmap.find(tmp) == mmap.end())
					result.at<cv::Vec3b>(re.locate.y, re.locate.x) = Vec3b(255, 255, 255);
				else if (mmap[tmp][0] > guassian_mul_value[centers_id])
				{
					out << re.locate.x << "," << re.locate.y << ","
						<< re.pixel.x << "," << re.pixel.y << ", " << -1 << endl;
				}
				else if (mmap[tmp][capture_pentile_g_value.size() - 1]
					< guassian_mul_value[centers_id])
				{
					out << re.locate.x << "," << re.locate.y << ","
						<< re.pixel.x << "," << re.pixel.y << ", " << 300 << endl;
				}
				else
				{
					int l = 0, r = capture_pentile_g_value.size() - 1, mid;
					while (l < r)
					{
						mid = (l + r) >> 1;
						if (mmap[tmp][mid] == guassian_mul_value[centers_id])
						{
							r = mid;
							break;
						}
						else if (mmap[tmp][mid] > guassian_mul_value[centers_id])
							r = mid;
						else
							l = mid + 1;
					}
					result.at<cv::Vec3b>(re.locate.y, re.locate.x)[1] =
						(guassian_mul_value[centers_id] *
						(capture_pentile_g_value[r] - capture_pentile_g_value[r - 1]) +
							capture_pentile_g_value[r - 1] * mmap[tmp][r] -
							capture_pentile_g_value[r] * mmap[tmp][r - 1])
						/ (mmap[tmp][r] - mmap[tmp][r - 1]);
					//cout << result.at<cv::Vec3b>(re.locate.y, re.locate.x)[1] << endl;
				}
				out << re.locate.x << "," << re.locate.y << "," << re.pixel.x << "," << re.pixel.y
					<< ", " << guassian_mul_value[centers_id] << ","
					<< (int)result.at<cv::Vec3b>(re.locate.y, re.locate.x)[1] << endl;
			}
			else
			{
				out << re.locate.x << "," << re.locate.y << "," << re.pixel.x << "," << re.pixel.y
					<< ", " << guassian_mul_value[centers_id] << endl;
			}
			centers_id++;
		}
	}
	out.close();
}

void compute_dumura_single_pic(vector<vector<vector<LED_info>>>& relationship,
	const vector<int>& capture_pentile_g_value,	const cv::Mat& img, 
	const char* input_prefix, const char* output_prefix, int width, int height)
{
	Performance p;
	vector<VectorXd> kernels;
	kernels.resize(sample_of_center*sample_of_center, VectorXd(9));
	cout << "[match guassian] matching single..." << endl;
	// blue 1.06 0.92
	// green 0.71, 0.85
	//double sigmax = 0.77, sigmay = 0.77; 
	Mat img_result = Mat(Size(3000, 2000), CV_8UC3, Scalar(0, 0, 0));
	vector<double> sigma_bgr({ 0.66, 0.55, 0.64 });
	//for (int i = 0; i < relationship.size(); i++)
	int select_rgb = 1;
	{
		/*double sigma_init, sigma, loss_old, loss = 1000000;
		bool flag = false;
		for (sigma_init=0.5; sigma_init < 1.5; sigma_init += 0.1)
		{
			loss_old = loss;
			loss = compute_sigma(data[i], kernels, sigma_init, sigma_init);
			if (!flag && loss_old - loss>0 )
			{
				flag = true;
			}
			else if(flag && loss_old - loss < 0)
			{
				break;
			}
		}
		sigma_init -= 0.1;
		cout << "init sigma: " << sigma_init << " ";
		flag = false;
		for (sigma = sigma_init - 0.1; sigma < sigma_init+0.1; sigma += 0.01)
		{
			loss_old = loss;
			loss = compute_sigma(data[i], kernels, sigma, sigma);
			if (!flag && loss_old - loss > 0)
			{
				flag = true;
			}
			else if (flag && loss_old - loss < 0)
			{
				break;
			}
		}
		sigma -= 0.01;
		cout <<"final sigma: "<< sigma << endl;*/

		get_kernels(kernels, sigma_bgr[select_rgb], sigma_bgr[select_rgb]);
		char result_file[MAX_PATH], input_txt[MAX_PATH], output_csv[MAX_PATH];
		sprintf(result_file, "%s/valid_result_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		sprintf(input_txt, "%s/gaussian_val.txt", input_prefix);
		sprintf(output_csv, "%s/valid_gaussian_val.txt", output_prefix);
		match_with_location2(img_result, result_file, input_txt, output_csv,
			relationship[select_rgb], capture_pentile_g_value,
			img, kernels); // , (RGB)select_rgb);
	}
	p.endAndPrint();

	cout << "total kernels count:" << kernels.size() << endl;
	char result_file[MAX_PATH];
	sprintf(result_file, "%s/valid_result.png", output_prefix);
	// rotate
	img_result = img_result(cv::Rect(0, 0, width, height));
	Mat rgb(Size(height, width), CV_8UC3, Scalar(0, 0, 0));
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			rgb.at<Vec3b>(width - 1 - x, y) = img_result.at<Vec3b>(y, x);
		}
	}
	imwrite(result_file, rgb);

	//Scalar scal(0, 0, 0);
	//Mat rgb(Size(height, width), CV_8UC3, scal);
	//for (int y = 0; y < height; y++)
	//{
	//	for (int x = 0; x < width; x++)
	//	{
	//		//rgb.at<Vec3b>(width - 1 - x, y)[0] =
	//			//img_result.at<Vec3b>(y, x)[0] == 0 ? 255 : (int)img_result.at<Vec3b>(y, x)[0];
	//		rgb.at<Vec3b>(width - 1 - x, y)[1] =
	//			img_result.at<Vec3b>(y, x)[1] == 0 ? 32 : (int)img_result.at<Vec3b>(y, x)[1];
	//		//rgb.at<Vec3b>(width - 1 - x, y)[2] = 
	//			//img_result.at<Vec3b>(y, x)[2] == 0 ? 255 : (int)img_result.at<Vec3b>(y, x)[2];
	//		/*rgb.at<Vec3b>(width - 1 - x, y) =
	//			color[img_result.at<Vec3b>(y, x)[1]];*/
	//	}
	//}
	//vector<int> need_x({ 60,1060 }), need_y({ 130, width - 130 });
	//for (auto y : need_y)
	//{
	//	for (auto x : need_x)
	//	{
	//		draw_cross(rgb, x, y, Vec3b(0, 32, 0));
	//	}
	//}

	//sprintf(result_file, "%s/result_rotate_90.bmp", output_prefix);
	//imwrite(result_file, rgb);
	//Mat pentile;
	//rgb2pentile(rgb, pentile);
	//sprintf(result_file, "%s/result_pentile.bmp", output_prefix);
	//imwrite(result_file, pentile);
}

void match_with_location_to_exr(vector<float>& result, const char *outfile,
	const vector<vector<LED_info>>& relationship,
	const Mat& img, const vector<VectorXd>& kernels, int width, int height)
{
	cout << "compute : " << outfile << endl;
	// flann
	Mat pic = img.clone();
	VectorXd center_region(9);
	flann::Matrix<int> indices;
	flann::Matrix<double> dists;
	vector<double> kd_tree_data;
	cvtColor(pic, pic, cv::COLOR_BGR2GRAY);
	for (auto d : relationship)
		for (auto dd : d)
		{
			center_region << pic.at<byte>(dd.pixel.y - 1, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y - 1, dd.pixel.x),
				pic.at<byte>(dd.pixel.y - 1, dd.pixel.x + 1),
				pic.at<byte>(dd.pixel.y, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y, dd.pixel.x),
				pic.at<byte>(dd.pixel.y, dd.pixel.x + 1),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x - 1),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x),
				pic.at<byte>(dd.pixel.y + 1, dd.pixel.x + 1);
			center_region.normalize();
			for (int mi = 0; mi < 9; mi++)
				kd_tree_data.push_back(center_region[mi]);
		}
	query_by_kdtree(kd_tree_data, kernels, indices, dists);

	set<int> index_set;
	VectorXd ones(9);
	ones << 1, 1, 1, 1, 1, 1, 1, 1, 1;
	int centers_id = 0;
	Point cur;
	for (int vj = 0; vj < relationship.size(); vj++)
	{
		//for (int vi = 0; vi < relationship[vj].size(); vi++)
		for(auto re: relationship[vj])
		{
			int index = indices[centers_id][0];
			MatrixXd A(kernels[index]);
			if (re.state == VALID && re.locate.x >= 0 && re.locate.x < width
				&& re.locate.y >= 0 && re.locate.y < height)
			{
				cur = re.pixel;
				center_region << pic.at<byte>(cur.y - 1, cur.x - 1),
					pic.at<byte>(cur.y - 1, cur.x),
					pic.at<byte>(cur.y - 1, cur.x + 1),
					pic.at<byte>(cur.y, cur.x - 1),
					pic.at<byte>(cur.y, cur.x),
					pic.at<byte>(cur.y, cur.x + 1),
					pic.at<byte>(cur.y + 1, cur.x - 1),
					pic.at<byte>(cur.y + 1, cur.x),
					pic.at<byte>(cur.y + 1, cur.x + 1);
				//guassian_mul_value[centers_id] = A.colPivHouseholderQr().solve(center_region)[0];
				//width - 1 - x, y  <=> y, x
				//result[((width - 1 - re.locate.x)*width + re.locate.y) * 3 + 1] = A.colPivHouseholderQr().solve(center_region)[0];
				result[(re.locate.y*width + re.locate.x) * 3 + 1] = A.colPivHouseholderQr().solve(center_region)[0] / 1000;
			}
			index_set.insert(index);
			centers_id++;
		}
	}
	cout << outfile << endl;
	save_exr_with_float(outfile, width, height, result);
}

void compute_intensity_to_exr(vector<vector<vector<LED_info>>>& relationship,
	const cv::Mat& img, const char* input_prefix, const char* output_prefix, int width, int height)
{
	Performance p;
	vector<VectorXd> kernels;
	kernels.resize(sample_of_center*sample_of_center, VectorXd(9));
	cout << "[match guassian] matching intensity to exr..." << endl;
	//Mat img_result = Mat(Size(width, height), CV_8UC3, Scalar(0, 0, 0));
	vector<float> result(width*height * 3, 0);
	vector<double> sigma_bgr({ 0.66, 0.55, 0.64 });
	//for (int i = 0; i < relationship.size(); i++)
	int select_rgb = 1;
	{
		/*double sigma_init, sigma, loss_old, loss = 1000000;
		bool flag = false;
		for (sigma_init=0.5; sigma_init < 1.5; sigma_init += 0.1)
		{
			loss_old = loss;
			loss = compute_sigma(data[i], kernels, sigma_init, sigma_init);
			if (!flag && loss_old - loss>0 )
			{
				flag = true;
			}
			else if(flag && loss_old - loss < 0)
			{
				break;
			}
		}
		sigma_init -= 0.1;
		cout << "init sigma: " << sigma_init << " ";
		flag = false;
		for (sigma = sigma_init - 0.1; sigma < sigma_init+0.1; sigma += 0.01)
		{
			loss_old = loss;
			loss = compute_sigma(data[i], kernels, sigma, sigma);
			if (!flag && loss_old - loss > 0)
			{
				flag = true;
			}
			else if (flag && loss_old - loss < 0)
			{
				break;
			}
		}
		sigma -= 0.01;
		sigma_bgr[select_rgb] = sigma;
		cout <<"final use sigma: "<< sigma_bgr[select_rgb] << endl; */

		get_kernels(kernels, sigma_bgr[select_rgb], sigma_bgr[select_rgb]);
		char result_file[MAX_PATH], input_txt[MAX_PATH], output_csv[MAX_PATH];
		sprintf(result_file, "%s_%s.exr",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		match_with_location_to_exr(result, result_file,
			relationship[select_rgb], img, kernels, width, height);//, (RGB)select_rgb);
	}
	p.endAndPrint();

	cout << "total kernels count:" << kernels.size() << endl;
	//char result_file[MAX_PATH];
	//sprintf(result_file, "%s/valid_result.png", output_prefix);
	//// rotate
	//img_result = img_result(cv::Rect(0, 0, width, height));
	//Mat rgb(Size(height, width), CV_8UC3, Scalar(0, 0, 0));
	//for (int y = 0; y < height; y++)
	//{
	//	for (int x = 0; x < width; x++)
	//	{
	//		rgb.at<Vec3b>(width - 1 - x, y) = img_result.at<Vec3b>(y, x);
	//	}
	//}
	//imwrite(result_file, rgb);
}