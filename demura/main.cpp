#include "location.h"
#include "gaussian_filter.h"
#include "preprocess.h"
#include "pentile.h"
#include "gmatch.h"
#include "diff.h"
#include <direct.h>
#include <fstream>
using namespace std;
using namespace cv;
using namespace Eigen;

void test_g32()
{
	int real_value = 80; //255.0 / 32;
	cout << "real_value: " << real_value << endl;
	Mat img = imread("G32.bmp") * real_value;
	GaussianBlur(img, img, Size(5, 5), 2, 2);

	cout << "rows*cols: " << img.rows << "*" << img.cols << endl;
	resize(img, img, Size(), 1 / 8.0, 1 / 8.0);
	cout << "rows*cols: " << img.rows << "*" << img.cols << endl;
	//cout << img.type() << endl;
	//cout << CV_8UC3 << endl;
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			if (img.at<Vec3b>(i, j)[0] != 255
				|| img.at<Vec3b>(i, j)[1] != 255
				|| img.at<Vec3b>(i, j)[2] != 255)
			{
				img.at<Vec3b>(i, j)[0] = 0;
				img.at<Vec3b>(i, j)[1] = 0;
				img.at<Vec3b>(i, j)[2] = 0;
			}
		}
	}
	imshow("G32", img);
	waitKey(0);
}

void test_pentile()
{
	Mat rgb = imread("./output/test_rgb.bmp"), pentile;
	rgb2pentile(rgb, pentile);
	imwrite("./output/pentile2.bmp", pentile);
	pentile2rgb(pentile, rgb);
	imwrite("./output/bask_test_rgb2.bmp", rgb);
}

void merge()
{
	Mat r = imread("./output/R16_result.png"), 
		g = imread("./output/B16_result.png"),
		b = imread("./output/B16_result.png"),
		bgr(r.size(), r.type());
	for (int y = 0; y < r.rows; y++)
	{
		for (int x = 0; x < r.cols; x++)
		{
			bgr.at<Vec3b>(y, x)[0] = b.at<Vec3b>(y, x)[0];
			bgr.at<Vec3b>(y, x)[1] = g.at<Vec3b>(y, x)[1];
			bgr.at<Vec3b>(y, x)[2] = r.at<Vec3b>(y, x)[2];
		}
	}
	imwrite("./output/merge.png", bgr);
}
#include <algorithm>
extern void draw_cross(Mat& img, int x, int y, Vec3b color);
void generate_compensate_value(const char* img_file)
{
	/*Mat img = imread(img_file), bgr[3],
		matr(img.size(), CV_8UC3),
		matg(img.size(), CV_8UC3),
		matb(img.size(), CV_8UC3),
		pentile_r, pentile_g, pentile_b;
	cout << img.type() <<","<< CV_8UC3 << endl;
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			matb.at<Vec3b>(y, x)[0] = img.at<Vec3b>(y, x)[0];
			matg.at<Vec3b>(y, x)[1] = img.at<Vec3b>(y, x)[1];
			matr.at<Vec3b>(y, x)[2] = img.at<Vec3b>(y, x)[2];
		}
	}

	rgb2pentile(matb, pentile_b);
	imwrite("./output/result_pentile_b.bmp", pentile_b);
	rgb2pentile(matg, pentile_g);
	imwrite("./output/result_pentile_g.bmp", pentile_g);
	rgb2pentile(matr, pentile_r);
	imwrite("./output/result_pentile_r.bmp", pentile_r);*/
	int h = 2436, w = 752 / 2 * 3;
	Mat img = /*imread(img_file)*/Mat(Size(w, h), CV_8UC3), pentile;
	vector<int> value{ 4,8,12,16,20,24,28,32,64 }, //{ 16, 32, 48, 64, 96, 128, 160, 192, 224, 250 },
		//({ 4,8,12,16,20,24,28,32 }), //({ 12,16,19,22,25,29,32 }),
		cross_value{ 32,32,32,32,32,32,16,16,16 };
	//{ 16 * 2, 32 * 2, 48 * 2, 64 * 2, 96 * 2, 255, 160 / 2, 192 / 2, 224 / 2, 250 / 2 };
	//{ 32,32,32,32,32,32,32,16 };
	vector<int> need_x({ 60,1060 }), need_y({ 130, h - 130 });
	for (int v = 0; v < value.size(); v++)
	{
		for (int y = 0; y < img.rows; y++)
		{
			for (int x = 0; x < img.cols; x++)
			{
				img.at<Vec3b>(y, x)[0] = value[v];
			}
		}

		for (auto y : need_y)
		{
			for (auto x : need_x)
			{
				draw_cross(img, x, y, Vec3b(cross_value[v], 0, 0));
			}
		}
		char path[MAX_PATH];
		sprintf_s(path, "./output/b_%d.bmp", value[v]);
		imwrite(path, img);
		rgb2pentile(img, pentile);
		sprintf_s(path, "./output/pentile_b_%d.bmp", value[v]);
		imwrite(path, pentile);
	}
	/*for (auto y : need_y)
	{
		for (auto x : need_x)
		{
			draw_cross(img, x, y, Vec3b(0, 64, 0));
		}
	}
	char path[MAX_PATH];
	sprintf_s(path, "./output/result_of_g_%d_.bmp", 32);
	imwrite(path, img);
	rgb2pentile(img, pentile);
	sprintf_s(path, "./output/result_of_g_%d_pentile.bmp", 32);
	imwrite(path, pentile);*/
}

void convert(const char* img_file)
{
	vector<Vec3b> color;
	{
		for(int i=0;i<23;i++)
			color.push_back(Vec3b(0, 0, 0));
		/*for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(i * 255 / 4, 0, 0));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(0, i * 255 / 4, 0));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(0, 0, i * 255 / 4));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(0, i * 255 / 4, i * 255 / 4));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(i * 255 / 4, 0, i * 255 / 4));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(i * 255 / 4, i * 255 / 4, 0));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(i * 255 / 4, i * 255 / 4, i * 255 / 4));
		}
		for (int i = 1; i <= 4; i++)
		{
			color.push_back(Vec3b(i * 255 / 8, i * 255 / 4, i * 255 / 8));
		}*/

		color.push_back(Vec3b(255, 0, 0)); //23
		color.push_back(Vec3b(255, 0, 0)); //24
		color.push_back(Vec3b(0, 255, 0)); //25
		color.push_back(Vec3b(0, 255, 0)); //26
		color.push_back(Vec3b(0, 0, 255)); //27
		color.push_back(Vec3b(0, 0, 255)); //28
		color.push_back(Vec3b(0, 255, 255));  //29
		color.push_back(Vec3b(255, 0, 255));  //30
		color.push_back(Vec3b(255, 255, 0));  //31
		color.push_back(Vec3b(255, 255, 255));//32
	}
	Mat img = imread(img_file), out(Size(img.cols, img.rows), CV_8UC3);
	cout << img.cols << " " << img.rows << endl;
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			out.at<Vec3b>(y, x) = color[img.at<Vec3b>(y, x)[1]];
		}
	}
	imwrite("./output/2.png", out);
}

void convert2(const char* img_file, const char* result_file)
{
	Vec3b b(255, 0, 0), g(0, 255, 0), r(0, 0, 255), ye(0, 255, 255),
		black(0,0,0);
	Mat img = imread(img_file),
		result = imread(result_file),
		out(Size(img.cols, img.rows), CV_8UC3, black);
	cout << img.cols << " " << img.rows << endl;
	for (int y = 0; y < img.rows; y++)
	{
		for (int x = 0; x < img.cols; x++)
		{
			// <24
			if (result.at<Vec3b>(y, x) == b)
				out.at<Vec3b>(y, x) = img.at<Vec3b>(y, x)[1] < 24 ? black : b;
			// [24,28)
			else if (result.at<Vec3b>(y, x) == g)
				out.at<Vec3b>(y, x) = 
				(img.at<Vec3b>(y, x)[1]>=24 && img.at<Vec3b>(y, x)[1] < 28)? black : 
				(img.at<Vec3b>(y, x)[1] < 24 ? r : g);
			// [24,28)
			else if (result.at<Vec3b>(y, x) == r)
				out.at<Vec3b>(y, x) = 
				(img.at<Vec3b>(y, x)[1] < 32) ? black : ye;
			/*else if (result.at<Vec3b>(y, x) == ye)
				out.at<Vec3b>(y, x) =
				(img.at<Vec3b>(y, x)[1] >= 30) ? black : Vec3b(255, 255, 255);*/
		}
	}
	imwrite("./output/is_correct_region_in_limit1.png", out);
}

void validation(const char* input_origin, const char* input_compute, 
	const char* output_file,
	ofstream& fout, const int offset = 0, const int scale = 1)
{
	Vec3b b(255, 0, 0), g(0, 255, 0), r(0, 0, 255), ye(0, 255, 255), black(0, 0, 0);
	Mat origin = imread(input_origin),
		result = imread(input_compute);
	if (origin.data == nullptr)
	{
		cout << "read error : " << input_origin << endl;
		return;
	}
	if (result.data == nullptr)
	{
		cout << "read error : " << input_compute << endl;
		return;
	}
	Mat	out(Size(origin.cols, origin.rows), CV_8UC3, black);
	//cout <<"size : "<< origin.cols << " " << origin.rows << endl;
	//int offset = 1;
	int mean = 0, tmp, cnt, mmax = 0, mmin = 0;
	for (int y = 0; y < origin.rows; y++)
	{
		for (int x = 0; x < origin.cols; x++)
		{
			//// <24
			//if (result.at<Vec3b>(y, x) == b)
			//	out.at<Vec3b>(y, x) = img.at<Vec3b>(y, x)[1] < 24 ? black : b;
			//// [24,28)
			//else if (result.at<Vec3b>(y, x) == g)
			//	out.at<Vec3b>(y, x) =
			//	(img.at<Vec3b>(y, x)[1] >= 24 && img.at<Vec3b>(y, x)[1] < 28) ? black :
			//	(img.at<Vec3b>(y, x)[1] < 24 ? r : g);
			//// [24,28)
			//else if (result.at<Vec3b>(y, x) == r)
			//	out.at<Vec3b>(y, x) =
			//	(img.at<Vec3b>(y, x)[1] < 32) ? black : ye;
			/*if(result.at<Vec3b>(y, x)[1])
				if (result.at<Vec3b>(y, x)[1] > origin.at<Vec3b>(y, x)[1] + offset)
					out.at<Vec3b>(y, x) = r;
				else if (result.at<Vec3b>(y, x)[1] < origin.at<Vec3b>(y, x)[1] - offset)
					out.at<Vec3b>(y, x) = b;
				else
					out.at<Vec3b>(y, x) = g;*/
			/*if (origin.at<Vec3b>(y, x)[1])
				if (origin.at<Vec3b>(y, x)[1] > 16)
					out.at<Vec3b>(y, x)[1] = origin.at<Vec3b>(y, x)[1]*8;*/
			if (result.at<Vec3b>(y, x)[1] > origin.at<Vec3b>(y, x)[1]
				&& result.at<Vec3b>(y, x)[1] - origin.at<Vec3b>(y, x)[1] > offset)
			{
				tmp = result.at<Vec3b>(y, x)[1] - origin.at<Vec3b>(y, x)[1];
				out.at<Vec3b>(y, x)[1] = tmp;
				mean += tmp;
				cnt++;
				mmax = max(tmp, mmax);
			}
			else if (result.at<Vec3b>(y, x)[1] < origin.at<Vec3b>(y, x)[1]
				&& origin.at<Vec3b>(y, x)[1] - result.at<Vec3b>(y, x)[1] > offset)
			{
				tmp = origin.at<Vec3b>(y, x)[1] - result.at<Vec3b>(y, x)[1];
				out.at<Vec3b>(y, x)[2] = tmp;
				mean += tmp;
				cnt++;
				mmin = min(-tmp, mmin);
			}
			/*else
				out.at<Vec3b>(y, x)[2] = 255;*/
		}
	}
	imwrite(output_file, out);
	fout << "min: " << mmin << endl
		<< "max: " << mmax << endl
		<< "mean: " << mean / (1.0*cnt) << endl << endl;
}

void validtation_batch(const char* input_origin, const char* input_compute,
	const char* output_file, const int offset = 0, const int scale = 0)
{
	_mkdir(output_file);
	//cout << argv[3] << endl;
	std::vector<cv::String> filenames1, filenames2;
	cv::String folder1 = input_origin, folder2 = input_compute;
	folder1 += "/*.bmp";
	folder2 += "/*.bmp";
	//cout << folder1 << endl;
	//cout << folder2 << endl;
	cv::glob(folder1, filenames1);
	cv::glob(folder2, filenames2);
	//cout << filenames1.size() << endl;
	string outtxt = output_file;
	outtxt += "/out.txt";
	ofstream out(outtxt);
	for (size_t i = 0; i < filenames1.size(); ++i)
	{
		out << filenames1[i] << endl << filenames2[i] << endl << endl;
		string tmp = output_file + filenames1[i].substr(filenames1[i].find_last_of("\\"));
		//cout << tmp << endl;
		validation(filenames1[i].c_str(), filenames2[i].c_str(), tmp.c_str(), out, offset, scale);
	}
	out.close();
	return ;
}

// 比较g4-8-12-16-20-24-28-32的demura值是否符合递增趋势
void tmp()
{
	//for (int pid = 1; pid <= 9; pid++)
	int pid = 2;
	{
		char path[MAX_PATH];
		int h = 2436, w = 752 / 2 * 3;
		Mat out = Mat(Size(h, w), CV_8UC3);
		int size = 8;
		vector<double> val(size);
		pair<int, int> tmp;
		sprintf_s(path, "E:/document/研二/demura/20200407_9块屏/output%i/gaussian_val.txt", pid);

		ifstream in(path);
		if (!in)
		{
			cout << "open fail: " << path << endl;
			return;
		}
		cout <<"process "<< pid << endl;
		int i;
		while (!in.eof())
		{
			in >> tmp.first >> tmp.second;
			//cout << d1 <<","<< d2 << endl;
			for (i = 0; i < size; i++)
			{
				in >> val[i];
			}
			if (tmp.first >= 0 && tmp.first < w && tmp.second >= 0 && tmp.second < h)
			{
				/*for (i = 1; i < size; i++)
				{
					if (val[i] < val[i - 1])
					{
						out.at<Vec3b>(tmp.first, tmp.second) = Vec3b(0, 0, 255);
						break;
					}
				}
				if (i == size)
					out.at<Vec3b>(tmp.first, tmp.second) = Vec3b(0, 255, 0);*/
				if (val[2] < val[1])
				{
					out.at<Vec3b>(tmp.first, tmp.second) = Vec3b(0, 0, 255);
				}
				else
					out.at<Vec3b>(tmp.first, tmp.second) = Vec3b(0, 255, 0);
			}
		}
		sprintf_s(path, "E:/document/研二/demura/20200407_9块屏/output%i/gaussian_val.png", pid);
		imwrite(path, out);
	}
}

// demura后的值-1
void decrease()
{
	char path[MAX_PATH], file[MAX_PATH];
	sprintf_s(path, "E:/document/研二/demura/20200407_9块屏/output1");
	sprintf_s(file, "%s/result_rotate_90.bmp", path);
	Mat img = imread(file), pentile;
	for (int y = 0; y < 400; y++)
	{
		
		for (int x = 0; x < 110; x++)
		{
			img.at<Vec3b>(y, x)[1]--;
		}
		for (int x = 190; x < 300; x++)
		{
			img.at<Vec3b>(y, x)[1] -= 2;
		}
		
		for (int x = 300; x < 500; x++)
		{
			img.at<Vec3b>(y, x)[1]--;
		}
		for (int x = 600; x < img.cols - 300; x++)
		{
			img.at<Vec3b>(y, x)[1] -= 2;
		}

		for (int x = img.cols - 300; x < img.cols - 190; x++)
		{
			img.at<Vec3b>(y, x)[1]--;
		}
		for (int x = img.cols - 110; x < img.cols; x++)
		{
			img.at<Vec3b>(y, x)[1] -= 2;
		}
	}

	int h = 2436, w = 752 / 2 * 3;
	vector<int> need_x({ 60,1060 }), need_y({ 130, h - 130 });
	for (auto y : need_y)
	{
		for (auto x : need_x)
		{
			draw_cross(img, x, y, Vec3b(0, 32, 0));
		}
	}
	sprintf_s(file, "%s/result_rotate_90_decrease.bmp", path);
	imwrite(file, img);
	rgb2pentile(img, pentile);
	sprintf_s(file, "%s/result_pentile_decrease.bmp", path);
	imwrite(file, pentile);
}

void add_special(char* input, char *output, char *output_bmp, int id, int val)
{
	//char path[MAX_PATH], file[MAX_PATH];
	//sprintf_s(path, "E:/document/研二/demura/20200407_9块屏/output1_sqrt");
	//sprintf_s(file, "%s/result_rotate_90.bmp", path);
	//Mat img = imread(file), pentile;
	Mat img = imread(input), pentile;
	int h = 2436, w = 752 / 2 * 3;
	//vector<int> need_x({ 60,1060 }), need_y({ 130, h - 130 });
	{
		// 1号屏
		int y = img.rows / 2, x = id * 100;
		/*for(int dy = -1; dy <=1; dy++)
			for(int dx = -1; dx <= 1; dx++)
				draw_cross(img, x+dx, y + dy, Vec3b(0, 32, 0));*/
		draw_cross(img, x, y, Vec3b(0, val, 0));
	}
	//sprintf_s(file, "%s/sqrt_with_mark.bmp", path);
	//imwrite(file, img);
	imwrite(output, img);
	rgb2pentile(img, pentile);
	//sprintf_s(file, "%s/sqrt_with_mark_pentile.bmp", path);
	//imwrite(file, pentile);
	imwrite(output_bmp, pentile);
}

// 加上scale: v->(v-16)*scale+16
void add_scale()
{
	float scale;
	Mat origin = imread("E:/document/研二/demura/20200407_9块屏/output1_sigma0.62/result_rotate_90.bmp"), 
		out(Size(origin.cols, origin.rows), CV_8UC3, Scalar(0, 0, 0)), pentile;
	char out_path[MAX_PATH];
	for (int i = 16; i <= 24; i++)
	{
		scale = i * 0.05;
		sprintf_s(out_path, "%s_scale_%.2f.bmp", 
			"E:/document/研二/demura/20200407_9块屏/output1_sigma0.62/scaled/result_rotate_90", scale);
		for (int y = 0; y < origin.rows; y++)
		{
			for (int x = 0; x < origin.cols; x++)
			{
				out.at<Vec3b>(y, x)[1] = (origin.at<Vec3b>(y, x)[1] - 16)*scale + 16;
			}
		}
		imwrite(out_path, out);
		sprintf_s(out_path, "%s_scale_%.2f.bmp",
			"E:/document/研二/demura/20200407_9块屏/output1_sigma0.62/scaled/pentile", scale);
		rgb2pentile(out, pentile);
		imwrite(out_path, pentile);

	}
}

// 读入一张图片，计算每个点的location，然后把亮度值输出(exr)
void intensity_to_exr(const char* inpath, const char* outpath,
	const char* infilename, const char* outfilename)
{
	_mkdir(outpath);
	const int pentile_width = 752, pentile_height = 2436;
	char b_file[MAX_PATH], g_file[MAX_PATH], r_file[MAX_PATH],
		//output
		mask_file[MAX_PATH], out_file[MAX_PATH];
	sprintf_s(b_file, "%s/%s", inpath, infilename);
	sprintf_s(g_file, "%s/%s", inpath, infilename);
	sprintf_s(r_file, "%s/%s", inpath, infilename);
	sprintf_s(mask_file, "%s/mask.png", inpath);
	sprintf_s(out_file, "%s/%s", outpath, outfilename);
	Mat mask = imread(mask_file, CV_8UC1);
	const RGB b = BLUE, g = GREEN, r = RED;
	vector<Mat> rgb;
	Mat bb = imread(b_file), gg = imread(g_file), rr = imread(r_file);
	rgb.push_back(bb);
	rgb.push_back(gg);
	rgb.push_back(rr);
	if (rgb[b].data == nullptr || rgb[g].data == nullptr || rgb[r].data == nullptr)
	{
		cout << "file read error" << endl << b_file << ", or " << g_file << ", or " << r_file << endl;
	}
	if (mask.data == nullptr)
		preprocess(rgb, 0, mask_file, mask);

	vector<vector<Point>> centers_error(3);
	vector<vector<vector<LED_info>>> relationship(3);
	{
		if (tmp_valid_find_location(rgb, mask, outpath, pentile_height, relationship))
			return ;
		cout << "--------------------" << endl;
		//compute_intensity_to_exr(relationship, gg, inpath, out_file,
			//pentile_height, pentile_width / 2 * 3);
		vector<int> val(3, 0), cnt(3, 0);
		for (auto r : relationship[1])
		{
			for (auto re : r)
			{
				if (re.state == VALID && re.locate.y>=0 
					&& re.locate.y< pentile_width/2*3)
				{
					if (re.locate.x >= 0 && re.locate.x < pentile_height / 3)
					{
						val[0] += gg.at<Vec3b>(re.pixel.y, re.pixel.x)[1];
						cnt[0]++;
					}
					else if (re.locate.x < pentile_height / 3 * 2)
					{
						val[1] += gg.at<Vec3b>(re.pixel.y, re.pixel.x)[1];
						cnt[1]++;
					}
					else if(re.locate.x < pentile_height)
					{
						val[2] += gg.at<Vec3b>(re.pixel.y, re.pixel.x)[1];
						cnt[2]++;
					}
				}
			}
		}
		for (int i = 0; i < val.size(); i++)
		{
			cout << val[i] << " " << cnt[i] << " ";
			cout << val[i] * 1.0 / cnt[i] << " " << endl;
		}
	}
	cout << "the end." << endl;
}

// 读入一系列图片(如g4-g32)，根据这些图片，计算出亮度-mura值对应关系，然后反向输出平均值对应的mura值
#define GEXPO_RANGE // 16_224 4_32
void compute_demura_value_use_range(char* inpath, char* outpath, 
	char* prefix, char* suffix, const vector<int>& ms)
{
	//String inpath("E:/coding/gray_game/demura/input2.2.2_pentile"),
	//	outpath("E:/coding/gray_game/demura/output");
	//String inpath("E:\\document\\研二\\demura\\20200407_9块屏\\output1"), // ("./input2.2_20200404"),
	//	outpath("E:\\document\\研二\\demura\\20200409-9块屏结果\\output1"), 
	//	prefix("2号屏pentile_g_"),
	//	postfix("-曝光时间750ms-增益3-施耐德镜头光圈5.7");
	//String valid_path("E:\\document\\研二\\demura\\20200409-9块屏结果\\output1"); 
	// ("./input2.2_20200404");
	
	cout << "inpath [argv1]" << inpath << endl
		<< "outpath [argv2]" << outpath << endl
		<< "prefix [argv3]" << prefix << endl
		<< "suffix [argv4]" << suffix << endl;
	_mkdir(outpath);
	int pentile_width = 752, pentile_height = 2436;
	int target_g_id;
	// 处理输入文件名字
	char b_file[MAX_PATH], g_file[MAX_PATH], r_file[MAX_PATH], range_file[7][MAX_PATH],
		//output
		mask_file[MAX_PATH], cross_file[MAX_PATH],
		output_csv[MAX_PATH];
	
	/// 1
	vector<int> capture_pentile_g_value({ 4,8,12,16,20,24,28,32 });
	//const vector<int> ms{ 2000,2000,1800,750,1200,800,600,500 };
	target_g_id = 3;
	sprintf_s(b_file, "%s/%s%02d%s%04d.bmp", inpath, prefix,
		capture_pentile_g_value[target_g_id], suffix, ms[target_g_id]);
	sprintf_s(g_file, "%s/%s%02d%s%04d.bmp", inpath, prefix,
		capture_pentile_g_value[target_g_id], suffix, ms[target_g_id]);
	sprintf_s(r_file, "%s/%s%02d%s%04d.bmp", inpath, prefix,
		capture_pentile_g_value[target_g_id], suffix, ms[target_g_id]);
	cout << g_file << endl;
	for (int i = 0; i < capture_pentile_g_value.size(); i++)
	{
		sprintf_s(range_file[i], "%s/%s%02d%s%04d.bmp", inpath, prefix,
			capture_pentile_g_value[i], suffix, ms[i]);
	}

	// mask
	sprintf_s(mask_file, "%s/mask.png", inpath);
	Mat mask = imread(mask_file, CV_8UC1);
	const RGB b = BLUE, g = GREEN, r = RED;
	vector<Mat> rgb;
	Mat bb = imread(b_file), gg = imread(g_file), rr = imread(r_file);
	rgb.push_back(bb);
	rgb.push_back(gg);
	rgb.push_back(rr);
	if (rgb[b].data == nullptr || rgb[g].data == nullptr || rgb[r].data == nullptr)
	{
		cout << "file read error" << endl << b_file << endl
			<< "or " << g_file << endl
			<< "or " << r_file << endl;
	}
	if (mask.data == nullptr)
		preprocess(rgb, 0, mask_file, mask);
	sprintf_s(output_csv, "%s/pentile_rgb_relationship.csv", outpath);

	cout << "--------------------" << endl;
	vector<vector<Point>> centers_error(3);
	vector<vector<vector<LED_info>>> relationship(3);
	{
		RGB select_rgb = BLUE;
		if (find_OLED_location_with_rgb_combination(rgb, mask,
			outpath, pentile_height, relationship))
		{
			cout << "something error!" << endl;
			return ;
		}
		cout << "--------------------" << endl;
		// to test 4,8,11,16,23,32,64 value
		vector<Mat> pic;
		{
			for (int i = 0; i < capture_pentile_g_value.size(); i++)
			{
				Mat tmp = imread(range_file[i]);
				if (tmp.data == NULL)
				{
					cout << range_file[i] << "read error." << endl;
					return ;
				}
				pic.push_back(tmp);
				cvtColor(*pic.rbegin(), *pic.rbegin(), COLOR_BGR2GRAY);
			}
		}
		compute_dumura(relationship, capture_pentile_g_value, ms,
			pic, target_g_id, outpath, pentile_height, pentile_width / 2 * 3);
	}
	cout << "the end." << endl;
	return ;
}

// demura后的拍摄图片，重新计算对应的亮度-mura关系下的mura值
void compute_single_validation(char* inpath, char* outpath, char* prefix, char* valid_path)
{
	cout << "inpath [argv1]" << inpath << endl
		<< "outpath [argv2]" << outpath << endl
		<< "prefix [argv3]" << prefix << endl
		<< "suffix [argv4]" << valid_path << endl;
	_mkdir(outpath);
	int pentile_width = 752, pentile_height = 2436;
	int target_g_id;
	
	vector<int> capture_pentile_g_value({ 4,8,12,16,20,24,28,32 });
	char valid_b_file[MAX_PATH], valid_g_file[MAX_PATH], valid_r_file[MAX_PATH],
		//output
		valid_mask_file[MAX_PATH];
	//int valid_ms = 1000;
	sprintf_s(valid_b_file, "%s/%s.bmp", valid_path, prefix);
	sprintf_s(valid_g_file, "%s/%s.bmp", valid_path, prefix);
	sprintf_s(valid_r_file, "%s/%s.bmp", valid_path, prefix);
	sprintf_s(valid_mask_file, "%s/mask.png", valid_path);
	Mat valid_mask = imread(valid_mask_file, CV_8UC1);
	const RGB b = BLUE, g = GREEN, r = RED;
	vector<Mat> rgb;
	Mat bb = imread(valid_b_file), gg = imread(valid_g_file), rr = imread(valid_r_file);
	rgb.push_back(bb);
	rgb.push_back(gg);
	rgb.push_back(rr);
	if (rgb[b].data == nullptr || rgb[g].data == nullptr || rgb[r].data == nullptr)
	{
		cout << "file read error" << endl << valid_b_file << ", or " << valid_g_file << ", or " << valid_r_file << endl;
	}
	if (valid_mask.data == nullptr)
		preprocess(rgb, 0, valid_mask_file, valid_mask);

	vector<vector<Point>> centers_error(3);
	vector<vector<vector<LED_info>>> relationship(3);
	{
		if (tmp_valid_find_location(rgb, valid_mask, outpath, pentile_height, relationship))
			return;
		cout << "--------------------" << endl;
		compute_dumura_single_pic(relationship, capture_pentile_g_value, gg,
			inpath, outpath, pentile_height, pentile_width / 2 * 3);

	}
	cout << "the end." << endl;
	return ;
}

//#define VALIDATION
int main(int argc, char* argv[])
{
	//draw_pattern2("./output/test2", pentile_width, pentile_height);
	//test_pentile();
	//merge();
	//generate_compensate_value("./output/result_rotate_90.bmp");
	//convert2("./output/result_rotate_90.bmp", "./output/valid_result.png");
	/*validation("E:/document/研二/demura/20200416-20200415/origin/valid_result.png",
		//"E:/document/研二/demura/20200409-9块屏结果/output9/valid_result.png", 
		//"E:/document/研二/demura/20200409-9块屏结果/output9/valid_result_campare2.png");
		"E:/document/研二/demura/20200416-20200415/decrease/valid_result.png",
		//"E:/document/研二/demura/20200407_9块屏/output9/result_rotate_90_large16.bmp");
		"E:/document/研二/demura/20200416-20200415/decrease/difference_with_origin.png");*/
	//tmp();
	//decrease();
	//add_special(argv[1], argv[2], argv[3], stoi(argv[4]), stoi(argv[5]));
	/*validation("E:/document/研二/demura/重复试验/关屏并重启屏/第一轮/1号屏750_origin_pentile-曝光时间750ms-增益3-施耐德新镜头光圈5.4.bmp",
		"E:/document/研二/demura/重复试验/不关屏/第二轮/1号屏750_origin_pentile-曝光时间750ms-增益3-施耐德新镜头光圈5.4.bmp",
		"E:/document/研二/demura/重复试验/validation/1号屏750_origin_pentile-曝光时间750ms-增益3-施耐德新镜头光圈5.4.bmp");*/
	//add_scale();

	//_mkdir(argv[2]);
	//std::vector<cv::String> filenames;
	//cv::String folder = argv[1];
	//folder += "/*.bmp";
	//cv::glob(folder, filenames);
	////cout << filenames1.size() << endl;
	//string outtxt = argv[2];
	//outtxt += "/out.txt";
	//ofstream out(outtxt);
	//for (size_t i = 0; i < filenames.size(); ++i)
	//{
	//	out << filenames[i] << endl << argv[5] << endl << endl;
	//	string tmp = argv[2] + filenames[i].substr(filenames[i].find_last_of("\\"));
	//	//cout << tmp << endl;
	//	validation(filenames[i].c_str(), argv[5], tmp.c_str(), out,
	//		atoi(argv[3]), atoi(argv[4]));
	//}
	//out.close();
	//return 0;

	/*vector<float> pic1, pic2;
	int w, h;
	cout << argv[1] << endl << argv[2] << endl;
	read_exr(argv[1], w, h, pic1);
	read_exr(argv[2], w, h, pic2);
	Mat	out(Size(w, h), CV_8UC3);
	int x1, y1, x2, y2;
	//vector<float> out(w*h * 3, 0);
	float a, mmin=0, mmax=0;
	for (int y = 100; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			a = pic1[(y*w + x) * 3 + 1] - pic2[(y*w + x) * 3 + 1];
			if (a > 0)
			{
				out.at<Vec3b>(y, x)[1] = __min(255, a * 1000);
				if (a > mmax)
				{
					x1 = x;
					y1 = y;
					mmax = a; //__max(a / pic1[(y*w + x) * 3 + 1], mmax);
				}
			}
			else if(a<0)
			{
				out.at<Vec3b>(y, x)[2] = __min(255, -a * 1000);
				if (a < mmin)
				{
					x2 = x;
					y2 = y;
					mmin = a;
				}
			}
		}
	}
	//save_exr_with_float(argv[3], w, h, out);
	imwrite(argv[3], out);
	cout << mmin << " " << mmax << endl
		<< x1 << "," << y1 << endl
		<< x2 << "," << y2 << endl;
	cout << "end." << endl;
	return 0;*/

	if (argv[1][0] == 'a')
		intensity_to_exr(argv[2], argv[3], argv[4], argv[5]);
	else if (argv[1][0] == 'b')
	{
		vector<int> ms{ 2000,2000,1800,750,1200,800,600,500 };
		if (argc == 5 + 8)
		{
			for (int i = 1; i <= 8; i++)
			{
				ms[i] = atoi(argv[i + 5]);
			}
		}
		compute_demura_value_use_range(argv[2], argv[3], argv[4], argv[5], ms);
	}
	else if (argv[1][0] == 'c')
		compute_single_validation(argv[2], argv[3], argv[4], argv[5]);
	else if (argv[1][0] == 'd')
		validtation_batch(argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]));
	return 0;
}

/*
亮度(demura后高斯的最小二乘拟合值)直接输出为exr
:: mode: a
:: para: input_path 
::		 output_path 
::		 input_file 
::		 output_exr_file_name_prefix
demura.exe a E:/document/研二/demura/重复试验/三号屏/第一轮 ^
E:/document/研二/demura/重复试验/三号屏/第一轮/result ^
3号屏_pentile_g_16.bmp曝光时间2000ms-增益5-施耐德新镜头光圈5.4.bm.bmp ^
result

通过八张图片的亮度得到亮度-mura值关系，线性插值得到demura值
:: mode: b
:: para: input_path
::		 output_path
::		 prefix
::		 suffix
demura.exe b E:/document/研二/demura/20200407_9块屏/1号屏 ^
E:/document/研二/demura/20200407_9块屏/output1_no_interpolation ^
1号屏pentile_g_ ^
-曝光时间750ms-增益3-施耐德镜头光圈5.7

验证demura后的图片在原【亮度-mura值关系】中的mura值
:: mode: c
:: para: input_path 
::		 output_path 
::		 prefix 
::		 valid_path
demura.exe c E:/document/研二/demura/20200407_9块屏/output1 ^
E:/document/研二/demura/20200409-9块屏结果/output1 ^
1号屏result_pentile-曝光时间750ms-增益3-施耐德镜头光圈5.7 ^
E:/document/研二/demura/20200409-9块屏结果/output1

两张原图直接相减
:: mode: d
:: para: input_path1 input_path2 output_path difference_offset scale
demura.exe d E:/document/研二/demura/重复试验/一号屏/关屏并重启屏/第一轮 E:/document/研二/demura/重复试验/一号屏/关屏并重启屏/第二轮 ^
E:/document/研二/demura/重复试验/一号屏/关屏并重启屏/第一轮/result_test 0 1
*/