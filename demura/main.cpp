#include "location.h"
#include "gaussian_filter.h"
#include "preprocess.h"
#include "pentile.h"
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
	Mat img = imread(img_file)/*(Size(w, h), CV_8UC3)*/, pentile;
	vector<int> value({ 4,8,12,16,20,24,28,32 }), //({ 12,16,19,22,25,29,32 }),
		cross_value({ 32,32,32,32,32,32,32,16 });
	vector<int> need_x({ 60,1060 }), need_y({ 130, h - 130 });
	for (int v = 0; v < value.size(); v++)
	{
		for (int y = 0; y < img.rows; y++)
		{
			for (int x = 0; x < img.cols; x++)
			{
				img.at<Vec3b>(y, x)[1] = value[v];
			}
		}

		for (auto y : need_y)
		{
			for (auto x : need_x)
			{
				draw_cross(img, x, y, Vec3b(0, cross_value[v], 0));
			}
		}
		char path[MAX_PATH];
		sprintf_s(path, "./output/g_%d.bmp", value[v]);
		imwrite(path, img);
		rgb2pentile(img, pentile);
		sprintf_s(path, "./output/pentile_g_%d.bmp", value[v]);
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

//#define VALIDATION
int main(int argc, char* argv[])
{
	/*String inpath("E:/coding/gray_game/demura/input2.2.2_pentile"),
		outpath("E:/coding/gray_game/demura/output");*/
	//String inpath("E:\\document\\研二\\demura\\20200407_9块屏\\2号屏"), // ("./input2.2_20200404"),
	//	outpath("E:\\document\\研二\\demura\\20200407_9块屏\\output2"), 
	//	prefix("2号屏pentile_g_"),
	//	postfix("-曝光时间750ms-增益3-施耐德镜头光圈5.7");
	String inpath(argv[1]),
		outpath(argv[2]),
		prefix(argv[3]),
		postfix(argv[4]);
	_mkdir(outpath.c_str());
	int pentile_width = 752, pentile_height = 2436;
	//draw_pattern2("./output/test2", pentile_width, pentile_height);
	//test_pentile();
	//merge();
	/*generate_compensate_value("./output/result_rotate_90.bmp");
	system("pause");
	return 0;*/
	/*convert2("./output/result_rotate_90.bmp", "./output/valid_result.png");
	system("pause");
	return 0;*/

#ifndef VALIDATION
	char b_file[MAX_PATH], g_file[MAX_PATH], r_file[MAX_PATH], range_file[7][MAX_PATH],
		//output
		mask_file[MAX_PATH], cross_file[MAX_PATH],
		output_csv[MAX_PATH];
	int ms = 720;
	sprintf_s(b_file, "%s/%s%d%s.bmp", inpath.c_str(), prefix.c_str(), 16, postfix.c_str());
	sprintf_s(g_file, "%s/%s%d%s.bmp", inpath.c_str(), prefix.c_str(), 16, postfix.c_str());
	sprintf_s(r_file, "%s/%s%d%s.bmp", inpath.c_str(), prefix.c_str(), 16, postfix.c_str());
	//vector<int> capture_pentile_g_value({ 4,8,11,16,23,32,64 });		  
	//vector<int> capture_pentile_g_value({ 12,16,19,22,25,29,32 });	  
	const vector<int> capture_pentile_g_value({ 4,8,12,16,20,24,28,32 });
	for (int i = 0; i < capture_pentile_g_value.size(); i++)
	{
		sprintf_s(range_file[i], "%s/%s%d%s.bmp", inpath.c_str(), prefix.c_str(), 
			capture_pentile_g_value[i], postfix.c_str());
	}
	sprintf_s(mask_file, "%s/mask.png", inpath.c_str(), ms);
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
	sprintf_s(output_csv, "%s/pentile_rgb_relationship.csv", outpath.c_str());

	vector<vector<Point>> centers_error(3);
	vector<vector<vector<LED_info>>> relationship(3);
	{
		RGB select_rgb = BLUE;
		if (find_OLED_location_with_rgb_combination(rgb, mask, 
			outpath.c_str(), pentile_height, relationship))
		{
			cout << "something error!" << endl;
			system("pause");
			return 0;
		}
		cout << "--------------------" << endl;
		// to test 4,8,11,16,23,32,64 value
		vector<Mat> pic;
		{
			for (int i=0;i< capture_pentile_g_value.size();i++)
			{
				Mat tmp = imread(range_file[i]);
				if (tmp.data == NULL)
				{
					cout << range_file[i] << "read error." << endl;
					system("pause");
					return 0;
				}
				pic.push_back(tmp);
				cvtColor(*pic.rbegin(), *pic.rbegin(), COLOR_BGR2GRAY);
			}
		}
		compute_dumura(relationship, capture_pentile_g_value,
			pic, 3, outpath.c_str(), pentile_height, pentile_width / 2 * 3);
	}
#else
	String valid_path("./input2.2_20200404");
	char valid_b_file[MAX_PATH], valid_g_file[MAX_PATH], valid_r_file[MAX_PATH],
		//output
		valid_mask_file[MAX_PATH];
	int ms = 960, valid_ms = 1000;
	sprintf_s(valid_b_file, "%s/2号屏result_pentile_limit1-曝光时间750ms-增益3-施耐德镜头光圈5.7.bmp", valid_path.c_str(), valid_ms);
	sprintf_s(valid_g_file, "%s/2号屏result_pentile_limit1-曝光时间750ms-增益3-施耐德镜头光圈5.7.bmp", valid_path.c_str(), valid_ms);
	sprintf_s(valid_r_file, "%s/2号屏result_pentile_limit1-曝光时间750ms-增益3-施耐德镜头光圈5.7.bmp", valid_path.c_str(), valid_ms);
	sprintf_s(valid_mask_file, "%s/mask.png", valid_path.c_str(), valid_ms);
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
		tmp_valid_find_location(rgb, valid_mask,
			outpath.c_str(), pentile_height, relationship);
		//int  cnt1 =0, cnt2 = 0, cnt3 = 0;
		//double sum1 = 0, sum2 = 0, sum3 = 0;
		//vector<int> rect1({649, 2880, 512, 3429}),
		//	rect3({ 3843, 4440, 308, 247 }),
		//	rect2({ 10371, 6284, 546, 546 });
		//cout << rgb[1].type() << CV_8UC3;
		//for(auto pp : relationship[1])
		//	for (auto p : pp)
		//	{
		//		if(p.state==VALID && rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1]!= 255)
		//			/*if (p.pixel.x >= 0 && p.pixel.x <= 1115)
		//			{
		//				sum1 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt1++;
		//			}
		//			else if (p.pixel.x > 1115 && p.pixel.x <= 10211)
		//			{
		//				sum2 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt2++;
		//			}
		//			else if (p.pixel.x > 10211 && p.pixel.x <= 11229)
		//			{
		//				sum3 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt3++;
		//			}*/
		//			if (p.pixel.x >= rect1[0] && p.pixel.x <= rect1[0]+ rect1[2]
		//				&& p.pixel.y >= rect1[1] && p.pixel.y <= rect1[1] + rect1[3])
		//			{
		//				sum1 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt1++;
		//			}
		//			else if (p.pixel.x >= rect2[0] && p.pixel.x <= rect2[0] + rect2[2]
		//				&& p.pixel.y >= rect2[1] && p.pixel.y <= rect2[1] + rect2[3])
		//			{
		//				sum2 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt2++;
		//			}
		//			else if (p.pixel.x >= rect3[0] && p.pixel.x <= rect3[0] + rect3[2]
		//				&& p.pixel.y >= rect3[1] && p.pixel.y <= rect3[1] + rect3[3])
		//			{
		//				sum3 += rgb[1].at<Vec3b>(p.pixel.y, p.pixel.x)[1];
		//				cnt3++;
		//			}
		//	}
		//cout << cnt1 << " " << cnt2 << " " << cnt3 << endl;
		//cout << sum1/ cnt1 << " " << sum2/ cnt2 << " " << sum3/ cnt3 << endl;
		cout << "--------------------" << endl;
		compute_dumura_single_pic(relationship, gg, 
			outpath.c_str(), pentile_height, pentile_width / 2 * 3);
	}
#endif
	//system("pause");
	return 0;
}