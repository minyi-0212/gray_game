#include "location.h"
#include "gaussian_filter.h"
#include "preprocess.h"
#include "pentile.h"
#include <direct.h>
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

/*
int xy, // x:0, y:1, no sigma compute:3
double from, double to, double another, double add, const char* prefix
*/
int main(int argc, char* argv[])
{
	String inpath("E:/coding/gray_game/demura/input2.2.2_pentile"),
		outpath("E:/coding/gray_game/demura/output");
	_mkdir(outpath.c_str());
	int pentile_height = 2436;
	//draw_pattern2("./output/test2", 752, pentile_height);
	//test_pentile();
	//merge();
	//system("pause");
	//return 0;

	char b_file[MAX_PATH], g_file[MAX_PATH], r_file[MAX_PATH], 
		//output
		mask_file[MAX_PATH], cross_file[MAX_PATH],
		output_csv[MAX_PATH];
	int ms = 1500;
	sprintf_s(b_file, "%s/b-%dms.bmp", inpath.c_str(), ms);
	sprintf_s(g_file, "%s/g-%dms.bmp", inpath.c_str(), ms);
	sprintf_s(r_file, "%s/r-%dms.bmp", inpath.c_str(), ms);
	sprintf_s(cross_file, "%s/test2_pentile-35000us.bmp", inpath.c_str());
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
	//sprintf_s(b_file, "%s/b-%dms.png", inpath.c_str(), ms);
	//imwrite(b_file, rgb[b]);
	if(mask.data == nullptr)
		preprocess(rgb, 0, mask_file, mask);
	sprintf_s(output_csv, "%s/pentile_rgb_relationship.csv", outpath.c_str());

	//vector<Point> centers_error;
	////map<int, VectorXd> centers;
	//vector<vector<Point>> centers_vec;
	//vector<vector<VectorXd>> data;
	vector<vector<Point>> centers_error(3);
	vector<vector<vector<Point>>> centers_vec(3);
	vector<vector<vector<VectorXd>>> data(3);
	vector<vector<pair<Point,Point>>> cross_points(3);
	//for (int i = 12000; i >= 7000; i -= 1000)
	{
		//cout<<endl<<"write to "<< output_csv << endl;
		RGB select_rgb = BLUE;
		/*find_OLED_location_with_mask(
			inputfile, "E:/coding/gray_game/demura/input2_pentile/mask_r.bmp",
			"E:/coding/gray_game/demura/output/B16_selected_points.png", 
			"E:/coding/gray_game/demura/output/B16_after_3x3_set.png",
			centers_vec, data, centers_error, select_rgb == GREEN);*/
		/*find_OLED_location_with_mask(
			"./input2/5.85_B16.bmp", "./input2/mask.png",
			"./output/B16_selected_points.png", "./output/B16_after_3x3_set.png",
			centers_vec, data, centers_error, select_rgb == GREEN);*/
		/*find_OLED_location_with_mask(
			inputfile, "E:/coding/gray_game/demura/input2_pentile_2/mask.png",
			"E:/coding/gray_game/demura/output/selected_points.png",
			"E:/coding/gray_game/demura/output/after_3x3_set.png",
			centers_vec, data, centers_error, select_rgb == GREEN);*/
		
		//find_OLED_cross(b_file, cross, output_selected_points_prefix, centers_vec, data, centers_error);
		find_OLED_location_with_rgb_combination(rgb, mask, cross_file,
			outpath.c_str(), pentile_height, centers_vec, data, centers_error, cross_points);
		cout << "--------------------" << endl;

		/*if (argc >= 7)
		{
			cout << "demura..." << endl;
			compute_dumura(centers_vec, data, centers_error, atoi(argv[1]),
				atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6], select_rgb,
				output_csv);
		}*/
	}
	system("pause");
	return 0;
}