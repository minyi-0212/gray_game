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
	Mat rgb = imread("./output_pentile/test_rgb.bmp"), pentile;
	rgb2pentile(rgb, pentile);
	imwrite("./output_pentile/pentile2.bmp", pentile);
	pentile2rgb(pentile, rgb);
	imwrite("./output_pentile/bask_test_rgb2.bmp", rgb);
}

void merge()
{
	Mat r = imread("./output_pentile/R16_result.png"), 
		g = imread("./output_pentile/B16_result.png"),
		b = imread("./output_pentile/B16_result.png"),
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
	imwrite("./output_pentile/merge.png", bgr);
}

/*
int xy, // x:0, y:1, no sigma compute:3
double from, double to, double another, double add, const char* prefix
*/
int main(int argc, char* argv[])
{
	//test_g32();
	_mkdir("./output_pentile");
	//draw_pattern2("./output_pentile/test2", 752, 2436);
	//test_pentile();
	//merge();
	//system("pause");
	//return 0;

	/*preprocess("./input2/5.85_B16.bmp", "./input2/5.85_B16.bmp", "./input2/5.85_B16.bmp",
		"./input2/mask.png");*/

	vector<Point> centers_error;
	//map<int, VectorXd> centers;
	vector<vector<Point>> centers_vec;
	vector<vector<VectorXd>> data;
	//for (int i = 12000; i >= 7000; i -= 1000)
	int i = 7000;
	{
		char path[MAX_PATH], inputfile[MAX_PATH], output_csv[MAX_PATH];
		sprintf_s(path, "E:/coding/gray_game/demura", i);
		sprintf_s(inputfile, "%s/input2_pentile/test_pentile_0_%dus.BMP", path, i);
		//sprintf_s(output_csv, "%s/output_pentile/%dus_pentile_rgb_relationship.csv",path, i);
		sprintf_s(output_csv, "%s/output_pentile/pentile_rgb_relationship.csv", path);
		//cout<<endl<<"write to "<< output_csv << endl;
		RGB select_rgb = BLUE;
		/*find_OLED_location_with_mask(
			inputfile, "E:/coding/gray_game/demura/input2_pentile/mask_r.bmp",
			"E:/coding/gray_game/demura/output_pentile/B16_selected_points.png", 
			"E:/coding/gray_game/demura/output_pentile/B16_after_3x3_set.png",
			centers_vec, data, centers_error, select_rgb == GREEN);*/
		/*find_OLED_location_with_mask(
			"./input2/5.85_B16.bmp", "./input2/mask.png",
			"./output_pentile/B16_selected_points.png", "./output_pentile/B16_after_3x3_set.png",
			centers_vec, data, centers_error, select_rgb == GREEN);*/
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