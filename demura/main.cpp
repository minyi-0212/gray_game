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
/*
int xy, // x:0, y:1, no sigma compute:3
double from, double to, double another, double add, const char* prefix
*/

void draw_box(int cx, int cy, int rx, int ry, Mat& p, Vec3b& color)
{
	int x = cx, y = cy;
	while (x < rx)
	{
		p.at<Vec3b>(y, x++) = color;
	}
	while (y < ry)
	{
		p.at<Vec3b>(y++, x) = color;
	}
	while (x > cx)
	{
		p.at<Vec3b>(y, x--) = color;
	}
	while (y > cy)
	{
		p.at<Vec3b>(y--, x) = color;
	}
}

void draw_pattern()
{
	/*Mat p(Size(2436, 752), CV_8UC3, Scalar(0, 0, 0));
	for (int y = 0; y<p.rows; )
	{
		for (int x = 0; x < p.cols; x++)
		{
			p.at<byte>(y, x) = x % 255;
		}
		if (y % 10 == 9)
			y += 11;
		else
			y++;
	}
	for (int x = 0; x < p.cols;)
	{
		for (int y = 0; y < p.rows; y++)
		{
			p.at<byte>(y, x) = 255;
		}
		if (x % 50 == 0)
			x++;
		else
			x += 49;
	}*/

	int cols = 2436, rows = 752, base = cols / 6;
	Mat p(Size(cols, rows), CV_8UC3, Scalar(0, 0, 0));
	byte tmp=0;
	for (int j = 3; j < rows-3; j+=2)
	{
		for (int i = 3; i < base; i++)
		{
			p.at<Vec3b>(j, i) = Vec3b(0, 0, tmp);
			p.at<Vec3b>(j, i + base) = Vec3b(0, tmp, 0);
			p.at<Vec3b>(j, i + base * 2) = Vec3b(tmp, 0, 0);
			tmp++;
		}
	}
	tmp = 0;
	for (int i = 3; i < base; i+=2)
	{
		for (int j = 3; j < rows-3; j++)
		{
			p.at<Vec3b>(j, i + base * 3) = Vec3b(0, 0, tmp);
			p.at<Vec3b>(j, i + base * 4) = Vec3b(0, tmp, 0);
			p.at<Vec3b>(j, i + base * 5) = Vec3b(tmp, 0, 0);
			tmp++;
		}
	}
	Vec3b white(255,255,255), black(0,0,0);
	draw_box(0, 0, cols-1, rows-1, p, white);
	draw_box(1, 1, cols-2, rows-2, p, black);
	draw_box(2, 2, cols-3, rows-3, p, black);
	imwrite("./output/pattern_rgb.png", p);

	ofstream out("./output/r.csv");
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			out << (int)p.at<Vec3b>(j, i)[2] << ",";
		}
		out << endl;
	}
	out.close();

	out.open("./output/g.csv");
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			out << (int)p.at<Vec3b>(j, i)[1] << ",";
		}
		out << endl;
	}
	out.close();

	out.open("./output/b.csv");
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			out << (int)p.at<Vec3b>(j, i)[0] << ",";
		}
		out << endl;
	}
	out.close();
}

void test_pentile()
{
	Mat rgb = imread("./output/test_rgb.bmp"), pentile;
	rgb2pentile(rgb, pentile);
	imwrite("./output/pentile2.bmp", pentile);
	pentile2rgb(pentile, rgb);
	imwrite("./output/bask_test_rgb2.bmp", rgb);
}

int main(int argc, char* argv[])
{
	//test_g32();
	_mkdir("./output");
	//draw_pattern();
	//test_pentile();

	/*preprocess("./input2/5.85_R16.bmp", "./input2/5.85_R16.bmp", "./input2/5.85_R16.bmp",
		"./input2/mask.png");*/

	//vector<Point> centers_error;
	////map<int, VectorXd> centers;
	//vector<vector<Point>> centers_vec;
	//vector<vector<VectorXd>> data;
	//bool is_green = false;
	//find_OLED_location_with_mask(centers_vec, data, centers_error, is_green);
	//cout << "--------------------" << endl;
	//if (argc >= 7)
	//{
	//	cout << "demura..." << endl;
	//	compute_dumura(centers_vec, data, centers_error, atoi(argv[1]),
	//		atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6]);
	//}

	system("pause");
	return 0;
}