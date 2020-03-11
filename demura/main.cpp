#include "location.h"
#include "gaussian_filter.h"
#include "preprocess.h"
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
/*
int xy, // x:0, y:1, no sigma compute:3
double from, double to, double another, double add, const char* prefix
*/
int main(int argc, char* argv[])
{
	//test_g32();
	_mkdir("./output");
	/*preprocess("./input2/5.85_G16.bmp", "./input2/5.85_G16.bmp", "./input2/5.85_G16.bmp",
		"./input2/mask.png");*/

	vector<Point> centers_error;
	//map<int, VectorXd> centers;
	vector<vector<Point>> centers_vec;
	vector<vector<VectorXd>> data;
	//find_OLED_location(centers_vec, data, centers_error);
	bool is_green = true;
	find_OLED_location_with_mask(centers_vec, data, centers_error, is_green);
	cout << "--------------------" << endl;
	if (argc >= 7)
	{
		cout << "demura..." << endl;
		compute_dumura(centers_vec, data, centers_error, atoi(argv[1]),
			atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6]);
	}
	system("pause");
	return 0;
}