#include "diff.h"
using namespace std;
using namespace cv;
void difference(cv::Mat a, const cv::Mat b)
{
	int r = a.rows, c = a.cols;
	if (r != b.rows || c != b.rows)
	{
		cout << "in difference(): size not equal." << endl;
		return;
	}
	Vec3b red(0, 0, 1), green(0, 1, 0), blue(1, 0, 0);
	int tmp;
	for (int y = 0; y < r; y++)
	{
		for (int x = 0; x < c; x++)
		{
			tmp = a.at<Vec3b>(y, x)[1] - b.at<Vec3b>(y, x)[1];
			if (tmp > 0)
				a.at<Vec3b>(y, x) = red * tmp;
			else
				a.at<Vec3b>(y, x) = blue * (-tmp);
		}
	}
}