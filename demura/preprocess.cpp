#include "preprocess.h"
#include <fstream>

using namespace std;
using namespace cv;
using namespace Eigen;

void preprocess(vector<Mat>& rgb,
	int low_limit, const char *outfile, Mat& mask)
{
	const int ks = 3, sigma = 2, /*low_limit = 255 / 100,*/ erode_ks = 7;
	Mat img_b = rgb[BLUE].clone(),
		img_g = rgb[GREEN].clone(),
		img_r = rgb[RED].clone();
	cvtColor(img_b, img_b, CV_BGR2GRAY);
	cvtColor(img_g, img_g, CV_BGR2GRAY);
	cvtColor(img_r, img_r, CV_BGR2GRAY);
	mask = Mat(img_b.size(), CV_8UC1, Scalar(0));
	GaussianBlur(img_b, img_b, Size(ks, ks), sigma, sigma);
	GaussianBlur(img_g, img_g, Size(ks, ks), sigma, sigma);
	GaussianBlur(img_r, img_r, Size(ks, ks), sigma, sigma);
	cout << "rows*cols: " << img_b.rows << "*" << img_b.cols << endl;
	for (int i = 0; i < img_b.rows; i++)
	{
		for (int j = 0; j < img_b.cols; j++)
		{
			if (img_b.at<byte>(i, j) > low_limit
				&& img_g.at<byte>(i, j) > low_limit
				&& img_r.at<byte>(i, j) > low_limit)
			{
				mask.at<byte>(i, j) = 255;
			}
			else
			{ 
				mask.at<byte>(i, j) = 0;
			}
		}
	}
	Mat element = getStructuringElement(MORPH_RECT, Size(erode_ks, erode_ks));
	erode(mask, mask, element);
	//GaussianBlur(mask, mask, Size(5, 5), 2, 2);
	cout << "[output mask]: " << outfile << endl;
	imwrite(outfile, mask);
}