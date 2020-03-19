#include "pentile.h"
using namespace std;
using namespace cv;
using namespace Eigen;

const int r = 2, g = 1, b = 0;
void rgb2pentile(const Mat& rgb, Mat& pentile)
{
	//Mat rgb = imread("./output/rgb2.png");
	cout <<"rgb2pentile :: (rows*cols)"<< rgb.rows << " * " << rgb.cols << " => ";
	//Mat pentile(Size(rgb.cols / 3 * 2, rgb.rows), rgb.type(), Scalar(0, 0, 0));
	pentile = Mat(Size(rgb.cols / 3 * 2, rgb.rows), rgb.type(), Scalar(0, 0, 0));
	cout << pentile.rows << " * " << pentile.cols << endl;
	for (int y = 0; y < pentile.rows; y++)
	{
		for (int x = 0; x < pentile.cols>>1; x++)
		{
			if (((x + y) & 1) == 0)
			{
				pentile.at<Vec3b>(y, x * 2)[r] = rgb.at<Vec3b>(y, x * 3)[r]; // r
				pentile.at<Vec3b>(y, x * 2)[g] = rgb.at<Vec3b>(y, x * 3)[g]; // g
				pentile.at<Vec3b>(y, x * 2)[b] = rgb.at<Vec3b>(y, x * 3 + 1)[b]; // b
				pentile.at<Vec3b>(y, x * 2 + 1)[r] = rgb.at<Vec3b>(y, x * 3 + 1)[g]; // g
				pentile.at<Vec3b>(y, x * 2 + 1)[g] = rgb.at<Vec3b>(y, x * 3 + 2)[r]; // r
				pentile.at<Vec3b>(y, x * 2 + 1)[b] = rgb.at<Vec3b>(y, x * 3 + 2)[g]; // g
			}
			else
			{
				pentile.at<Vec3b>(y, x * 2)[r] = rgb.at<Vec3b>(y, x * 3)[b]; // r
				pentile.at<Vec3b>(y, x * 2)[g] = rgb.at<Vec3b>(y, x * 3)[g]; // g
				pentile.at<Vec3b>(y, x * 2)[b] = rgb.at<Vec3b>(y, x * 3 + 1)[r]; // b
				pentile.at<Vec3b>(y, x * 2 + 1)[r] = rgb.at<Vec3b>(y, x * 3 + 1)[g]; // g
				pentile.at<Vec3b>(y, x * 2 + 1)[g] = rgb.at<Vec3b>(y, x * 3 + 2)[b]; // r
				pentile.at<Vec3b>(y, x * 2 + 1)[b] = rgb.at<Vec3b>(y, x * 3 + 2)[g]; // g
			}
		}
	}
	//imwrite("./output/pattern_pentile2.png", pentile);
}

void pentile2rgb(const Mat& pentile, Mat& rgb)
{
	//Mat pentile = imread("./output/pentile2.bmp");
	cout << "pentile2rgb :: (rows*cols)" << pentile.rows << " * " << pentile.cols << " => ";
	//Mat rgb(Size(pentile.cols / 2 * 3, pentile.rows), pentile.type(), Scalar(0, 0, 0));
	rgb= Mat(Size(pentile.cols / 2 * 3, pentile.rows), pentile.type(), Scalar(0, 0, 0));
	cout << rgb.rows << " * " << rgb.cols << endl;
	for (int y = 0; y < pentile.rows; y++)
	{
		for (int x = 0; x < pentile.cols >> 1; x++)
		{
			if (((x + y) & 1) == 1)
			{
				rgb.at<Vec3b>(y, x * 3)[r] = pentile.at<Vec3b>(y, x * 2)[r]; // r
				rgb.at<Vec3b>(y, x * 3)[g] = pentile.at<Vec3b>(y, x * 2)[g]; // g
				rgb.at<Vec3b>(y, x * 3)[b] = pentile.at<Vec3b>(y, x * 2)[b];
				rgb.at<Vec3b>(y, x * 3 + 1)[b] = pentile.at<Vec3b>(y, x * 2)[b]; // b
				rgb.at<Vec3b>(y, x * 3 + 1)[g] = pentile.at<Vec3b>(y, x * 2 + 1)[r]; // g
				rgb.at<Vec3b>(y, x * 3 + 1)[r] = pentile.at<Vec3b>(y, x * 2 + 1)[g];
				rgb.at<Vec3b>(y, x * 3 + 2)[r] = pentile.at<Vec3b>(y, x * 2 + 1)[g]; // r
				rgb.at<Vec3b>(y, x * 3 + 2)[g] = pentile.at<Vec3b>(y, x * 2 + 1)[b]; // g
				rgb.at<Vec3b>(y, x * 3 + 2)[b] = pentile.at<Vec3b>(y, x * 2 + 2)[r];
			}
			else
			{
				rgb.at<Vec3b>(y, x * 3)[b] = pentile.at<Vec3b>(y, x * 2)[r]; // b
				rgb.at<Vec3b>(y, x * 3)[g] = pentile.at<Vec3b>(y, x * 2)[g]; // g
				rgb.at<Vec3b>(y, x * 3)[r] = pentile.at<Vec3b>(y, x * 2)[b];
				rgb.at<Vec3b>(y, x * 3 + 1)[r] = pentile.at<Vec3b>(y, x * 2)[b]; // r
				rgb.at<Vec3b>(y, x * 3 + 1)[g] = pentile.at<Vec3b>(y, x * 2 + 1)[r]; // g
				rgb.at<Vec3b>(y, x * 3 + 1)[b] = pentile.at<Vec3b>(y, x * 2 + 1)[g];
				rgb.at<Vec3b>(y, x * 3 + 2)[b] = pentile.at<Vec3b>(y, x * 2 + 1)[g]; // b
				rgb.at<Vec3b>(y, x * 3 + 2)[g] = pentile.at<Vec3b>(y, x * 2 + 1)[b]; // g
				rgb.at<Vec3b>(y, x * 3 + 2)[r] = pentile.at<Vec3b>(y, x * 2 + 2)[r];
			}
		}
	}

	/*int penh = pentile.rows, penw = pentile.cols,h = penh, w = penw / 2 * 3;
	vector<BYTE> penimg, img;
	penimg.resize(penw * penh * 3);
	img.resize(w * h * 3);
	
	penimg.assign((BYTE*)pentile.datastart, (BYTE*)pentile.dataend);
	for (int y = 0; y < penh; y++)
		for (int x = 0; x < penw / 2; x++)
		{
			const int penidx = (x + y * penw / 2);
			const int imgidx = (x + y * w / 3);

			BYTE r[3], g[3], b[3];
			for (int i = 0; i < 3; i++)
				r[i] = g[i] = b[i] = 0;

			if ((x + y) % 2 == 0)
			{
				r[0] = penimg[penidx * 6 + 0];
				g[0] = penimg[penidx * 6 + 1];
				b[1] = penimg[penidx * 6 + 2];
				g[1] = penimg[penidx * 6 + 3];
				r[2] = penimg[penidx * 6 + 4];
				g[2] = penimg[penidx * 6 + 5];
			}
			else {
				b[0] = penimg[penidx * 6 + 0];
				g[0] = penimg[penidx * 6 + 1];
				r[1] = penimg[penidx * 6 + 2];
				g[1] = penimg[penidx * 6 + 3];
				b[2] = penimg[penidx * 6 + 4];
				g[2] = penimg[penidx * 6 + 5];
			}

			img[imgidx * 9 + 0] = r[0]; //penimg[penidx*6+0];
			img[imgidx * 9 + 1] = g[0]; //penimg[penidx*6+1];
			img[imgidx * 9 + 2] = b[0];
			img[imgidx * 9 + 3] = r[1];

			img[imgidx * 9 + 4] = g[1]; //penimg[penidx*6+3];
			img[imgidx * 9 + 5] = b[1]; //penimg[penidx*6+2];
			img[imgidx * 9 + 6] = r[2]; //penimg[penidx*6+4];
			img[imgidx * 9 + 7] = g[2]; //penimg[penidx*6+5];
			img[imgidx * 9 + 8] = b[2];
		}
	memcpy(rgb.data, img.data(), img.size() * sizeof(BYTE));*/
	//imwrite("./output/rgb2.png", rgb);
}

void find_pentile_rgb_relationship()
{

}