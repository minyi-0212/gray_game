#include "pentile.h"
#include <fstream>
using namespace std;
using namespace cv;
using namespace Eigen;

const int r = RED, g = GREEN, b = BLUE;
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
	byte tmp = 0;
	for (int j = 3; j < rows - 3; j += 2)
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
	for (int i = 3; i < base; i += 2)
	{
		for (int j = 3; j < rows - 3; j++)
		{
			p.at<Vec3b>(j, i + base * 3) = Vec3b(0, 0, tmp);
			p.at<Vec3b>(j, i + base * 4) = Vec3b(0, tmp, 0);
			p.at<Vec3b>(j, i + base * 5) = Vec3b(tmp, 0, 0);
			tmp++;
		}
	}
	Vec3b white(255, 255, 255), black(0, 0, 0);
	draw_box(0, 0, cols - 1, rows - 1, p, white);
	draw_box(1, 1, cols - 2, rows - 2, p, black);
	draw_box(2, 2, cols - 3, rows - 3, p, black);
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

vector<int> cross_range({ -5,-4,-3,-2,-1,0,1,2,3,4,5 });
void draw_cross(Mat& img, int x, int y, Vec3b color)
{
	for (int xi = 0; xi < cross_range.size(); xi++)
		img.at<Vec3b>(y, x + cross_range[xi]) = color;
	for (int yi = 0; yi < cross_range.size(); yi++)
		img.at<Vec3b>(y + cross_range[yi], x) = color;
}

void draw_pattern2(const char* prefix, int penw, int penh)
{
	//void gen_test_image(std::vector<BYTE> &img, int &w, int &h)
	vector<BYTE> img;
	int w = penw / 2 * 3, h = penh;
	Mat rgb = Mat(Size(w, h), CV_8UC3, Scalar(0, 0, 0));
	img.resize(w * h * 3);
	Vec3b white(255, 255, 255),
		white_one(1, 1, 1),
		blue(255, 0, 0), green(0, 255, 0), red(0, 0, 255),
		blue_one(1, 0, 0), green_one(0, 1, 0), red_one(0, 0, 1);
	{
		// x±ß¿ò
		for (int y = 0; y <= 1; y++)
			for (int x = 0; x < w; x++)
				rgb.at<Vec3b>(y, x) = white;
		for (int y = h - 1; y >= h - 2; y--)
			for (int x = 0; x < w; x++)
				rgb.at<Vec3b>(y, x) = white;

		// y±ß¿ò
		for (int x = 0; x <= 1; x++)
			for (int y = 0; y < h; y++)
				rgb.at<Vec3b>(y, x) = white;
		for (int x = w - 1; x >= w - 5; x--)
			for (int y = 0; y < h; y++)
				rgb.at<Vec3b>(y, x) = white;
		/*for (int x = 20; x < w; x += 20)
			for (int y = 0; y < h; y++)
				if (y < 20 || y >= h - 20)
					rgb.at<Vec3b>(y, x) = white;*/

		int startx = 43 + 12;
		int starty = 130 + 12;
		int gapx = 12;
		int gapy = 12;
		int x0 = startx;
		int y0 = starty;
		/*{
			y0 += 100;
			int stride_width = 400;
			for (int y = 0; y < stride_width; y++)
				for (int x = 0; x < 256*4; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = (255-x/4)*white_one;

			y0 += stride_width + gapy;
			for (int y = 0; y < 400; y++)
				for (int x = 0; x < 256 * 4; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = x/4*white_one;

			y0 += stride_width + gapy + 340;
			stride_width = 500;
			for (int y = 0; y < 256*3; y++)
				for (int x = 0; x < stride_width; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = y / 3 * white_one;

			x0 += stride_width + gapx;
			for (int y = 0; y < 256 * 3; y++)
				for (int x = 0; x < stride_width; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = (255-y / 3) * white_one;
		}*/

		{
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = blue_one * x;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + 256 + gapx) = blue_one * (255 - x);
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 2) = blue_one * y;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 3) = blue_one * (255 - y);

			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256 * 4; x++)
					rgb.at<Vec3b>(y0 + 256 + gapy + y, x0 + x) = x / 4 * blue_one;
		}

		{
			y0 += (256 + gapy) * 2;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = green_one * x;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + 256 + gapx) = green_one * (255 - x);
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 2) = green_one * y;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 3) = green_one * (255 - y);

			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256 * 4; x++)
					rgb.at<Vec3b>(y0 + 256 + gapy + y, x0 + x) = x / 4 * green_one;
		}

		{
			y0 += (256 + gapy) * 2 + 20;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = red_one * x;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + 256 + gapx) = red_one * (255 - x);
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 2) = red_one * y;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 3) = red_one * (255 - y);

			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256 * 4; x++)
					rgb.at<Vec3b>(y0 + 256 + gapy + y, x0 + x) = x / 4 * red_one;
		}

		{
			y0 += (256 + gapy) * 2;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x) = red_one * x;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + 256 + gapx) = green_one * (255 - x);
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 2) = blue_one * y;
			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256; x++)
					rgb.at<Vec3b>(y0 + y, x0 + x + (256 + gapx) * 3) = 
					(x/3 > 60 ? red_one :(x/3>30?green_one :blue_one)) * (255 - y);

			for (int y = 0; y < 256; y++)
				for (int x = 0; x < 256 * 4; x++)
					rgb.at<Vec3b>(y0 + 256 + gapy + y, x0 + x) = 
					(255 - y) * (x / 3 > 240 ? red_one : (x / 3 > 120 ? green_one : blue_one));
		}
	}

	// cross
	{
		ofstream out("cross_point_xy.txt");
		vector<int> need_y, need_x({ 30,60,90, 530,560,590, 1030,1060,1090 });
		//for (int y = 30; y <= 150; y += 50)
		int y = 130;
		{
			need_y.push_back(y);
			for (int x = 30; x < w; x += 500)
			{
				draw_cross(rgb, x, y + 1, blue);
				draw_cross(rgb, x + 30, y, green);
				draw_cross(rgb, x + 60, y, red);
			}
		}
		//for (int y = h - 30; y >= h - 150; y -= 50)
		y = h - 130;
		{
			need_y.push_back(y);
			for (int x = 30; x < w; x += 500)
			{
				draw_cross(rgb, x, y + 1, blue);
				draw_cross(rgb, x + 30, y, green);
				draw_cross(rgb, x + 60, y, red);
			}
		}
		y = 1216;
		//for (int y = 1217 - 50; y <= 1217 + 50; y += 50)
		{
			need_y.push_back(y);
			for (int x = 30; x < w; x += 500)
			{
				draw_cross(rgb, x, y + 1, blue);
				draw_cross(rgb, x + 30, y, green);
				draw_cross(rgb, x + 60, y, red);
			}
		}
		//draw_cross(rgb, 100, 100, white);

		/*for (auto y : need_y)
			for (int x = 0; x < w; x++)
				if (x < 20 || x >= w - 20)
					rgb.at<Vec3b>(y, x) = white;
		for (auto x : need_x)
			for (int y = 0; y< h; y++)
				if (y < 20 || y >= h - 20)
					rgb.at<Vec3b>(y, x) = white;*/
		for (auto y : need_y)
			out << y << " ";
		out << endl;
		for (auto x : need_x)
			out << x << " ";
		out.close();
	}

	char outfile[MAX_PATH];
	Mat rgb_tmp = Mat(Size(w, h), CV_8UC3, Scalar(0, 0, 0));
	/*memcpy(rgb_tmp.data, img.data(), img.size() * sizeof(BYTE));
	sprintf_s(outfile, "%s_rgb_origin.bmp", prefix);
	imwrite(outfile, rgb_tmp);*/
	sprintf_s(outfile, "%s_rgb.bmp", prefix);
	rgb /= 8;
	imwrite(outfile, rgb);

	sprintf_s(outfile, "%s_pentile.bmp", prefix);
	Mat pentile;
	rgb2pentile(rgb, pentile);
	imwrite(outfile, pentile);

	rgb_tmp = Mat(Size(w, h), CV_8UC3, Scalar(0, 0, 16));
	rgb2pentile(rgb_tmp, pentile);
	sprintf_s(outfile, "%s_r.bmp", prefix);
	imwrite(outfile, pentile);

	rgb_tmp = Mat(Size(w, h), CV_8UC3, Scalar(0, 16, 0));
	rgb2pentile(rgb_tmp, pentile);
	sprintf_s(outfile, "%s_g.bmp", prefix);
	imwrite(outfile, pentile);

	rgb_tmp = Mat(Size(w, h), CV_8UC3, Scalar(16, 0, 0));
	rgb2pentile(rgb_tmp, pentile);
	sprintf_s(outfile, "%s_b.bmp", prefix);
	imwrite(outfile, pentile);
}