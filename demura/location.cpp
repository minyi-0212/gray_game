#include "location.h"
#include <fstream>
#include <stack>

using namespace std;
using namespace cv;
using namespace Eigen;

Mat img;
int IMG_WIDTH, IMG_HEIGHT;
//vector<int> t_update_small({ -2, -1, 0, 1, 2 });
//vector<int> t_update_large({ -2, -1, 0, 1, 2 });
//vector<int> t_small({ -2, -1, 0, 1, 2 });
//vector<int> t_large({ -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,-8 });

vector<int> t_update_small({ -1, 0, 1 });
vector<int> t_update_large({ -1, 0, 1 });
vector<int> t_small({ -2, -1, 0, 1, 2 });
vector<int> t_large({ -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8 });
vector<int> t_right_large({ -5,-4,-3,-2,-1,0,1,2,3,4,5 });
void update(Point& p)
{
	Point tmp(p);
	for (int i = 0; i < t_update_small.size(); i++)
	{
		for (int j = 0; j < t_update_large.size(); j++)
		{
			if (img.at<byte>(p.y, p.x)
				< img.at<byte>(tmp.y + t_update_large[j], tmp.x + t_update_small[i]))
			{

				p.y = tmp.y + t_update_large[j];
				p.x = tmp.x + t_update_small[i];
			}
			/*if (img.at<byte>(p.y, p.x) < 
				img.at<byte>(p.y + t_update_large[j], p.x + t_update_small[i]))
			{
				p.y += t_update_large[j];
				p.x += t_update_small[i];
			}*/
		}
	}

	/*for (int i = 0; i < t_small.size(); i++)
	{
		if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y, tmp.x + t_small[i]))
		{
			p.x = tmp.x + t_small[i];
		}
	}*/
}

void update(const Mat m,  Point& p)
{
	Point tmp(p);
	for (int i = 0; i < t_update_small.size(); i++)
	{
		for (int j = 0; j < t_update_large.size(); j++)
		{
			if (m.at<byte>(p.y, p.x)
				< m.at<byte>(tmp.y + t_update_large[j], tmp.x + t_update_small[i]))
			{

				p.y = tmp.y + t_update_large[j];
				p.x = tmp.x + t_update_small[i];
			}
		}
	}
}

void optimize_conner(vector<Point>& corner, int i)
{
	Point tmp = corner[i];
	//cout << corner[i] << " -> ";
	for (int ty = 0; ty < t_small.size(); ty++)
	{
		for (int tx = 0; tx < t_large.size(); tx++)
		{
			//cout << (int)img.at<byte>(corner[i].y + t_small[ty], corner[i].x + t_large[tx]) << " ";
			if (img.at<byte>(tmp.y, tmp.x)
				< img.at<byte>(corner[i].y + t_small[ty], corner[i].x + t_large[tx]))
			{
				tmp.x = corner[i].x + t_large[tx];
				tmp.y = corner[i].y + t_small[ty];
			}
		}
	}
	corner[i] = tmp;
	//cout << corner[i] << endl;

	tmp = corner[i+1];
	for (int ty = 0; ty < t_large.size(); ty++)
	{
		for (int tx = 0; tx < t_small.size(); tx++)
		{
			//cout << (int)img.at<byte>(corner[i].y + t19[ty], corner[i].x + t5[tx]) << " ";
			if (img.at<byte>(tmp.y, tmp.x)
				< img.at<byte>(corner[i+1].y + t_large[ty], corner[i+1].x + t_small[tx]))
			{
				tmp.x = corner[i+1].x + t_small[tx];
				tmp.y = corner[i+1].y + t_large[ty];
			}
		}
	}
	corner[i+1] = tmp;

	update(corner[i]);
	update(corner[i+1]);
	if ((corner[i].x - corner[i + 1].x)*(corner[i].x - corner[i + 1].x) +
		(corner[i].y - corner[i + 1].y)*(corner[i].y - corner[i + 1].y) < 8)
	{
		if (img.at<byte>(corner[i].y, corner[i].x)
			< img.at<byte>(corner[i+1].y, corner[i+1].x))
		{
			corner[i] = corner[i + 1];
		}
		else
		{
			corner[i+1] = corner[i];
		}
	}
}

void update_line_head_max(Point& p)
{
	Point tmp = p;
	for (int tx = 0; tx < t_right_large.size(); tx++)
	{
		if (img.at<byte>(p.y, p.x)
			< img.at<byte>(tmp.y, tmp.x + t_right_large[tx]))
		{
			p.x = tmp.x + t_right_large[tx];
		}
	}
}

void update_line_head_max(const Mat m, Point& p)
{
	Point tmp = p;
	for (int tx = 0; tx < t_right_large.size(); tx++)
	{
		if (m.at<byte>(p.y, p.x)
			< m.at<byte>(tmp.y, tmp.x + t_right_large[tx]))
		{
			p.x = tmp.x + t_right_large[tx];
		}
	}
}

//void find_OLED_location(map<int, VectorXd>& centers)
void find_OLED_location(vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error)
{
	Performance p;
	img = imread("./input2/5.85_B16.bmp");
	//imwrite("./input2/B16.png", img);
	Mat img_copy = img.clone();
	cvtColor(img, img, CV_BGR2GRAY);
	IMG_WIDTH = img.cols;
	IMG_HEIGHT = img.rows;
	int i0 = -1, i1, j0 = -1, j1;

	// find corner - not use
	//{
	//	Mat img_copy = img * 80;
	//	//GaussianBlur(img_copy, img_copy, Size(5, 5), 2, 2);
	//	/*for (int i = 0; i < img_copy.rows; i++)
	//	{
	//		for (int j = 0; j < img_copy.cols; j++)
	//		{
	//			if (img_copy.at<byte>(i, j) != 255
	//				|| img_copy.at<byte>(i, j)[1] != 255
	//				|| img_copy.at<byte>(i, j)[2] != 255)
	//			{
	//				img_copy.at<byte>(i, j) = 0;
	//				img_copy.at<byte>(i, j)[1] = 0;
	//				img_copy.at<byte>(i, j)[2] = 0;
	//			}
	//		}
	//	}*/
	//	int tmp = img_copy.cols / 2, i0_old = -1;
	//	for (int i = 0; i < img_copy.rows; i++)
	//	{
	//		if (img_copy.at<byte>(i, tmp) == Vec3b(255, 255, 255)
	//			/*&& img_copy.at<byte>(i, tmp + 1) == Vec3b(255, 255, 255)
	//			&& img_copy.at<byte>(i, tmp - 1) == Vec3b(255, 255, 255)*/)
	//		{
	//			if (i0 == -1) i0 = i;
	//			i1 = i;
	//		}
	//	}
	//	tmp = img_copy.rows / 2;
	//	for (int j = 0; j < img_copy.cols; j++)
	//	{
	//		if (img_copy.at<byte>(tmp, j) == Vec3b(255, 255, 255)
	//			/*&& img_copy.at<byte>(tmp - 1, j) == Vec3b(255, 255, 255)
	//			&& img_copy.at<byte>(tmp + 1, j) == Vec3b(255, 255, 255)*/)
	//		{
	//			if (j0 == -1) j0 = j;
	//			j1 = j;
	//		}
	//	}
	//	cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
	//	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
	//	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255,0,0), 1);
	//	imwrite("./output/G32_cpoy.png", img_copy);
	//	//waitKey(0);
	//}

	Performance pp;
	// find corner
	{
		const int low_limit = 5, initail_interval = 6;
		int tmp = img.cols / 2;
		for (int i = 0; i < img.rows; i++)
		{
			if (img.at<byte>(i, tmp) > low_limit
				|| img.at<byte>(i, tmp - initail_interval) > low_limit)
			{
				i0 = i;
				break;
			}
		}
		for (int i = img.rows - 1; i > 0; i--)
		{
			if (img.at<byte>(i, tmp) > low_limit
				|| img.at<byte>(i, tmp - initail_interval) > low_limit)
			{
				i1 = i;
				break;
			}
		}
		tmp = img.rows / 2;
		for (int j = 0; j < img.cols; j++)
		{
			if (img.at<byte>(tmp, j) > low_limit
				&&img.at<byte>(tmp - initail_interval, j) > low_limit)
			{
				j0 = j;
				break;
			}
		}
		for (int j = img.cols - 1; j > 0; j--)
		{
			if (img.at<byte>(tmp, j) > low_limit
				&&img.at<byte>(tmp - initail_interval, j) > low_limit)
			{
				j1 = j;
				break;
			}
		}
	}

	// computer the count of points in one line
	/*{
		Mat img_copy = img;
		float k = (1011 - 1009)*1.0 / (9603 - 337),
			b = 1009 - 337 * k;
		cout << k << " " << b << endl;
		cout << -b / k << endl;
		vector<Point> cnt_points;
		const int low_limit = 6;
		ofstream out("./output/one_row.csv");
		for (int j = 0; j < img_copy.cols; j++)
		{
			if (img_copy.at<byte>(k*j + b, j) > low_limit)
			{
				if (cnt_points.size() == 0
					|| j - cnt_points[cnt_points.size() - 1].x > 5)
				{
					Point p(j, k*j + b), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y, tmp.x + t_small[i]))
						{
							p.x = tmp.x + t_small[i];
						}
					}
					update(p);
					out << p.x << endl;
					cnt_points.push_back(p);
					img_copy.at<byte>(p.y, p.x) = Vec3b(0, 255, 0);
				}
			}
		}
		out.close();
		cout << "size of point in one line: " << cnt_points.size() << endl;
		imwrite("./output/G32_one_row.png", img_copy);
	}*/

	// computer the count of points in one bar
	/*{
		Mat img_copy = img;
		//float k = (6220 - 1008)*1.0 / (334 - 337),
		//	b = 1008 - 337 * k;
		float k = (6223 - 1012)*1.0 / (329 - 332),
			b = 1012 - 332 * k;
		cout << k << " " << b << endl;
		cout << -b / k << endl;
		vector<Point> cnt_points;
		ofstream out("./output/one_column.csv");
		for (int y = 0; y < img_copy.rows; y++)
		{
			if (img_copy.at<byte>(y, (y - b) / k) > 7)
			{
				if (cnt_points.size() == 0
					|| y - cnt_points[cnt_points.size() - 1].y > 5)
				{
					//cnt_points.push_back(Point((y - b) / k, y));
					//img_copy.at<byte>(y, (y - b) / k) = Vec3b(0, 255, 0);
					Point p((y - b) / k, y), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y + t_small[i], tmp.x))
						{
							p.y = tmp.y + t_small[i];
						}
					}
					update(p);
					out << p.y << endl;
					cnt_points.push_back(p);
					img_copy.at<byte>(p.y, p.x) = Vec3b(0, 255, 0);
				}
			}
		}
		out.close();
		cout << "size of point in one line: " << cnt_points.size() << endl;
		imwrite("./output/G32_one_column.png", img_copy * 8);
	}*/

	// draw corner
	/*cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
	Mat img_copy = img;
	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255, 0, 0), 1);
	imwrite("./output/G32_cpoy.png", img_copy);*/

	// draw corner
	/*cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
	Mat img_copy = img;
	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255, 0, 0), 1);
	imwrite("./output/G32_cpoy.png", img_copy);*/

	// get more detailed corner
	/*  0 бнбн 2
	   1  бнбн  3

	   5  бнбн  7
		4 бнбн 6
	*/
	vector<Point> corner({
		Point(j0, i0), Point(j0, i0),
		Point(j1, i0), Point(j1, i0),
		Point(j0, i1), Point(j0, i1),
		Point(j1, i1), Point(j1, i1)
		});
	for (int i = 0; i < corner.size(); i += 2)
	{
		optimize_conner(corner, i);
		cout << corner[i] << endl;
		cout << corner[i+1] << endl;
	}

	// draw corner
	/*Mat img_copy = img;
	line(img_copy, corner[0], corner[2], Scalar(0, 0, 255), 1);
	line(img_copy, corner[1], corner[5], Scalar(0, 0, 255), 1);
	line(img_copy, corner[4], corner[6], Scalar(0, 0, 255), 1);
	line(img_copy, corner[7], corner[3], Scalar(0, 0, 255), 1);
	for (int i = 0; i < corner.size(); i += 2)
	{
		img_copy.at<byte>(corner[i].y, corner[i].x) = Vec3b(255, 0, 255);
		img_copy.at<byte>(corner[i+1].y, corner[i+1].x) = Vec3b(0, 255, 255);
	}
	imwrite("./output/G32_corner.png", img_copy);*/
	pp.endAndPrint();

	//1

//	{
//		cout << "-----------------" << endl;
//		vector<Point> centers;
//		Point begin = corner[0], end_line = corner[2], end_bar = corner[4];
//		float real_value = 255.0 / 32,
//			addi = (end_bar.y - 5 - begin.y) / 720.0,
//			addj = (end_line.x - begin.x) / 960.0;
//		Mat img_copy = img * real_value;
////#pragma parallel omp for
////		for (float i = corner[0].y; i <= corner[2].y; i += addi)
////		//for (float i = i0; i < i1 + 3; i += addi)
////		//int i = i0;
////		{
////			//cout << (int)((i-i0)/addi + 1) % 2 * addj / 2 << endl;
////			//for (float j = j0 + (int)((i - i0) / addi + 1) % 2 * addj; j <= j1; j += addj)
////			for (float j = corner[1].x; j <= corner[1].x; j += addj)
////			{
////				Point p(j, i);
////				img_copy.at<byte>(i, j) = Vec3b(0, 255, 255);
////				update(p);
////				//update(p);
////				centers.push_back(p);
////			}
////		}
//		cout << corner[0] << endl
//			<< corner[2] << endl
//			<< corner[4] << endl;
//		float k_line = (begin.y - end_line.y)*1.0 / (begin.x - end_line.x),
//			k_bar = (begin.y - end_bar.y + 5)*1.0 / (begin.x - end_bar.x + 5);
//		cout << k_line << " " << k_bar << endl;
////#pragma parallel omp for
//		for (float i = begin.y; i <= end_bar.y; i += addi)
////		int i = begin.y + addi;
//		{
//			//for (float j = begin.x; j <= end_line.x+1; j += addj)
//			float j = begin.x;
//			{
//				//float k = (1011 - 1009)*1.0 / (9603 - 337);
//				//float b = 1009 - 337 * k;
//				//cout << k << endl;
//
//				float b = i - begin.x * k_line;
//				Point p(j, k_line*j+b);
//				img_copy.at<byte>(p.y, p.x) = Vec3b(0, 255, 255);
//				update(p);
//				centers.push_back(p);
//			}
//		}
//		p.endAndPrint();
//
//		cout << "cornor index : " << i0 << " " << j0 << " " << i1 << " " << j1 << endl;
//		cout << "internal: " << addi << " " << addj << endl;
//		cout << "point size: " << endl << 1920 * 1080 << endl << centers.size() << endl;
//		/*img_copy.at<byte>(i0, j0) = Vec3b(0, 255, 0);
//		img_copy.at<byte>(i1, j1) = Vec3b(0, 255, 0);
//		img_copy.at<byte>(i1, j0) = Vec3b(0, 255, 0);
//		img_copy.at<byte>(i0, j1) = Vec3b(0, 255, 0);*/
//		for (auto p : centers)
//		{
//			img_copy.at<byte>(p.y, p.x) = Vec3b(0, 0, 255);
//		}
//
//		imwrite("./output/G32_position.png", img_copy);
//	}

	// 2

//	{
//		//cout << "-----------------" << endl;
//		//vector<Point> centers, line_head;
//		vector<Point> line_head;
//		{
//			float k1 = (corner[1].y - corner[5].y)*1.0 / (corner[1].x - corner[5].x),
//				b1 = corner[1].y - corner[1].x * k1,
//				b2 = corner[0].y - corner[0].x * k1,
//				add = (corner[5].y - corner[1].y) / 720.0;
//			//cout << "add in y:"<< add << endl;
//			for (float y1 = corner[1].y, y2 = corner[0].y;
//				y1 <= corner[5].y, y2 <= corner[6].y; y1 += add, y2 += add)
//			{
//				line_head.push_back(Point((y2 - b2) / k1, y2));
//				line_head.push_back(Point((y1 - b1) / k1, y1));
//			}
//		}
//		{
//			float k = (corner[0].y - corner[2].y)*1.0 / (corner[0].x - corner[2].x),
//				add = (corner[2].x - corner[0].x) / 960.0;
//			//cout << "add in x:" << add << endl;
//#pragma parallel omp for
//			for (auto p : line_head)
//			{
//				for (int x = p.x; x < corner[2].x; x += add)
//				{
//					float b = p.y - p.x*k;
//					Point p(x, x*k + b);
//					update(p);
//					update(p);
//					centers.push_back(p);
//				}
//			}
//		}
//		p.endAndPrint();
//
//		float real_value = 255.0 / 32;
//		Mat img_copy = img * real_value;
//		for (auto p : line_head)
//		{
//			img_copy.at<byte>(p.y, p.x) = Vec3b(0, 255, 255);
//		}
//		for (auto p : centers)
//		{
//			img_copy.at<byte>(p.y, p.x) = Vec3b(0, 0, 255);
//		}
//		//imwrite("./output/G32_position.png", img_copy);
//	}

	// 3
	{
		vector<Point> line_head;
		{
			vector<vector<Point>> line_head_two(2);
			const int head_low_limit = 5, head_initail_interval = 6;
			float k = (corner[1].y - corner[5].y)*1.0 / (corner[1].x - corner[5].x),
				b1 = corner[1].y - corner[1].x * k,
				b2 = corner[0].y - corner[0].x * k;
			/*for (float y = corner[1].y; y <= corner[5].y; y++)
			{
				if (img.at<byte>(y, (y - b1) / k) > head_low_limit &&
					(line_head.size() == 0 || y - line_head[line_head.size() - 1].y > head_initail_interval))
				{
					Point p((y - b1) / k, y), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y + t_small[i], tmp.x))
						{
							p.y = tmp.y + t_small[i];
						}
					}
					update(p);
					line_head.push_back(p);
				}
			}
			int tmp = line_head.size();
			for (float y = corner[0].y; y <= corner[4].y; y++)
			{
				if (img.at<byte>(y, (y - b2) / k) > head_low_limit &&
					(line_head.size() == tmp || y - line_head[line_head.size() - 1].y > head_initail_interval))
				{
					Point p((y - b2) / k, y), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y + t_small[i], tmp.x))
						{
							p.y = tmp.y + t_small[i];
						}
					}
					update(p);
					line_head.push_back(p);
				}
			}*/
			for (float y = corner[1].y; y <= corner[5].y; y++)
			{
				if (img.at<byte>(y, (y - b1) / k) > head_low_limit &&
					(line_head_two[0].size() == 0 
						|| y - line_head_two[0][line_head_two[0].size() - 1].y > head_initail_interval))
				{
					Point p((y - b1) / k, y), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y + t_small[i], tmp.x))
						{
							p.y = tmp.y + t_small[i];
						}
					}
					update(p);
					line_head_two[0].push_back(p);
				}
			}
			for (float y = corner[0].y; y <= corner[4].y; y++)
			{
				if (img.at<byte>(y, (y - b2) / k) > head_low_limit &&
					(line_head_two[1].size() == 0 
						|| y - line_head_two[1][line_head_two[1].size() - 1].y > head_initail_interval))
				{
					Point p((y - b2) / k, y), tmp(p);
					for (int i = 0; i < t_small.size(); i++)
					{
						if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y + t_small[i], tmp.x))
						{
							p.y = tmp.y + t_small[i];
						}
					}
					update(p);
					line_head_two[1].push_back(p);
				}
			}
			if (line_head_two[0].size() != line_head_two[1].size())
			{
				cout << "something error :: find the head point, with size "
					<< line_head_two[0].size()<< "!=" << line_head_two[1].size() << endl;
			}
			line_head.resize(line_head_two[0].size()*2);
			for (int i = 0; i < line_head_two[0].size(); i++)
			{
				line_head[2 * i] = line_head_two[0][i];
				line_head[2 * i + 1] = line_head_two[1][i];
			}
		}

		centers_vec.resize(line_head.size());
		data.resize(line_head.size());
		set<int> center_index_set;
		{
			float k = (corner[0].y - corner[2].y)*1.0 / (corner[0].x - corner[2].x);
			VectorXd center_region(9), center_error_region(9);
			center_error_region << -1, -1, -1, -1, -1, -1, -1, -1, -1;
			Point last;
			int tmp_interval;
			for (int hi = 0; hi < line_head.size(); hi++)
			{
				//centers.push_back(line_head[hi]);
				center_region << img.at<byte>(line_head[hi].y - 1, line_head[hi].x - 1),
					img.at<byte>(line_head[hi].y - 1, line_head[hi].x),
					img.at<byte>(line_head[hi].y - 1, line_head[hi].x + 1),
					img.at<byte>(line_head[hi].y, line_head[hi].x - 1),
					img.at<byte>(line_head[hi].y, line_head[hi].x),
					img.at<byte>(line_head[hi].y, line_head[hi].x + 1),
					img.at<byte>(line_head[hi].y + 1, line_head[hi].x - 1),
					img.at<byte>(line_head[hi].y + 1, line_head[hi].x),
					img.at<byte>(line_head[hi].y + 1, line_head[hi].x + 1);
				//centers[line_head[hi].y*IMG_WIDTH + line_head[hi].x] = center_region;
				//centers.insert({ line_head[hi].y*IMG_WIDTH + line_head[hi].x, center_region });
				center_index_set.insert(line_head[hi].y*IMG_WIDTH + line_head[hi].x);
				centers_vec[hi].push_back(line_head[hi]);
				data[hi].push_back(center_region);

				last = line_head[hi];
				tmp_interval = 8;
				for (int x = line_head[hi].x++; x < corner[3].x + 2; x++)
				{
					float b = line_head[hi].y - line_head[hi].x*k;
					if (/*img.at<byte>(k*x + b, x) > low_limit
						&&*/ x - last.x >= tmp_interval)
					{
						Point p(x, k*x + b), tmp(p);
						for (int i = 0; i < t_small.size(); i++)
						{
							if (img.at<byte>(p.y, p.x) < img.at<byte>(tmp.y, tmp.x + t_small[i]))
							{
								p.x = tmp.x + t_small[i];
							}
						}
						update(p);
						update(p);
						update(p);
						//if (img.at<byte>(p.y, p.x) > 12)
						if (img.at<byte>(p.y, p.x) < 12
							|| p.x - last.x < 6 || p.x - last.x > 11
							|| p.y - last.y > 2
							//|| centers.find(p.y*IMG_WIDTH + p.x) != centers.end())
							|| center_index_set.find(p.y*IMG_WIDTH + p.x) != center_index_set.end())
						{
							centers_error.push_back(tmp);
							/*cout << tmp << " , " << p <<" , " << last << " , "
								<< (p.x - last.x < 6) <<","
								<< (p.x - last.x > 11) << ","
								<< (p.y - last.y > 2) << ","
								<< (centers.find(p.y*IMG_WIDTH + p.x) != centers.end())<< " "
								<< tmp_interval  << endl;*/

							centers_vec[hi].push_back(tmp);
							data[hi].push_back(center_error_region);

							p = tmp;
						}
						else
						{
							center_region << img.at<byte>(p.y - 1, p.x - 1),
								img.at<byte>(p.y - 1, p.x),
								img.at<byte>(p.y - 1, p.x + 1),

								img.at<byte>(p.y, p.x - 1),
								img.at<byte>(p.y, p.x),
								img.at<byte>(p.y, p.x + 1),

								img.at<byte>(p.y + 1, p.x - 1),
								img.at<byte>(p.y + 1, p.x),
								img.at<byte>(p.y + 1, p.x + 1);
							//centers.insert({ p.y*IMG_WIDTH + p.x, center_region });
							center_index_set.insert(p.y*IMG_WIDTH + p.x);
							centers_vec[hi].push_back(p);
							data[hi].push_back(center_region);
							tmp_interval = abs(p.x - last.x - tmp_interval) < 2 ? p.x - last.x : tmp_interval;
						}
						//centers.insert({ p.y*IMG_WIDTH + p.x, center_region });
						last = p;
					}
				}
			}
		}
		p.endAndPrint();

		/*int cnt = 0;
		for (int i = 0; i < centers.size(); i++)
		{
			for (int j = 0; j < centers.size(); j++)
				if (i != j && (abs(centers[j].x - centers[i].x)<=3) 
					&& (abs(centers[i].y - centers[j].y)<=3))
				{
					if(abs(centers[i].y - centers[j].y) != 0)
						cout << i<<","<<j<<" " << centers[i] <<","<<centers[j]<< endl;
					cnt++;
				}
		}
		cout <<"again cnt: "<< cnt << endl;*/
		
		for (int vj=0;vj<centers_vec.size();vj++)
		{
			for (int vi = 0; vi < centers_vec[vj].size(); vi++)
			{
				img_copy.at<Vec3b>(centers_vec[vj][vi].y, centers_vec[vj][vi].x) = Vec3b(0, 255, 0);
			}
		}
		for (auto ce : centers_error)
		{
			img_copy.at<Vec3b>(ce.y, ce.x) = Vec3b(0, 0, 255);
		}
		cout << "size of head points: " << line_head.size() << endl;
		cout << 960 * 720 * 2 << ", size of points: " << center_index_set.size() << endl;
		cout << "size of error points: " << centers_error.size() << endl;
		imwrite("./output/G32_pick_point.png", img_copy * 8);

		img_copy = Mat(img_copy.size(), CV_8UC3, Scalar(0,0,0));
		vector<int> region({ -1,0,1 });
		for (int vj = 0; vj < centers_vec.size(); vj++)
		{
			for (int vi = 0; vi < centers_vec[vj].size(); vi++)
			{
				for (int i = 0; i < region.size(); i++)
				{
					for (int j = 0; j < region.size(); j++)
					{
						int x = centers_vec[vj][vi].x + region[j],
							y = centers_vec[vj][vi].y + region[i];
						byte t = img.at<byte>(y, x);
						img_copy.at<Vec3b>(y, x) = Vec3b(t, t, t);
					}
				}
			}
		}
		for (auto ce : centers_error)
		{
			for (int i = 0; i < region.size(); i++)
			{
				for (int j = 0; j < region.size(); j++)
				{
					byte t = img.at<byte>(ce.y, ce.x);
					img_copy.at<Vec3b>(ce.y, ce.x) = Vec3b(t, t, t);
				}
			}
		}
		imwrite("./output/G32_after_set.png", img_copy);
		//cvtColor(img_copy, img_copy, CV_BGR2GRAY);
		//imwrite("./output/G32_after_set.bmp", img_copy);
	}
}

void find_OLED_location_with_mask(
	const char *infile_name, const char *mask_file_name, 
	const char *outfile_selected_point_name, const char *outfile_set_3x3_region_name,
	vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error, bool green)
{
	Performance p;
	img = imread(infile_name);
	Mat mask = imread(mask_file_name, CV_8UC1);
	if (img.data == NULL || mask.data == NULL)
	{
		cout << "image or mask read error." << endl;
		return ;
	}
	imwrite("./input2_pentile/real.png", img);
	Mat img_copy = img.clone();
	cvtColor(img, img, CV_BGR2GRAY);
	IMG_WIDTH = img.cols;

	Performance pp;
	vector<Point> line_head;
	{
		Point p_start(-1,-1);
		const int head_low_limit = 5;
		int head_initail_interval_x = 4, head_initail_interval_y = 4;
		// find first point
		for (int y = 0; y < mask.rows; y++)
		{
			for (int x = 0; x < mask.cols; x++)
			{
				if (mask.at<byte>(y, x) > 0)
				{
					p_start.x = x;
					p_start.y = y;
					break;
				}
			}
			if (p_start.x != -1)
				break;
		}
		//cout << p_start << endl;
		update(p_start);
		for (int x = 0; x < mask.cols; x++)
		{
			if (mask.at<byte>(p_start.y, x) > 0
				&& img.at<byte>(p_start.y, x) > head_low_limit)
			{
				p_start.x = x;
				break;
			}
		}
		update_line_head_max(p_start);
		update(p_start);
		//cout << p_start << endl;
		//line_head.push_back(p_start);

		// find head point
		int i1; // , tmp = img.cols / 2;
		for (int i = img.rows - 1; i > 0; i--)
		{
			if (mask.at<byte>(i, p_start.x) > 0 
				&& (img.at<byte>(i, p_start.x) > head_low_limit
				|| img.at<byte>(i, p_start.x - head_initail_interval_x) > head_low_limit))
			{
				i1 = i;
				break;
			}
		}
		cout <<"line range: "<< p_start.y <<"->"<< i1 << endl;
		int xx = p_start.x;
		while(p_start.y < i1)
		{
			if(line_head.size()==0 || p_start.y - line_head[line_head.size()-1].y > 3)
				line_head.push_back(p_start);
			if(green)
				p_start.x = xx;
			else
			{
				p_start.x += head_initail_interval_x;
				head_initail_interval_x *= -1;
			}
			p_start.y += head_initail_interval_y;
			update_line_head_max(p_start);
			update(p_start);
		}
	}
	{
		stack<Point> stack_of_points;
		Point cur, last;
		const int low_limit = 3;
		int  interval_x = 18; //10;
		if (green)
			interval_x = 4;
		centers_vec.resize(line_head.size());
		data.resize(line_head.size());
		for (int hi = 0; hi < line_head.size(); hi++)
		//int hi = 0;
		{
			VectorXd center_region(9), center_error_region(9);
			center_error_region << -1, -1, -1, -1, -1, -1, -1, -1, -1;
			// left
			cur = line_head[hi];
			stack_of_points.push(cur);
			while (mask.at<byte>(cur.y, cur.x) > 0
				/*&& mask.at<byte>(cur.y, cur.x - interval_x) > 0*/)
			{
				last = cur;
				cur.x -= interval_x;
				if (mask.at<byte>(cur.y, cur.x) == 0)
					break;
				update(cur);
				update(cur);
				if (mask.at<byte>(cur.y, cur.x) == 0)
					break;
				stack_of_points.push(cur);
				if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x -= interval_x;
				}
			}
			while (!stack_of_points.empty())
			{
				last = cur;
				cur = stack_of_points.top();
				stack_of_points.pop();
				if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x += interval_x;
					centers_error.push_back(cur);

					centers_vec[hi].push_back(cur);
					data[hi].push_back(center_error_region);
				}
				else
				{
					centers_vec[hi].push_back(cur);
					center_region << img.at<byte>(cur.y - 1, cur.x - 1),
						img.at<byte>(cur.y - 1, cur.x),
						img.at<byte>(cur.y - 1, cur.x + 1),

						img.at<byte>(cur.y, cur.x - 1),
						img.at<byte>(cur.y, cur.x),
						img.at<byte>(cur.y, cur.x + 1),

						img.at<byte>(cur.y + 1, cur.x - 1),
						img.at<byte>(cur.y + 1, cur.x),
						img.at<byte>(cur.y + 1, cur.x + 1);
					data[hi].push_back(center_region);
				}
			}

			// right
			cur = line_head[hi];
			while (1 /*mask.at<byte>(cur.y, cur.x) > 0
				&& mask.at<byte>(cur.y, cur.x + interval_x) > 0*/)
			{
				last = cur;
				cur.x += interval_x;
				if (mask.at<byte>(cur.y, cur.x) == 0) {
					break;
				}
				update(cur);
				update(cur);
				if (mask.at<byte>(cur.y, cur.x) == 0) {
					break;
				}
				if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x += interval_x;
					centers_error.push_back(cur);

					centers_vec[hi].push_back(cur);
					data[hi].push_back(center_error_region);
				}
				else
				{
					centers_vec[hi].push_back(cur);
					center_region << img.at<byte>(cur.y - 1, cur.x - 1),
						img.at<byte>(cur.y - 1, cur.x),
						img.at<byte>(cur.y - 1, cur.x + 1),

						img.at<byte>(cur.y, cur.x - 1),
						img.at<byte>(cur.y, cur.x),
						img.at<byte>(cur.y, cur.x + 1),

						img.at<byte>(cur.y + 1, cur.x - 1),
						img.at<byte>(cur.y + 1, cur.x),
						img.at<byte>(cur.y + 1, cur.x + 1);
					data[hi].push_back(center_region);
				}
			}
		}
	}
	pp.endAndPrint();

	{
		cout << "head size: " << line_head.size() << endl;
		int sum = 0;
		cout << centers_vec.size() << "*" << centers_vec[510].size() << endl;
		for (auto p : centers_vec)
		{
			sum += p.size();
			for (auto pp : p)
			{
				img_copy.at<Vec3b>(pp.y, pp.x) = Vec3b(0, 255, 0);
			}
		}
		for (auto p : line_head)
			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 0, 255);
		cout << "total ceneter points size : " << sum << endl;
		cout << "centers_error size : " << centers_error.size() << endl;
		for (auto p : centers_error)
		{
			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
			//cout << p << endl;
		}
		imwrite(outfile_selected_point_name, img_copy);
	}
	{
		img_copy = Mat(img_copy.size(), CV_8UC3, Scalar(0, 0, 0));
		vector<int> region({ -1,0,1 });
		for (int vj = 0; vj < centers_vec.size(); vj++)
		{
			for (int vi = 0; vi < centers_vec[vj].size(); vi++)
			{
				for (int i = 0; i < region.size(); i++)
				{
					for (int j = 0; j < region.size(); j++)
					{
						int x = centers_vec[vj][vi].x + region[j],
							y = centers_vec[vj][vi].y + region[i];
						byte t = img.at<byte>(y, x);
						img_copy.at<Vec3b>(y, x) = Vec3b(t, t, t);
					}
				}
			}
		}
		imwrite(outfile_set_3x3_region_name, img_copy);
	}
	p.endAndPrint();
}

void get_box_corner_points(const Mat& mask, vector<Point>& corners, int low_limit)
{
	int ys, ye, xs, xe;
	int tmp = mask.cols / 2;
	for (int y = 0; y < mask.rows; y += 200)
	{
		if (mask.at<byte>(y, tmp) > low_limit)
		{
			ys = y;
			break;
		}
	}
	for (int y = mask.rows - 1; y > 0; y -= 200)
	{
		if (mask.at<byte>(y, tmp) > low_limit)
		{
			ye = y;
			break;
		}
	}
	tmp = mask.rows / 4 * 3;
	for (int x = 0; x < mask.cols; x += 200)
	{
		if (mask.at<byte>(tmp, x) > low_limit)
		{
			xs = x;
			break;
		}
	}
	for (int x = mask.cols - 1; x > 0; x -= 200)
	{
		if (mask.at<byte>(tmp, x) > low_limit)
		{
			xe = x;
			break;
		}
	}
	corners.resize(9);
	//// 0 1
	//// 2 3
	//corners[0].x = corners[2].x = xs;
	//corners[1].x = corners[3].x = xe;
	//corners[0].y = corners[1].y = ys;
	//corners[2].y = corners[3].y = ye;
	
	// 0 1 2
	// 3 4 5
	// 6 7 8
	
	corners[0].x = corners[3].x = corners[6].x = xs;
	corners[1].x = corners[4].x = corners[7].x = (xe + xs) / 2;
	corners[2].x = corners[5].x = corners[8].x = xe;
	corners[0].y = corners[1].y = corners[2].y = ys;
	corners[3].y = corners[4].y = corners[5].y = (ye + ys) / 2;
	corners[6].y = corners[7].y = corners[8].y = ye;
}

void find_OLED_cross(const char *img_file, const char *cross_file,
	const char *output_select_point_file,
	vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error)
{
	Performance p;
	Performance pp;
	img = imread(img_file);
	Mat cross = imread(cross_file),
		tmp, mask(img.size(), CV_8UC1);
	Mat img_copy = img.clone();
	cvtColor(img, img, CV_BGR2GRAY);
	if (img.data == NULL || cross.data == NULL)
	{
		cout << "image or cross file read error." << endl;
		return;
	}
	
	const int cross_low_limit = 30, select_range = 500, select_interval = 9;

	// find corner 0,1; 2,3
	vector<Point> corners;
	{
		const int ks = 7, sigma = 3;
		GaussianBlur(img, tmp, Size(ks, ks), sigma, sigma);
		/*for (int i = 0; i < tmp.rows; i++)
		{
			for (int j = 0; j < tmp.cols; j++)
			{
				if (tmp.at<byte>(i, j) > scale)
				{
					mask.at<byte>(i, j) = 255;
				}
				else
				{
					mask.at<byte>(i, j) = 0;
				}
			}
		}*/
		threshold(tmp, mask, 10, 255, THRESH_BINARY);
		get_box_corner_points(mask, corners, 0);
		tmp = cross.clone();
		for (auto p : corners)
		{
			//tmp.at<Vec3b>(range[0], range[2]) = Vec3b(0, 255, 0);
			//tmp.at<Vec3b>(range[1], range[3]) = Vec3b(0, 255, 0);
			circle(tmp, p, 2, Scalar(0, 255, 0), 2);
			//circle(tmp, p, 500, Scalar(0, 255, 0), 2);
		}
	}
	
	// find cross
	vector<vector<Point>> cross_points(corners.size());
	{
		int x0 = corners[0].x, y0 = corners[0].y;
		Point p_tmp;
		for(int ci=0;ci< corners.size();ci++)
			for (int y = corners[ci].y - select_range; y <= corners[ci].y + select_range; y++)
			{
				for (int x = corners[ci].x - select_range; x <= corners[ci].x + select_range; x++)
				{
					if (cross.at<Vec3b>(y, x)[0] > cross_low_limit
						&& cross.at<Vec3b>(y + select_interval, x)[0] > cross_low_limit
						&& cross.at<Vec3b>(y, x + select_interval)[0] > cross_low_limit
						&& cross.at<Vec3b>(y - select_interval, x)[0] > cross_low_limit
						&& cross.at<Vec3b>(y, x - select_interval)[0] > cross_low_limit
						&& cross.at<Vec3b>(y + select_interval, x + select_interval)[0] < cross_low_limit /2)
					{
						//circle(tmp, Point(x0,y0), 2, Scalar(0, 255, 0), 2);
						p_tmp.x = x;
						p_tmp.y = y;
						tmp.at<Vec3b>(p_tmp.y, p_tmp.x) = Vec3b(0, 255, 255);
						for (int i = -2; i <= 2; i++)
						{
							for (int j = -2; j <= 2; j++)
							{
								if (cross.at<Vec3b>(p_tmp.y, p_tmp.x)[0]
									< cross.at<Vec3b>(y + i, x + j)[0])
								{
									p_tmp.y = y + i;
									p_tmp.x = x + j;
								}
							}
						}
						tmp.at<Vec3b>(p_tmp.y, p_tmp.x) = Vec3b(0, 0, 255);
						x += 48;
						cross_points[ci].push_back(p_tmp);
					}
				}
			}
		for (auto cp : cross_points)
		{
			cout << cp.size() << endl;
		}
	}
	
	// find center points
	{
		IMG_WIDTH = img.cols;
		vector<Point> line_head;
		{
			Point p_start(-1, -1);
			const int head_low_limit = 5;
			int head_initail_interval_x = 2, head_initail_interval_y = 2;
			// find first point
			for (int y = 0; y < mask.rows; y++)
			{
				for (int x = 0; x < mask.cols; x++)
				{
					if (mask.at<byte>(y, x) > 0)
					{
						p_start.x = x;
						p_start.y = y;
						break;
					}
				}
				if (p_start.x != -1)
					break;
			}
			//cout << p_start << endl;
			update(p_start);
			for (int x = 0; x < mask.cols; x++)
			{
				if (mask.at<byte>(p_start.y, x) > 0
					&& img.at<byte>(p_start.y, x) > head_low_limit)
				{
					p_start.x = x;
					break;
				}
			}
			update_line_head_max(p_start);
			update(p_start);
			//cout << p_start << endl;
			//line_head.push_back(p_start);

			// find head point
			int i1; // , tmp = img.cols / 2;
			for (int i = img.rows - 1; i > 0; i--)
			{
				if (mask.at<byte>(i, p_start.x) > 0
					&& (img.at<byte>(i, p_start.x) > head_low_limit
						|| img.at<byte>(i, p_start.x - head_initail_interval_x) > head_low_limit))
				{
					i1 = i;
					break;
				}
			}
			cout << "line range: " << p_start.y << "->" << i1 << endl;
			int xx = p_start.x;
			while (p_start.y < i1)
			{
				if (line_head.size() == 0 || p_start.y - line_head[line_head.size() - 1].y > 0)
					line_head.push_back(p_start);
				//if (green)
					//p_start.x = xx;
				//else
				//{
					p_start.x += head_initail_interval_x;
					head_initail_interval_x *= -1;
				//}
				p_start.y += head_initail_interval_y;
				//update_line_head_max(p_start);
				update(p_start);
			}
		}
		{
			stack<Point> stack_of_points;
			Point cur, last;
			const int low_limit = 3;
			int  interval_x = 4; //18; //10;
			//if (green)
				//interval_x = 4;
			centers_vec.resize(line_head.size());
			data.resize(line_head.size());
			for (int hi = 0; hi < line_head.size(); hi++)
				//int hi = 0;
			{
				VectorXd center_region(9), center_error_region(9);
				center_error_region << -1, -1, -1, -1, -1, -1, -1, -1, -1;
				// left
				cur = line_head[hi];
				stack_of_points.push(cur);
				while (mask.at<byte>(cur.y, cur.x) > 0
					/*&& mask.at<byte>(cur.y, cur.x - interval_x) > 0*/)
				{
					last = cur;
					cur.x -= interval_x;
					if (mask.at<byte>(cur.y, cur.x) == 0)
						break;
					update(cur);
					//update(cur);
					if (mask.at<byte>(cur.y, cur.x) == 0)
						break;
					stack_of_points.push(cur);
					if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
					{
						cur = last;
						cur.x -= interval_x;
					}
				}
				while (!stack_of_points.empty())
				{
					last = cur;
					cur = stack_of_points.top();
					stack_of_points.pop();
					if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
					{
						cur = last;
						cur.x += interval_x;
						centers_error.push_back(cur);

						centers_vec[hi].push_back(cur);
						data[hi].push_back(center_error_region);
					}
					else
					{
						centers_vec[hi].push_back(cur);
						center_region << img.at<byte>(cur.y - 1, cur.x - 1),
							img.at<byte>(cur.y - 1, cur.x),
							img.at<byte>(cur.y - 1, cur.x + 1),

							img.at<byte>(cur.y, cur.x - 1),
							img.at<byte>(cur.y, cur.x),
							img.at<byte>(cur.y, cur.x + 1),

							img.at<byte>(cur.y + 1, cur.x - 1),
							img.at<byte>(cur.y + 1, cur.x),
							img.at<byte>(cur.y + 1, cur.x + 1);
						data[hi].push_back(center_region);
					}
				}

				// right
				cur = line_head[hi];
				while (1 /*mask.at<byte>(cur.y, cur.x) > 0
					&& mask.at<byte>(cur.y, cur.x + interval_x) > 0*/)
				{
					last = cur;
					cur.x += interval_x;
					if (mask.at<byte>(cur.y, cur.x) == 0) {
						break;
					}
					update(cur);
					//update(cur);
					if (mask.at<byte>(cur.y, cur.x) == 0) {
						break;
					}
					if (img.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
					{
						cur = last;
						cur.x += interval_x;
						centers_error.push_back(cur);

						centers_vec[hi].push_back(cur);
						data[hi].push_back(center_error_region);
					}
					else
					{
						centers_vec[hi].push_back(cur);
						center_region << img.at<byte>(cur.y - 1, cur.x - 1),
							img.at<byte>(cur.y - 1, cur.x),
							img.at<byte>(cur.y - 1, cur.x + 1),

							img.at<byte>(cur.y, cur.x - 1),
							img.at<byte>(cur.y, cur.x),
							img.at<byte>(cur.y, cur.x + 1),

							img.at<byte>(cur.y + 1, cur.x - 1),
							img.at<byte>(cur.y + 1, cur.x),
							img.at<byte>(cur.y + 1, cur.x + 1);
						data[hi].push_back(center_region);
					}
				}
			}
		}

		{
			cout << "head size: " << line_head.size() << endl;
			int sum = 0;
			cout << centers_vec.size() << "*" << centers_vec[510].size() << endl;
			for (auto p : centers_vec)
			{
				sum += p.size();
				for (auto pp : p)
				{
					img_copy.at<Vec3b>(pp.y, pp.x) = Vec3b(0, 255, 0);
				}
			}
			for (auto p : line_head)
				img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 0, 255);
			cout << "total ceneter points size : " << sum << endl;
			cout << "centers_error size : " << centers_error.size() << endl;
			for (auto p : centers_error)
			{
				img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
				//cout << p << endl;
			}
			imwrite(output_select_point_file, img_copy);
		}
		//{
		//	img_copy = Mat(img_copy.size(), CV_8UC3, Scalar(0, 0, 0));
		//	vector<int> region({ -1,0,1 });
		//	for (int vj = 0; vj < centers_vec.size(); vj++)
		//	{
		//		for (int vi = 0; vi < centers_vec[vj].size(); vi++)
		//		{
		//			for (int i = 0; i < region.size(); i++)
		//			{
		//				for (int j = 0; j < region.size(); j++)
		//				{
		//					int x = centers_vec[vj][vi].x + region[j],
		//						y = centers_vec[vj][vi].y + region[i];
		//					byte t = img.at<byte>(y, x);
		//					img_copy.at<Vec3b>(y, x) = Vec3b(t, t, t);
		//				}
		//			}
		//		}
		//	}
		//	imwrite(outfile_set_3x3_region_name, img_copy);
		//}
	}
	pp.endAndPrint();
	imwrite("./output/mask.png", mask);
	imwrite("./output/find_cross.png", tmp);
	p.endAndPrint();
}

// change from find_OLED_location_with_mask
void process(const Mat& rgb, const Mat& mask,
	const char *outfile_selected_point_name,
	const char *outfile_set_3x3_region_name, 
	vector<pair<Point, Point>> cross_points, // start from index 9
	vector<vector<Point>>& centers_vec,
	vector<vector<VectorXd>>& data,
	vector<Point>& centers_error, int start_cross, bool is_green)
{
	Performance p;
	Mat img_copy = rgb.clone(),
		image = rgb.clone();
	cvtColor(image, image, CV_BGR2GRAY);
	line(img_copy, cross_points[0].first, cross_points[6].first, Scalar(0, 255, 0), 1);
	imwrite(outfile_set_3x3_region_name, img_copy);
	cout << "[find all points]in process..." << endl;

	Performance pp;
	//vector<Point> line_head;
	//{
	//	Point p_start(-1, -1);
	//	const int head_low_limit = 5;
	//	int head_initail_interval_x = 4, head_initail_interval_y = 4;
	//	// find first point
	//	for (int y = 0; y < mask.rows; y++)
	//	{
	//		for (int x = 0; x < mask.cols; x++)
	//		{
	//			if (mask.at<byte>(y, x) > 0)
	//			{
	//				p_start.x = x;
	//				p_start.y = y;
	//				break;
	//			}
	//		}
	//		if (p_start.x != -1)
	//			break;
	//	}
	//	update(image, p_start);
	//	for (int x = 0; x < mask.cols; x++)
	//	{
	//		if (mask.at<byte>(p_start.y, x) > 0
	//			&& image.at<byte>(p_start.y, x) > head_low_limit)
	//		{
	//			p_start.x = x;
	//			break;
	//		}
	//	}
	//	update_line_head_max(image, p_start);
	//	update(image, p_start);
	//	// find head point
	//	int i1;
	//	for (int i = image.rows - 1; i > 0; i--)
	//	{
	//		if (mask.at<byte>(i, p_start.x) > 0
	//			&& (image.at<byte>(i, p_start.x) > head_low_limit
	//				|| image.at<byte>(i, p_start.x - head_initail_interval_x) > head_low_limit))
	//		{
	//			i1 = i;
	//			break;
	//		}
	//	}
	//	cout << "line range: " << p_start.y << "->" << i1 << endl;
	//	int xx = p_start.x;
	//	while (p_start.y < i1)
	//	{
	//		if (line_head.size() == 0 || p_start.y - line_head[line_head.size() - 1].y > 3)
	//			line_head.push_back(p_start);
	//		if (is_green)
	//			p_start.x = xx;
	//		else
	//		{
	//			p_start.x += head_initail_interval_x;
	//			head_initail_interval_x *= -1;
	//		}
	//		p_start.y += head_initail_interval_y;
	//		update_line_head_max(image, p_start);
	//		update(image, p_start);
	//	}
	//}
	//const int start_cross = 9;
	{
		stack<Point> stack_of_points;
		Point cur, last;
		const int low_limit = 3;
		int  interval_x = 9;
		if (is_green)
			interval_x = 4;
		//centers_vec.resize(line_head.size());
		centers_vec.resize(cross_points.size() - start_cross);
		data.resize(cross_points.size() - start_cross);
		for (int hi = 0; hi < centers_vec.size(); hi++)
		{
			VectorXd center_region(9), center_error_region(9);
			center_error_region << -1, -1, -1, -1, -1, -1, -1, -1, -1;
			// left
			cur = cross_points[hi+ start_cross].first;
			stack_of_points.push(cur);
			while (mask.at<byte>(cur.y, cur.x) > 0)
			{
				last = cur;
				cur.x -= interval_x;
				if (mask.at<byte>(cur.y, cur.x) == 0)
					break;
				update(image, cur);
				update(image, cur);
				if (mask.at<byte>(cur.y, cur.x) == 0)
					break;
				stack_of_points.push(cur);
				if (image.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x -= interval_x;
				}
			}
			while (!stack_of_points.empty())
			{
				last = cur;
				cur = stack_of_points.top();
				stack_of_points.pop();
				if (image.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x += interval_x;
					centers_error.push_back(cur);
					centers_vec[hi].push_back(cur);
					data[hi].push_back(center_error_region);
				}
				else
				{
					centers_vec[hi].push_back(cur);
					center_region << image.at<byte>(cur.y - 1, cur.x - 1),
						image.at<byte>(cur.y - 1, cur.x),
						image.at<byte>(cur.y - 1, cur.x + 1),
						image.at<byte>(cur.y, cur.x - 1),
						image.at<byte>(cur.y, cur.x),
						image.at<byte>(cur.y, cur.x + 1),
						image.at<byte>(cur.y + 1, cur.x - 1),
						image.at<byte>(cur.y + 1, cur.x),
						image.at<byte>(cur.y + 1, cur.x + 1);
					data[hi].push_back(center_region);
				}
			}
			// right
			cur = cross_points[hi + start_cross].first;
			while (1)
			{
				last = cur;
				cur.x += interval_x;
				if (mask.at<byte>(cur.y, cur.x) == 0) {
					break;
				}
				update(image, cur);
				update(image, cur);
				if (mask.at<byte>(cur.y, cur.x) == 0) {
					break;
				}
				if (image.at<byte>(cur.y, cur.x) < low_limit || abs(cur.y - last.y) > 1)
				{
					cur = last;
					cur.x += interval_x;
					centers_error.push_back(cur);
					centers_vec[hi].push_back(cur);
					data[hi].push_back(center_error_region);
				}
				else
				{
					centers_vec[hi].push_back(cur);
					center_region << image.at<byte>(cur.y - 1, cur.x - 1),
						image.at<byte>(cur.y - 1, cur.x),
						image.at<byte>(cur.y - 1, cur.x + 1),
						image.at<byte>(cur.y, cur.x - 1),
						image.at<byte>(cur.y, cur.x),
						image.at<byte>(cur.y, cur.x + 1),
						image.at<byte>(cur.y + 1, cur.x - 1),
						image.at<byte>(cur.y + 1, cur.x),
						image.at<byte>(cur.y + 1, cur.x + 1);
					data[hi].push_back(center_region);
				}
			}
		}
	}
	pp.endAndPrint();

	{
		cout << "head size: " << cross_points.size() - 9 << endl;
		int sum = 0;
		cout << centers_vec.size() << "*" << centers_vec[510].size() << endl;
		for (auto p : centers_vec)
		{
			sum += p.size();
			for (auto pp : p)
			{
				img_copy.at<Vec3b>(pp.y, pp.x) = Vec3b(0, 255, 0);
			}
		}
		for (auto p : cross_points)
			img_copy.at<Vec3b>(p.first.y, p.first.x) = Vec3b(0, 0, 255);
		cout << "total ceneter points size : " << sum << endl;
		cout << "centers_error size : " << centers_error.size() << endl;
		for (auto p : centers_error)
		{
			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
		}
		imwrite(outfile_selected_point_name, img_copy);
	}
	{
		/*img_copy = Mat(img_copy.size(), CV_8UC3, Scalar(0, 0, 0));
		vector<int> region({ -1,0,1 });
		for (int vj = 0; vj < centers_vec.size(); vj++)
		{
			for (int vi = 0; vi < centers_vec[vj].size(); vi++)
			{
				for (int i = 0; i < region.size(); i++)
				{
					for (int j = 0; j < region.size(); j++)
					{
						int x = centers_vec[vj][vi].x + region[j],
							y = centers_vec[vj][vi].y + region[i];
						byte t = image.at<byte>(y, x);
						img_copy.at<Vec3b>(y, x) = Vec3b(t, t, t);
					}
				}
			}
		}
		imwrite(outfile_set_3x3_region_name, img_copy);*/
	}
	p.endAndPrint();
}

void location_one_column(const Mat& rgb, const Mat& mask,
	vector<pair<Point, Point>>& relation, int region_size, int interval, bool is_green)
{
	Mat img_copy = rgb.clone(), image = rgb.clone();
	cvtColor(image, image, CV_BGR2GRAY);
	//line(img_copy,relation[0].first, relation[6].first, Scalar(0,255,0), 1);
	//imwrite("./output/tmp.png", img_copy);

	double k = (relation[6].first.y - relation[0].first.y) / (relation[6].first.x - relation[0].first.x);
	if (relation[6].first.y == relation[0].first.y)
		k = 0;
	double b = relation[6].first.y - k * relation[6].first.x;
	Point p1(relation[0].first.x, relation[0].first.y), p2;
	//cout << p1 << endl;
	// down
	if (is_green)
	{
		int x = relation[0].second.x, y = relation[0].second.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y++) });
			//img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p1.y += region_size;
			p1.x = (p1.y - b) / k;
		}
		// up
		p1.y = relation[0].first.y - region_size;
		p1.x = (p1.y - b) / k;
		y = relation[0].second.y-1;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y--) });
			//img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p1.y -= region_size;
			p1.x = (p1.y - b) / k;
		}
	}
	else
	{
		int x = relation[0].second.x, y = relation[0].second.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y++) });
			//img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p2.x = p1.x + region_size;
			p2.y = p1.y + region_size;
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y++) });
			//img_copy.at<Vec3b>(p2.y, p2.x) = Vec3b(0, 255, 0);
			p1.y += interval;
			p1.x = (p1.y - b) / k;
		}
		if (mask.at<byte>((*relation.rbegin()).first.y, (*relation.rbegin()).first.x) == 0)
			relation.pop_back();
		// up
		p1.y = relation[0].first.y - interval;
		p1.x = (p1.y - b) / k;
		y = relation[0].second.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			y -= 2;
			update(image, p1);
			relation.push_back({ p1, Point(x, y) });
			//img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p2.x = p1.x + region_size;
			p2.y = p1.y + region_size;
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y + 1) });
			//img_copy.at<Vec3b>(p2.y, p2.x) = Vec3b(0, 255, 0);
			p1.y -= interval;
			p1.x = (p1.y - b) / k;
		}
		p2.x = p1.x + region_size;
		p2.y = p1.y + region_size;
		if (mask.at<byte>(p2.y, p2.x) != 0)
		{
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y - 1) });
		}
	}
	//imwrite("./output/tmp.png", img_copy);
}

void find_OLED_location_with_rgb_combination(
	vector<Mat>& rgb, Mat& mask,
	const char *cross_file,
	const char *output_selected_points_prefix,
	int pentile_width,
	vector<vector<vector<Point>>>& centers_vec,
	vector<vector<vector<VectorXd>>>& data,
	vector<vector<Point>>& centers_error,
	vector<vector<pair<Point,Point>>> cross_points)
{

	IMG_WIDTH = rgb[0].cols;
	IMG_HEIGHT = rgb[0].rows;
	
	{
		const int cross_low_limit = 30, select_range = 500, select_interval = 9;
		Mat cross = imread(cross_file), tmp = cross.clone();
		if (cross.data == NULL)
		{
			cout << "cross file read error: " << cross_file << endl;
			return;
		}
		// find corner 0,1,2; 3,4,5; 6,7,8
		vector<Point> corners;
		{
			get_box_corner_points(mask, corners, 0);
			tmp = cross.clone();
			for (auto p : corners)
			{
				circle(tmp, p, 2, Scalar(0, 255, 0), 2);
			}
		}

		// find cross
		//vector<vector<Point>> cross_points(corners.size());
		set<pair<int,int>> cross_points_set;
		// this can read from cross_point_xy.txt
		vector<int> need_y{ 130, 1216, 2306 }, need_x({ 30,60,90, 530,560,590, 1030,1060,1090 });
		vector<Point> locate_xy;
		/*in the capture file
		  g 0		1		2
		  b 3		4		5
		  r 6		7		8

		  g 9		10		11
		  b 12		13		14
		  r 15		16		17
		  
		  g 18		19		20
		  b 21		22		23
		  r 24		25		26
		*/
		for (int j = 0; j < need_x.size(); j++)
		{
			for (int i = need_y.size() - 1; i >= 0; i--)
			{
				locate_xy.push_back(Point(pentile_width - need_y[i] - (j % 3 == 0 ? 2 : 1), need_x[j]));
			}
		}
		{
			int x0 = corners[0].x, y0 = corners[0].y;
			Point p_tmp;
			for (int ci = 0; ci < corners.size(); ci++)
				for (int y = corners[ci].y - select_range; y <= corners[ci].y + select_range; y++)
				{
					for (int x = corners[ci].x - select_range; x <= corners[ci].x + select_range; x++)
					{
						if (cross.at<Vec3b>(y, x)[0] > cross_low_limit
							&& cross.at<Vec3b>(y + select_interval, x)[0] > cross_low_limit
							&& cross.at<Vec3b>(y, x + select_interval)[0] > cross_low_limit
							&& cross.at<Vec3b>(y - select_interval, x)[0] > cross_low_limit
							&& cross.at<Vec3b>(y, x - select_interval)[0] > cross_low_limit
							&& cross.at<Vec3b>(y + select_interval, x + select_interval)[0] < cross_low_limit / 5)
						{
							//circle(tmp, Point(x0,y0), 2, Scalar(0, 255, 0), 2);
							p_tmp.x = x;
							p_tmp.y = y;
							//tmp.at<Vec3b>(p_tmp.y, p_tmp.x) = Vec3b(0, 255, 255);
							// update(p_tmp)
							for (int i = -2; i <= 2; i++)
							{
								for (int j = -2; j <= 2; j++)
								{
									if (cross.at<Vec3b>(p_tmp.y, p_tmp.x)[0]
										< cross.at<Vec3b>(y + i, x + j)[0])
									{
										p_tmp.y = y + i;
										p_tmp.x = x + j;
									}
								}
							}
							//tmp.at<Vec3b>(p_tmp.y, p_tmp.x) = Vec3b(0, 0, 255);
							x += 48;
							//cross_points[ci].push_back(p_tmp);
							cross_points_set.insert({ p_tmp.y, p_tmp.x });
						}
					}
				}
			cout << cross_points_set.size() << endl;
			int index = -1;
			cross_points.resize(3);
			for (auto cp : cross_points_set)
			{
				//cout << cp.size() << endl;
				Vec3b color_tmp(0, 0, 0);
				//index = (index + 1) % 9;
				cross_points[index % 9 / 3].push_back({ Point(cp.second, cp.first),locate_xy[++index] });
				color_tmp[index % 9 / 3] = 255;
				tmp.at<Vec3b>(cp.first, cp.second) = color_tmp;
				cout << cp.first << "," << cp.second <<"->"<< locate_xy[index] << endl;
			}
			char out_cross[MAX_PATH];
			sprintf(out_cross, "%s/find_cross.png", output_selected_points_prefix);
			imwrite(out_cross, tmp);
		}
		for (int i = 0; i < rgb.size(); i++)
		//int i = 1;
		{
			location_one_column(rgb[i], mask, cross_points[i], 4, 9, i == GREEN);
		}
	}

	for (int i = 0; i < rgb.size(); i++)
	//int i = 1;
	{
		cout << (i == RED ? "r" : (i == BLUE ? "b" : "g")) << endl;
		char selected[MAX_PATH], region3x3[MAX_PATH];
		sprintf(selected, "%s/select_points_%s.png",
			output_selected_points_prefix, i == RED ? "r" : (i == BLUE ? "b" : "g"));
		sprintf(region3x3, "%s/region_3x3_%s.png",
			output_selected_points_prefix, i == RED ? "r" : (i == BLUE ? "b" : "g"));
		cout <<"[output] :"<< selected << endl;
		cout << output_selected_points_prefix << endl;
		process(rgb[i], mask, selected, region3x3, cross_points[i], 
			centers_vec[i], data[i], centers_error[i], 9, i == GREEN);
		// first 9 points are cross points
	}
}
