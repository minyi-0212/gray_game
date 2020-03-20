#include "location.h"
#include <fstream>
#include <stack>

using namespace std;
using namespace cv;
using namespace Eigen;

Mat img;
int IMG_WIDTH;
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
	//	imwrite("./output_pentile/G32_cpoy.png", img_copy);
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
		ofstream out("./output_pentile/one_row.csv");
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
		imwrite("./output_pentile/G32_one_row.png", img_copy);
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
		ofstream out("./output_pentile/one_column.csv");
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
		imwrite("./output_pentile/G32_one_column.png", img_copy * 8);
	}*/

	// draw corner
	/*cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
	Mat img_copy = img;
	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255, 0, 0), 1);
	imwrite("./output_pentile/G32_cpoy.png", img_copy);*/

	// draw corner
	/*cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
	Mat img_copy = img;
	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255, 0, 0), 1);
	imwrite("./output_pentile/G32_cpoy.png", img_copy);*/

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
	imwrite("./output_pentile/G32_corner.png", img_copy);*/
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
//		imwrite("./output_pentile/G32_position.png", img_copy);
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
//		//imwrite("./output_pentile/G32_position.png", img_copy);
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
		imwrite("./output_pentile/G32_pick_point.png", img_copy * 8);

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
		imwrite("./output_pentile/G32_after_set.png", img_copy);
		//cvtColor(img_copy, img_copy, CV_BGR2GRAY);
		//imwrite("./output_pentile/G32_after_set.bmp", img_copy);
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