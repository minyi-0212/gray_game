#include "location.h"
#include <fstream>
using namespace std;
using namespace cv;
using namespace Eigen;

Mat img;
int IMG_WIDTH;
vector<int> t_update_small({ -2, -1, 0, 1, 2 });
vector<int> t_update_large({ -2, -1, 0, 1, 2 });
vector<int> t_small({ -2, -1, 0, 1, 2 });
//vector<int> t_large({ -9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9 });
vector<int> t_large({ -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,-8 });
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

//void find_OLED_location(unordered_map<int, VectorXd>& centers)
void find_OLED_location(unordered_map<int, VectorXd>& centers)
{
	Performance p;
	img = imread("G32.bmp");
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
	/*  0 ���� 2
	   1  ����  3

	   5  ����  7
		4 ���� 6
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
		const int low_limit = 5, initail_interval = 5;
		vector<Point> line_head;
		{
			float k = (corner[1].y - corner[5].y)*1.0 / (corner[1].x - corner[5].x),
				b1 = corner[1].y - corner[1].x * k,
				b2 = corner[0].y - corner[0].x * k;
			for (float y = corner[1].y; y <= corner[5].y; y++)
			{
				if (img.at<byte>(y, (y - b1) / k) > low_limit &&
					(line_head.size() == 0 || y - line_head[line_head.size() - 1].y > initail_interval))
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
				if (img.at<byte>(y, (y - b2) / k) > low_limit &&
					(line_head.size() == tmp || y - line_head[line_head.size() - 1].y > initail_interval))
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
			}

			k = (corner[0].y - corner[2].y)*1.0 / (corner[0].x - corner[2].x);
			VectorXd center_region(9);
			Point last;
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
				centers.insert({ line_head[hi].y*IMG_WIDTH + line_head[hi].x, center_region });
				last = line_head[hi];
				for (int x = line_head[hi].x++; x < corner[3].x; x++)
				{
					float b = line_head[hi].y - line_head[hi].x*k;
					if (img.at<byte>(k*x + b, x) > low_limit
						&& x - last.x > initail_interval)
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
						/*update(p);
						update(p);*/
						center_region << img.at<byte>(p.y - 1, p.x - 1),
							img.at<byte>(p.y - 1, p.x),
							img.at<byte>(p.y - 1, p.x + 1),

							img.at<byte>(p.y, p.x - 1),
							img.at<byte>(p.y, p.x),
							img.at<byte>(p.y, p.x + 1),

							img.at<byte>(p.y + 1, p.x - 1),
							img.at<byte>(p.y + 1, p.x),
							img.at<byte>(p.y + 1, p.x + 1);
						if(img.at<byte>(p.y, p.x) > 12)
							//centers[p.y*IMG_WIDTH+p.x] = center_region;
							centers.insert({ p.y*IMG_WIDTH + p.x, center_region });

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
		
		for (auto p : centers)
		{
			img_copy.at<Vec3b>(p.first / IMG_WIDTH, p.first % IMG_WIDTH) = Vec3b(0, 255, 0);
		}
		cout << "size of head points: " << line_head.size() << endl;
		cout << 960*720*2 << ", size of points: " << centers.size() << endl;
		imwrite("./output/G32_pick_point.png", img_copy * 8);

		img_copy = Mat::zeros(img_copy.size(), CV_8UC1);
		vector<int> region({ -1,0,1 });
		for (auto p : centers)
		{
			for (int i = 0; i < region.size(); i++)
			{
				for (int j = 0; j < region.size(); j++)
				{
					int x = p.first % IMG_WIDTH + region[j],
						y = p.first / IMG_WIDTH + region[i];
					img_copy.at<byte>(y, x) = img.at<byte>(y, x);
				}
			}
		}
		imwrite("./output/G32_after_set.png", img_copy);
		imwrite("./output/G32_after_set.bmp", img_copy);
	}
}