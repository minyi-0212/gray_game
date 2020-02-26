#include "location.h"

using namespace std;
using namespace cv;

Mat img;
vector<int> t_update_small({ -2, -1, 0, 1, 2 });
vector<int> t_update_large({ -2, -1, 0, 1, 2 });
vector<int> t_small({ -2, -1, 0, 1, 2 });
//vector<int> t_large({ -9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9 });
vector<int> t_large({ -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,-8 });
void update(Point& p)
{
	Point tmp = p;
	for (int i = 0; i < t_update_small.size(); i++)
	{
		for (int j = 0; j < t_update_large.size(); j++)
		{
			/*if (img.at<Vec3b>(p.y, p.x)[0] < img.at<Vec3b>(tmp.y+ t_update[j], tmp.x+ t_update[i])[0])
			{
				p.y = tmp.y + t_update[j];
				p.x = tmp.x + t_update[i];
			}*/
			if (img.at<Vec3b>(p.y, p.x)[0] < 
				img.at<Vec3b>(p.y + t_update_large[j], p.x + t_update_small[i])[0])
			{
				p.y += t_update_large[j];
				p.x += t_update_small[i];
			}
		}
	}

	/*for (int i = 0; i < t_small.size(); i++)
	{
		if (img.at<Vec3b>(p.y, p.x)[0] < img.at<Vec3b>(tmp.y, tmp.x + t_small[i])[0])
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
			//cout << (int)img.at<Vec3b>(corner[i].y + t_small[ty], corner[i].x + t_large[tx])[0] << " ";
			if (img.at<Vec3b>(tmp.y, tmp.x)[0]
				< img.at<Vec3b>(corner[i].y + t_small[ty], corner[i].x + t_large[tx])[0])
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
			//cout << (int)img.at<Vec3b>(corner[i].y + t19[ty], corner[i].x + t5[tx])[0] << " ";
			if (img.at<Vec3b>(tmp.y, tmp.x)[0]
				< img.at<Vec3b>(corner[i+1].y + t_large[ty], corner[i+1].x + t_small[tx])[0])
			{
				tmp.x = corner[i+1].x + t_small[tx];
				tmp.y = corner[i+1].y + t_large[ty];
			}
		}
	}
	corner[i+1] = tmp;

	if ((corner[i].x - corner[i + 1].x)*(corner[i].x - corner[i + 1].x) +
		(corner[i].y - corner[i + 1].y)*(corner[i].y - corner[i + 1].y) < 8)
	{
		if (img.at<Vec3b>(corner[i].y, corner[i].x)[0]
			< img.at<Vec3b>(corner[i+1].y, corner[i+1].x)[0])
		{
			corner[i] = corner[i + 1];
		}
		else
		{
			corner[i+1] = corner[i];
		}
	}
}

void find_OLED_location(vector<Point>& centers)
{
	Performance p;
	img = imread("G32.bmp");
	int i0 = -1, i1, j0 = -1, j1;

	// find corner - not use
	//{
	//	Mat img_copy = img * 80;
	//	//GaussianBlur(img_copy, img_copy, Size(5, 5), 2, 2);
	//	/*for (int i = 0; i < img_copy.rows; i++)
	//	{
	//		for (int j = 0; j < img_copy.cols; j++)
	//		{
	//			if (img_copy.at<Vec3b>(i, j)[0] != 255
	//				|| img_copy.at<Vec3b>(i, j)[1] != 255
	//				|| img_copy.at<Vec3b>(i, j)[2] != 255)
	//			{
	//				img_copy.at<Vec3b>(i, j)[0] = 0;
	//				img_copy.at<Vec3b>(i, j)[1] = 0;
	//				img_copy.at<Vec3b>(i, j)[2] = 0;
	//			}
	//		}
	//	}*/
	//	int tmp = img_copy.cols / 2, i0_old = -1;
	//	for (int i = 0; i < img_copy.rows; i++)
	//	{
	//		if (img_copy.at<Vec3b>(i, tmp) == Vec3b(255, 255, 255)
	//			/*&& img_copy.at<Vec3b>(i, tmp + 1) == Vec3b(255, 255, 255)
	//			&& img_copy.at<Vec3b>(i, tmp - 1) == Vec3b(255, 255, 255)*/)
	//		{
	//			if (i0 == -1) i0 = i;
	//			i1 = i;
	//		}
	//	}
	//	tmp = img_copy.rows / 2;
	//	for (int j = 0; j < img_copy.cols; j++)
	//	{
	//		if (img_copy.at<Vec3b>(tmp, j) == Vec3b(255, 255, 255)
	//			/*&& img_copy.at<Vec3b>(tmp - 1, j) == Vec3b(255, 255, 255)
	//			&& img_copy.at<Vec3b>(tmp + 1, j) == Vec3b(255, 255, 255)*/)
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
			if (img.at<Vec3b>(i, tmp)[0] > low_limit
				|| img.at<Vec3b>(i, tmp - initail_interval)[0] > low_limit)
			{
				i0 = i;
				break;
			}
		}
		for (int i = img.rows - 1; i > 0; i--)
		{
			if (img.at<Vec3b>(i, tmp)[0] > low_limit
				|| img.at<Vec3b>(i, tmp - initail_interval)[0] > low_limit)
			{
				i1 = i;
				break;
			}
		}
		tmp = img.rows / 2;
		for (int j = 0; j < img.cols; j++)
		{
			if (img.at<Vec3b>(tmp, j)[0] > low_limit
				&&img.at<Vec3b>(tmp - initail_interval, j)[0] > low_limit)
			{
				j0 = j;
				break;
			}
		}
		for (int j = img.cols - 1; j > 0; j--)
		{
			if (img.at<Vec3b>(tmp, j)[0] > low_limit
				&&img.at<Vec3b>(tmp - initail_interval, j)[0] > low_limit)
			{
				j1 = j;
				break;
			}
		}
	}

	// computer the count of points in one line
	/*Mat img_copy = img;
	float k = (1011 - 1009)*1.0 / (9603 - 337),
		b = 1009 - 337 * k;
	cout << k << " " << b << endl;
	cout <<  - b / k << endl;
	vector<Point> cnt_points;
	for (int j = 0; j < img_copy.cols; j++)
	{
		if (img_copy.at<Vec3b>(k*j + b, j)[0] > low_limit)
		{
			if (cnt_points.size() == 0
				|| j - cnt_points[cnt_points.size() - 1].x > 5)
			{
				cnt_points.push_back(Point(j, k*j + b));
				img_copy.at<Vec3b>(k*j + b, j) = Vec3b(0, 255, 0);
			}
		}
	}
	cout << "size of point in one line: " << cnt_points.size() << endl;
	imwrite("./output/G32_test.png", img_copy);*/

	// computer the count of points in one bar
	/*Mat img_copy = img;
	float k = (6220 - 1008)*1.0 / (334 - 337),
		b = 1008 - 337 * k;
	cout << k << " " << b << endl;
	cout <<  - b / k << endl;
	vector<Point> cnt_points;
	for (int y = 0; y < img_copy.rows; y++)
	{
		if (img_copy.at<Vec3b>(y, (y-b)/k)[0] > 10)
		{
			if (cnt_points.size() == 0
				|| y - cnt_points[cnt_points.size() - 1].y > 6)
			{
				cnt_points.push_back(Point((y - b) / k, y));
				img_copy.at<Vec3b>(y, (y - b) / k) = Vec3b(0, 255, 0);
			}
		}
	}
	cout << "size of point in one line: " << cnt_points.size() << endl;
	imwrite("./output/G32_test.png", img_copy);*/

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
	}

	// draw corner
	/*Mat img_copy = img;
	line(img_copy, corner[0], corner[2], Scalar(0, 0, 255), 1);
	line(img_copy, corner[1], corner[5], Scalar(0, 0, 255), 1);
	line(img_copy, corner[4], corner[6], Scalar(0, 0, 255), 1);
	line(img_copy, corner[7], corner[3], Scalar(0, 0, 255), 1);
	for (int i = 0; i < corner.size(); i += 2)
	{
		img_copy.at<Vec3b>(corner[i].y, corner[i].x) = Vec3b(255, 0, 255);
		img_copy.at<Vec3b>(corner[i+1].y, corner[i+1].x) = Vec3b(0, 255, 255);
	}
	imwrite("./output/G32_corner.png", img_copy);*/
	pp.endAndPrint();

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
////				img_copy.at<Vec3b>(i, j) = Vec3b(0, 255, 255);
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
//				img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
//				update(p);
//				centers.push_back(p);
//			}
//		}
//		p.endAndPrint();
//
//		cout << "cornor index : " << i0 << " " << j0 << " " << i1 << " " << j1 << endl;
//		cout << "internal: " << addi << " " << addj << endl;
//		cout << "point size: " << endl << 1920 * 1080 << endl << centers.size() << endl;
//		/*img_copy.at<Vec3b>(i0, j0) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i1, j1) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i1, j0) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i0, j1) = Vec3b(0, 255, 0);*/
//		for (auto p : centers)
//		{
//			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 0, 255);
//		}
//
//		imwrite("./output/G32_position.png", img_copy);
//	}

	{
		//cout << "-----------------" << endl;
		//vector<Point> centers, line_head;
		vector<Point> line_head;
		{
			float k1 = (corner[1].y - corner[5].y)*1.0 / (corner[1].x - corner[5].x),
				b1 = corner[1].y - corner[1].x * k1,
				b2 = corner[0].y - corner[0].x * k1,
				add = (corner[5].y - corner[1].y) / 720.0;
			//cout << "add in y:"<< add << endl;
			for (float y1 = corner[1].y, y2 = corner[0].y;
				y1 <= corner[5].y, y2 <= corner[6].y; y1 += add, y2 += add)
			{
				line_head.push_back(Point((y2 - b2) / k1, y2));
				line_head.push_back(Point((y1 - b1) / k1, y1));
			}
		}
		{
			float k = (corner[0].y - corner[2].y)*1.0 / (corner[0].x - corner[2].x),
				add = (corner[2].x - corner[0].x) / 960.0;
			//cout << "add in x:" << add << endl;
#pragma parallel omp for
			for (auto p : line_head)
			{
				for (int x = p.x; x < corner[2].x; x += add)
				{
					float b = p.y - p.x*k;
					Point p(x, x*k + b);
					update(p);
					update(p);
					centers.push_back(p);
				}
			}
		}
		p.endAndPrint();

		float real_value = 255.0 / 32;
		Mat img_copy = img * real_value;
		for (auto p : line_head)
		{
			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
		}
		for (auto p : centers)
		{
			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 0, 255);
		}
		//imwrite("./output/G32_position.png", img_copy);
	}
}

//void find_OLED_location()
//{
//	Performance p;
//	img = imread("G32.bmp");
//	int i0 = -1, i1, j0 = -1, j1;
//
//	// find corner - not use
//	//{
//	//	Mat img_copy = img * 80;
//	//	//GaussianBlur(img_copy, img_copy, Size(5, 5), 2, 2);
//	//	/*for (int i = 0; i < img_copy.rows; i++)
//	//	{
//	//		for (int j = 0; j < img_copy.cols; j++)
//	//		{
//	//			if (img_copy.at<Vec3b>(i, j)[0] != 255
//	//				|| img_copy.at<Vec3b>(i, j)[1] != 255
//	//				|| img_copy.at<Vec3b>(i, j)[2] != 255)
//	//			{
//	//				img_copy.at<Vec3b>(i, j)[0] = 0;
//	//				img_copy.at<Vec3b>(i, j)[1] = 0;
//	//				img_copy.at<Vec3b>(i, j)[2] = 0;
//	//			}
//	//		}
//	//	}*/
//	//	int tmp = img_copy.cols / 2, i0_old = -1;
//	//	for (int i = 0; i < img_copy.rows; i++)
//	//	{
//	//		if (img_copy.at<Vec3b>(i, tmp) == Vec3b(255, 255, 255)
//	//			/*&& img_copy.at<Vec3b>(i, tmp + 1) == Vec3b(255, 255, 255)
//	//			&& img_copy.at<Vec3b>(i, tmp - 1) == Vec3b(255, 255, 255)*/)
//	//		{
//	//			if (i0 == -1) i0 = i;
//	//			i1 = i;
//	//		}
//	//	}
//	//	tmp = img_copy.rows / 2;
//	//	for (int j = 0; j < img_copy.cols; j++)
//	//	{
//	//		if (img_copy.at<Vec3b>(tmp, j) == Vec3b(255, 255, 255)
//	//			/*&& img_copy.at<Vec3b>(tmp - 1, j) == Vec3b(255, 255, 255)
//	//			&& img_copy.at<Vec3b>(tmp + 1, j) == Vec3b(255, 255, 255)*/)
//	//		{
//	//			if (j0 == -1) j0 = j;
//	//			j1 = j;
//	//		}
//	//	}
//	//	cout << i0 << " " << i1 << " " << j0 << " " << j1 << endl;
//	//	line(img_copy, Point(j0, i0), Point(j1, i1), Scalar(0, 0, 255), 1);
//	//	line(img_copy, Point(330, 1009), Point(9610, 6231), Scalar(255,0,0), 1);
//	//	imwrite("./output/G32_cpoy.png", img_copy);
//	//	//waitKey(0);
//	//}
//
//	Performance pp;
//	// find corner
//	{
//		const int low_limit = 5, initail_interval = 6;
//		int tmp = img.cols / 2;
//		for (int i = 0; i < img.rows; i++)
//		{
//			if (img.at<Vec3b>(i, tmp)[0] > low_limit
//				|| img.at<Vec3b>(i, tmp - initail_interval)[0] > low_limit)
//			{
//				i0 = i;
//				break;
//			}
//		}
//		for (int i = img.rows - 1; i > 0; i--)
//		{
//			if (img.at<Vec3b>(i, tmp)[0] > low_limit
//				|| img.at<Vec3b>(i, tmp - initail_interval)[0] > low_limit)
//			{
//				i1 = i;
//				break;
//			}
//		}
//		tmp = img.rows / 2;
//		for (int j = 0; j < img.cols; j++)
//		{
//			if (img.at<Vec3b>(tmp, j)[0] > low_limit
//				&&img.at<Vec3b>(tmp - initail_interval, j)[0] > low_limit)
//			{
//				j0 = j;
//				break;
//			}
//		}
//		for (int j = img.cols - 1; j > 0; j--)
//		{
//			if (img.at<Vec3b>(tmp, j)[0] > low_limit
//				&&img.at<Vec3b>(tmp - initail_interval, j)[0] > low_limit)
//			{
//				j1 = j;
//				break;
//			}
//		}
//	}
//	pp.endAndPrint();
//
//	{
//		cout << "-----------------" << endl;
//		vector<Point> centers;
//		float real_value = 255.0 / 32,
//			addi = 7,
//			addj = 9;
//		Mat img_copy = img * real_value;
//#pragma parallel omp for
//		for (float i = i0; i <= i1; i += addi)
//			//int i = i0;
//		{
//			//cout << (int)((i-i0)/addi + 1) % 2 * addj / 2 << endl;
//			//for (float j = j0 + (int)((i - i0) / addi + 1) % 2 * addj; j <= j1; j += addj)
//			//for (float j = corner[1].x; j <= corner[1].x; j += addj)
//			for (float j = j0; j <= j1; j += addj)
//			{
//				Point p(j, i);
//				img_copy.at<Vec3b>(i, j) = Vec3b(0, 255, 255);
//				update(p);
//				//update(p);
//				centers.push_back(p);
//			}
//		}
//		p.endAndPrint();
//
//		cout << "cornor index : " << i0 << " " << j0 << " " << i1 << " " << j1 << endl;
//		cout << "internal: " << addi << " " << addj << endl;
//		cout << "point size: " << centers.size() << endl;
//		/*img_copy.at<Vec3b>(i0, j0) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i1, j1) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i1, j0) = Vec3b(0, 255, 0);
//		img_copy.at<Vec3b>(i0, j1) = Vec3b(0, 255, 0);*/
//		for (auto p : centers)
//		{
//			img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 0, 255);
//		}
//
//		imwrite("./output/G32_position.png", img_copy);
//	}
//}