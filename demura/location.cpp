#include "location.h"
#include <fstream>
#include <stack>
#include <algorithm>

using namespace std;
using namespace cv;
using namespace Eigen;

int IMG_WIDTH, IMG_HEIGHT;
const vector<int> t_update_small({ -1, 0, 1 });
const vector<int> t_update_large({ -1, 0, 1 });
void update(const Mat& m, Point& p)
{
	Point tmp(p);
	for (int i = 0; i < t_update_small.size(); i++)
	{
		for (int j = 0; j < t_update_large.size(); j++)
		{
			if (m.at<byte>(p.y, p.x)
				<= m.at<byte>(tmp.y + t_update_large[j], tmp.x + t_update_small[i]))
			{

				p.y = tmp.y + t_update_large[j];
				p.x = tmp.x + t_update_small[i];
			}
		}
	}
}

void get_mask_corner_points(const Mat& mask, vector<Point>& corners, int low_limit)
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
	// 0 1
	// 2 3
	corners.resize(4);
	corners[0].x = corners[2].x = xs;
	corners[1].x = corners[3].x = xe;
	corners[0].y = corners[1].y = ys;
	corners[2].y = corners[3].y = ye;

	//// 0 1 2
	//// 3 4 5
	//// 6 7 8
	//corners.resize(9);
	//corners[0].x = corners[3].x = corners[6].x = xs;
	//corners[1].x = corners[4].x = corners[7].x = (xe + xs) / 2;
	//corners[2].x = corners[5].x = corners[8].x = xe;
	//corners[0].y = corners[1].y = corners[2].y = ys;
	//corners[3].y = corners[4].y = corners[5].y = (ye + ys) / 2;
	//corners[6].y = corners[7].y = corners[8].y = ye;
}

void location_one_column(const Mat& rgb, const Mat& mask,
	vector<LED_info>& relation, int region_size, double interval, bool is_green,
	const int s = 0, const int e = 2)
{
	Mat img_copy = rgb.clone(), image = rgb.clone();
	cvtColor(image, image, CV_BGR2GRAY);
	//int s = 1, e = 3;
	Point line_from(relation[s].pixel), line_to(relation[e].pixel),
		locate_from(relation[s].locate), locate_to(relation[e].locate);
	double k = (line_from.y == line_to.y)? 0 : (line_from.y - line_to.y) / (line_from.x - line_to.x),
		b = line_to.y - k * line_to.x;
	relation.clear();

	Point p1, p2;
	int x = locate_from.x, y = locate_from.y;
	if (is_green)
	{
		// up
		p1.y = line_from.y - region_size;
		p1.x = (p1.y - b) / k;
		y = locate_from.y - 1;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y--) , LINE_ONE });
			img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p1.y -= region_size;
			p1.x = (p1.y - b) / k;
		}
		reverse(relation.begin(), relation.end());
		// down
		p1.y = line_from.y;
		p1.x = (p1.y - b) / k;
		y = locate_from.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y++) , LINE_ONE });
			img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p1.y += region_size;
			p1.x = (p1.y - b) / k;
		}
	}
	else
	{
		// up
		p1.y = line_from.y - interval;
		p1.x = (p1.y - b) / k;
		y = locate_from.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			y -= 2;
			update(image, p1);
			img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p2.x = p1.x + region_size;
			p2.y = p1.y + region_size;
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y + 1), LINE_ONE });
			relation.push_back({ p1, Point(x, y), LINE_ONE });
			img_copy.at<Vec3b>(p2.y, p2.x) = Vec3b(0, 255, 0);
			p1.y -= interval;
			p1.x = (p1.y - b) / k;
		}
		p2.x = p1.x + region_size;
		p2.y = p1.y + region_size;
		if (mask.at<byte>(p2.y, p2.x) != 0)
		{
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y - 1), LINE_ONE });
		}
		reverse(relation.begin(), relation.end());
		// down
		p1.y = line_from.y;
		p1.x = (p1.y - b) / k;
		y = locate_from.y;
		while (mask.at<byte>(p1.y, p1.x) != 0)
		{
			update(image, p1);
			relation.push_back({ p1, Point(x, y++), LINE_ONE });
			img_copy.at<Vec3b>(p1.y, p1.x) = Vec3b(0, 255, 255);
			p2.x = p1.x + region_size;
			p2.y = p1.y + region_size;
			update(image, p2);
			relation.push_back({ p2, Point(x + 1, y++) , LINE_ONE });
			img_copy.at<Vec3b>(p2.y, p2.x) = Vec3b(0, 255, 0);
			p1.y += interval;
			p1.x = (p1.y - b) / k;
		}
		if (mask.at<byte>((*relation.rbegin()).pixel.y, (*relation.rbegin()).pixel.x) == 0)
			relation.pop_back();
	}
	/*cout << endl << "each line:" << endl;
	for (auto p : relation)
		cout <<p.pixel<<"->"<< p.locate << endl;*/
	//if(is_green)
		imwrite("./output/tmp.png", img_copy);
}

vector<int> locate_interval({ 2,1,2 });
const int OLED_pick_region = 8;
bool is_error(const Mat& image, const Point& p)
{
	return image.at<byte>(p.y, p.x) < image.at<byte>(p.y - 1, p.x - 1)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y - 1, p.x)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y - 1, p.x + 1)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y, p.x - 1)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y, p.x + 1)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y + 1, p.x - 1)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y + 1, p.x)
		|| image.at<byte>(p.y, p.x) < image.at<byte>(p.y + 1, p.x + 1);
}

void vote(vector<pair<Point, int>>& v, const int x, const int y)
{
	int i = 0;
	for (; i < v.size(); i++)
	{
		if (x == v[i].first.x && y == v[i].first.y)
		{
			v[i].second++;
			break;
		}
	}
	if (i == v.size())
		v.push_back({Point(x,y), 1});
}

void find_point_in_region(const Mat& image, const Mat& mask,
	const double interval_x, const double interval_y, const double locate_interval,
	const int index, vector<Point>& locate, vector<Point>& start_p,
	vector<vector<LED_info>>& centers_vec, set<pair<int, int>>& points_set,
	const int low_limit, const int lor = 1, bool is_green = true) // lor 1:right,-1:left
{
	//cout << lor * interval_x;
	Point cur, last, start;
	start = start_p[0];
	vector<pair<Point, int>> vote_result;
	//cout << start << " ";
	for (int y = 0; y < OLED_pick_region && index + y < centers_vec.size(); y++)
	{
		//cout <<" --"<< y << " ";
		cur.y = start.y + y * interval_y;
		for (int x = 0; x < OLED_pick_region; x++)
		{
			cur.x = start.x + lor * x * interval_x + (is_green ? 0 : (y & 1)*interval_x / 2);
			//cout << cur << " ";
			if (mask.at<byte>(cur.y, cur.x) == 0) {
				centers_vec[index + y].push_back({ cur, locate[y], INVALID });
				locate[y].x += lor * locate_interval;
				continue;
			}
			//cout << cur << "->";
			last = cur;
			update(image, cur);
			update(image, cur);
			//cout << cur << " ";
			if (mask.at<byte>(cur.y, cur.x) == 0
				|| points_set.count({ cur.x, cur.y })
				|| is_error(image, cur)
				|| image.at<byte>(cur.y, cur.x) < low_limit) {
				cur = last;
				//cout << last << "->" << cur << " is invalid" << endl;
				centers_vec[index + y].push_back({ cur, locate[y], INVALID });
			}
			else
			{
				if (x == 0 && y == 0) 
					start = cur;
				//if (cur.x == 9913 && cur.y == 3123)
					//cout << index << endl;
				centers_vec[index + y].push_back({ cur, locate[y], /*(x == 0) ? FIRST :*/ VALID });
				if (is_green)
					vote(vote_result,
						cur.x + lor * (OLED_pick_region - x)*interval_x,
						cur.y - y * interval_y);
				else
					vote(vote_result,
						cur.x + lor * (OLED_pick_region - x)*interval_x - y % 2 * interval_x / 2,
						cur.y - y * interval_y);
			}
			points_set.insert({ cur.x, cur.y });
			locate[y].x += lor * locate_interval;
		}
	}
	
	if (vote_result.size())
	{
		sort(vote_result.begin(), vote_result.end(), [](const pair<Point, int>& a, const pair<Point, int>& b) {
			return a.second > b.second;
		});
		start_p[0] = vote_result[0].first;
		for (int y = 1; y < OLED_pick_region && index + y < centers_vec.size(); y++)
		{
			start_p[y] = start_p[0];
			start_p[y].y += y * interval_y;
		}
	}
	else
	{
		for (int y = 0; y < OLED_pick_region; y++)
		{
			start_p[y] = { -1, -1 };
		}
	}
}

void location_all_points(const Mat& rgb, const Mat& mask,
	const char *outfile_selected_point_name,
	const char *outfile_set_3x3_region_name,
	const vector<LED_info>& cross_points,
	vector<vector<LED_info>>& centers_vec,
	RGB select_rgb)
{
	cout << "[find all points] in process total line " << cross_points.size() << "..." << endl;
	Performance p;
	Mat img_copy = rgb.clone(), image = rgb.clone();
	cvtColor(image, image, CV_BGR2GRAY);

	Performance pp;
	{
		Point  cur, last;
		set<pair<int, int>> points_set;
		vector<Point> left(OLED_pick_region), locate(OLED_pick_region);
		const double low_limit = 1, interval_x = select_rgb == GREEN ? 4.5 : 9,
			interval_y = 4.5;
		centers_vec.resize(cross_points.size());
		bool flag;
		for (int hi = 0; hi < centers_vec.size(); hi += OLED_pick_region)
		//int hi = 252;
		{
			//cout <<"hi: "<< hi << endl;
			// right
			for (int li = 0; li < OLED_pick_region && hi + li < cross_points.size(); li++)
			{
				left[li] = cross_points[hi + li].pixel; // line [hi, hi+OLED_pick_region)
				locate[li] = cross_points[hi + li].locate;
			}
			while (1)
			{
				flag = false;
				for (int li = 0; li < OLED_pick_region&& hi + li < cross_points.size(); li++)
				{
					//cout << hi + li <<":"<< left[li].x << ", ";
					//if (flag || mask.at<byte>(left[li].y, left[li].x) != 0)
					if (flag || left[li].x != -1)
						flag = true;
				}
				//cout << endl;
				if (!flag)
					break;
				find_point_in_region(image, mask, interval_x, interval_y, locate_interval[select_rgb],
					hi, locate, left, centers_vec, points_set, low_limit, 1, select_rgb==GREEN);
			}
			//cout << "left" << endl;
			// left
			for (int li = 0; li < OLED_pick_region && hi + li < cross_points.size(); li++)
			{
				left[li] = cross_points[hi + li].pixel; // line [hi, hi+OLED_pick_region)
				left[li].x -= interval_x;
				locate[li] = cross_points[hi + li].locate;
				locate[li].x -= locate_interval[select_rgb];
			}
			while (1)
			{
				flag = false;
				for (int li = 0; li < OLED_pick_region && hi + li < cross_points.size(); li++)
				{
					//if (flag || mask.at<byte>(left[li].y, left[li].x) != 0)
					if (flag || left[li].x != -1)
						flag = true;
				}
				if (!flag)
					break;
				find_point_in_region(image, mask, interval_x, interval_y, locate_interval[select_rgb],
					hi, locate, left, centers_vec, points_set, low_limit, -1, select_rgb == GREEN);
			}
		}
	}
	pp.endAndPrint();

	{
		cout << "head size: " << cross_points.size() << endl;
		int sum = 0;
		cout << "rough total size: " << centers_vec.size() << "*" << centers_vec[510].size() << endl;
		for (auto p : centers_vec)
		{
			sum += p.size();
			for (auto pp : p)
			{
				if (pp.state == INVALID)
					img_copy.at<Vec3b>(pp.pixel.y, pp.pixel.x) = Vec3b(0, 255, 0);
				else if (pp.state == VALID)
					img_copy.at<Vec3b>(pp.pixel.y, pp.pixel.x) = Vec3b(0, 0, 255);
				else if (pp.state == FIRST)
					img_copy.at<Vec3b>(pp.pixel.y, pp.pixel.x) = Vec3b(255, 0, 0);
			}
		}
		//for (auto p : cross_points)
			//img_copy.at<Vec3b>(p.pixel.y, p.pixel.x) = Vec3b(0, 0, 255);
		cout << "total ceneter points size : " << sum << endl;
		//cout << "centers_error size : " << centers_error.size() << endl;
		//for (auto p : centers_error)
		//{
			//img_copy.at<Vec3b>(p.y, p.x) = Vec3b(0, 255, 255);
		//}
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

int get_pixel(const Mat& img, int y, int x)
{
	Point tmp(x, y);
	//cout << tmp << "->";
	update(img, tmp);
	//cout << tmp << endl;
	return img.at<byte>(tmp.y, tmp.x);
}

void find_cross(const Mat& img, const Mat& mask, 
	const char *output_prefix, const int pentile_width,
	const RGB select_rgb, vector<LED_info>& cross_points)
{
	/*const double upscale = 0.6, low_scale = 0.8;
	const int cross_low_limit = 150, select_range = 600;*/
	const double upscale = 0.6, low_scale = 0.6;
	const int cross_low_limit = 90, select_range = 600;
	const double select_interval = 9;
	Mat cross = img.clone(), tmp = cross.clone();
	cvtColor(cross, cross, COLOR_BGR2GRAY);
	// find corner 0,1; 2,3 / 0,1,2; 3,4,5; 6,7,8 
	vector<Point> corners;
	get_mask_corner_points(mask, corners, 0);
	for (auto p : corners)
	{
		circle(tmp, p, 2, Scalar(0, 255, 0), 2);
	}
	// find cross
	vector<Point> cross_points_vec(corners.size(), Point(-1, -1));
	{
		int x0 = corners[0].x, y0 = corners[0].y, cross_value;
		Point p_tmp;
		for (int ci = 0; ci < corners.size(); ci++)
			for (int y = corners[ci].y - select_range; y <= corners[ci].y + select_range; y++)
			{
				for (int x = corners[ci].x - select_range; x <= corners[ci].x + select_range; x++)
				{
					cross_value = cross.at<byte>(y, x);
					p_tmp.x = x;
					p_tmp.y = y;
					/*if (p_tmp.x == 10671 && p_tmp.y == 6623)
						cout << is_error(cross, p_tmp) << " "
						<< (cross.at<byte>(y, x) > cross_low_limit) << " "
						<< (get_pixel(cross, y + select_interval, x) > cross_value* upscale) << " "
						<< get_pixel(cross, y + select_interval, x)<<" "<< cross_value* upscale << " "
						<< (get_pixel(cross, y, x + select_interval) > cross_value* upscale) << " "
						<< (get_pixel(cross, y - select_interval, x) > cross_value* upscale) << " "
						<< (get_pixel(cross, y, x - select_interval) > cross_value* upscale) << " "
						<< (get_pixel(cross, y + 2 * select_interval, x) > cross_value* upscale) << " "
						<< (get_pixel(cross, y, x + 2 * select_interval) > cross_value* upscale) << " "
						<< (get_pixel(cross, y - 2 * select_interval, x) > cross_value* upscale) << " "
						<< (get_pixel(cross, y, x - 2 * select_interval) > cross_value* upscale) << " "
						<< (get_pixel(cross, y + select_interval, x + select_interval) < cross_value*low_scale) << " "
						<< (get_pixel(cross, y + select_interval, x - select_interval) < cross_value*low_scale) << " "
						<< (get_pixel(cross, y - select_interval, x + select_interval) < cross_value*low_scale) << " "
						<< (get_pixel(cross, y - select_interval, x - select_interval) < cross_value*low_scale) << " "
						<< endl;*/
					if (!is_error(cross, p_tmp) 
						&& cross.at<byte>(y, x) > cross_low_limit
						&& ((get_pixel(cross, y + select_interval, x) > cross_value* upscale
							&& get_pixel(cross, y, x + select_interval) > cross_value* upscale
							&& get_pixel(cross, y - select_interval, x) > cross_value* upscale
							&& get_pixel(cross, y, x - select_interval) > cross_value* upscale)
							|| (get_pixel(cross, y + 2 * select_interval, x) > cross_value* upscale
								&& get_pixel(cross, y, x + 2 * select_interval) > cross_value* upscale
								&& get_pixel(cross, y - 2 * select_interval, x) > cross_value* upscale
								&& get_pixel(cross, y, x - 2 * select_interval) > cross_value* upscale))
						&& (get_pixel(cross, y + select_interval, x + select_interval) < cross_value*low_scale
							&& get_pixel(cross, y + select_interval, x - select_interval) < cross_value*low_scale
							&& get_pixel(cross, y - select_interval, x + select_interval) < cross_value*low_scale
							&& get_pixel(cross, y - select_interval, x - select_interval) < cross_value*low_scale))
					{
						x += 48;
						update(cross, p_tmp);
						//cout <<ci << " " << p_tmp << endl;
						cross_points_vec[ci] = cross_points_vec[ci].x == -1 ? p_tmp : 
							(cross.at<byte>(p_tmp.y, p_tmp.x) > 
								cross.at<byte>(cross_points_vec[ci].y, cross_points_vec[ci].x) ?
								p_tmp: cross_points_vec[ci]);
					}
				}
			}

		// the pentile x&y, this can read from cross_point_xy.txt
		///vector<int> need_y{ 130, 1216, 2306 }, need_x({ 30,60,90, 530,560,590, 1030,1060,1090 });
		vector<Point> cross_locate_xy;
		vector<int> need_y{ 130, 2306 }, need_x({ 60 ,1060 });
		for (int j = 0; j < need_x.size(); j++)
		{
			for (int i = need_y.size() - 1; i >= 0; i--)
			{
				cross_locate_xy.push_back(Point(pentile_width - need_y[i] - (select_rgb == 0 ? 2 : 1), need_x[j]));
			}
		}

		int index = -1;
		Vec3b color_tmp(0, 0, 0);
		color_tmp[select_rgb] = 255;
		for (auto cp : cross_points_vec)
		{
			cross_points.push_back({ cp, cross_locate_xy[++index], CROSS });
			cout << cp.y << "," << cp.x << "->" << cross_locate_xy[index] << endl;
			if(cp.x != -1)
				//tmp.at<Vec3b>(cp.y, cp.x) = color_tmp;
				circle(tmp, cp, 2, Scalar(0, 255, 0), 2);
		}
		char out_cross[MAX_PATH];
		sprintf(out_cross, "%s/find_cross_%s.png", output_prefix,
			select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		imwrite(out_cross, tmp);
	}
}

int find_OLED_location_with_rgb_combination(
	const std::vector<cv::Mat>& rgb_image, const cv::Mat& mask,
	const char *output_prefix, const int pentile_width,
	std::vector<std::vector<std::vector<LED_info>>>& centers_vec)
{
	vector<vector<LED_info>> cross_points(3);
	for (int select_rgb = 0; select_rgb < rgb_image.size(); select_rgb++)
	//int select_rgb = 0;
	{
		cout << "--------------------" << endl;
		cout << "locate " << (select_rgb == RED ? "red" : (select_rgb == BLUE ? "blue" : "green")) << endl;
		// find cross
		{
			if (select_rgb == 2)
			{
				cross_points[select_rgb].push_back({ Point(1021, 2151), Point(128, 60), CROSS });
				cross_points[select_rgb].push_back({ Point(-1, -1), Point(2304, 60), CROSS });
				cross_points[select_rgb].push_back({ Point(1010, 6594), Point(129, 1060), CROSS });
				cross_points[select_rgb].push_back({ Point(-1, -1), Point(2304, 1060), CROSS });
			}
			else
			find_cross(rgb_image[select_rgb], mask, output_prefix, pentile_width,
				(RGB)select_rgb, cross_points[select_rgb]);
			/*cross_points[select_rgb].push_back({ Point(1027, 2217), Point(129, 60), CROSS });
			cross_points[select_rgb].push_back({ Point(-1, -1), Point(2305, 60), CROSS });
			cross_points[select_rgb].push_back({ Point(1088, 6627), Point(129, 1060), CROSS });
			cross_points[select_rgb].push_back({ Point(-1, -1), Point(2305, 1060), CROSS });*/
			// find one point in each line(represent each line)
			if (cross_points[select_rgb][0].pixel.x != -1 
				&& cross_points[select_rgb][2].pixel.x != -1)
				location_one_column(rgb_image[select_rgb], mask, cross_points[select_rgb],
					cross_points[select_rgb].size(), 9, select_rgb == GREEN, 0, 2);
			else if (cross_points[select_rgb][1].pixel.x != -1
				&& cross_points[select_rgb][3].pixel.x != -1)
				location_one_column(rgb_image[select_rgb], mask, cross_points[select_rgb],
					cross_points[select_rgb].size(), 9, select_rgb == GREEN, 1, 3);
			else
			{
				cout << "not find cross." << endl;
				return -1;
			}
		}

		char outfile_selected_point_name[MAX_PATH], outfile_set_3x3_region_name[MAX_PATH];
		sprintf(outfile_selected_point_name, "%s/select_points_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		sprintf(outfile_set_3x3_region_name, "%s/region_3x3_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		cout << "[output] :" << endl << outfile_selected_point_name << endl
			<< outfile_set_3x3_region_name << endl;
		location_all_points(rgb_image[select_rgb], mask,
			outfile_selected_point_name, outfile_set_3x3_region_name,
			cross_points[select_rgb], centers_vec[select_rgb], (RGB)select_rgb);

		char csv_name[MAX_PATH];
		sprintf(csv_name, "%s/origin_%s.csv", output_prefix, 
			select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		ofstream out(csv_name);
		for (auto line : centers_vec[select_rgb])
		{
			for (auto p : line)
			{
				out << p.pixel.x << "," << p.pixel.y << ","
					<< p.locate.x << "," << p.locate.y << ",,";
			}
			out << endl;
		}
		out.close();
	}
	return 0;
}

int tmp_valid_find_location(const std::vector<cv::Mat>& rgb_image, const cv::Mat& mask,
	const char *output_prefix, const int pentile_width,
	std::vector<std::vector<std::vector<LED_info>>>& centers_vec)
{
	vector<vector<LED_info>> cross_points(3);
	//for (int select_rgb = 0; select_rgb < rgb_image.size(); select_rgb++)
	int select_rgb = 1;
	{
		// find cross
		find_cross(rgb_image[select_rgb], mask, output_prefix, pentile_width,
			(RGB)select_rgb, cross_points[select_rgb]);
		// find one point in each line(represent each line)
		if (cross_points[select_rgb][0].pixel.x != -1
			&& cross_points[select_rgb][2].pixel.x != -1)
			location_one_column(rgb_image[select_rgb], mask, cross_points[select_rgb],
				cross_points[select_rgb].size(), 9, select_rgb == GREEN, 0, 2);
		else if (cross_points[select_rgb][1].pixel.x != -1
			&& cross_points[select_rgb][3].pixel.x != -1)
			location_one_column(rgb_image[select_rgb], mask, cross_points[select_rgb],
				cross_points[select_rgb].size(), 9, select_rgb == GREEN, 1, 3);
		else
		{
			cout << "not find cross." << endl;
			return -1;
		}

		cout << "locate " << (select_rgb == RED ? "red" : (select_rgb == BLUE ? "blue" : "green")) << endl;
		char outfile_selected_point_name[MAX_PATH], outfile_set_3x3_region_name[MAX_PATH];
		sprintf(outfile_selected_point_name, "%s/select_points_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		sprintf(outfile_set_3x3_region_name, "%s/region_3x3_%s.png",
			output_prefix, select_rgb == RED ? "r" : (select_rgb == BLUE ? "b" : "g"));
		cout << "[output] :" << endl << outfile_selected_point_name << endl
			<< outfile_set_3x3_region_name << endl;
		location_all_points(rgb_image[select_rgb], mask,
			outfile_selected_point_name, outfile_set_3x3_region_name,
			cross_points[select_rgb], centers_vec[select_rgb], (RGB)select_rgb);

		char valid_csv[MAX_PATH];
		sprintf(valid_csv, "%s/valid_result.csv", output_prefix);
		ofstream out(valid_csv);
		for (auto line : centers_vec[select_rgb])
		{
			for (auto p : line)
			{
				out << p.pixel.x << "," << p.pixel.y << ","
					<< p.locate.x << "," << p.locate.y << ",,";
			}
			out << endl;
		}
		out.close();
	}
}
