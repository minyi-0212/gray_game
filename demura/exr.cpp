#include "exr.h"
#include <windows.h>
#include <fstream>
using namespace std;
using namespace Imf;
using namespace Imf_2_2;
bool save_exr_with_rgb(const char fileName[], int width, int height, 
	const Imf::Array<Imf::Rgba> &pixels, const int compression)
{
	// add compression
	Header header(width, height, 1, IMATH_NAMESPACE::V2f(0, 0),  1, INCREASING_Y, 
		(Compression)compression);
	/*
	enum Compression
	{
		NO_COMPRESSION  = 0,	// no compression
		RLE_COMPRESSION = 1,	// run length encoding
		ZIPS_COMPRESSION = 2,	// zlib compression, one scan line at a time
		ZIP_COMPRESSION = 3,	// zlib compression, in blocks of 16 scan lines
		PIZ_COMPRESSION = 4,	// piz-based wavelet compression
		PXR24_COMPRESSION = 5,	// lossy 24-bit float compression
		B44_COMPRESSION = 6,	// lossy 4-by-4 pixel block compression,
						// fixed compression rate
		B44A_COMPRESSION = 7,	// lossy 4-by-4 pixel block compression,
						// flat fields are compressed more
		DWAA_COMPRESSION = 8,       // lossy DCT based compression, in blocks
									// of 32 scanlines. More efficient for partial
									// buffer access.
		DWAB_COMPRESSION = 9,       // lossy DCT based compression, in blocks
									// of 256 scanlines. More efficient space
									// wise and faster to decode full frames
									// than DWAA_COMPRESSION.
		NUM_COMPRESSION_METHODS	// number of different compression methods
	};*/
	//Imf::RgbaOutputFile file(fileName, width, height, Imf::WRITE_RGBA);
	Imf::RgbaOutputFile file(fileName, header, Imf::WRITE_RGBA);
	//cout << file.compression() << endl;
	file.setFrameBuffer(pixels, 1, width);
	try
	{
		file.writePixels(height);
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << std::endl;
		return false;
	}

	return true;
}

bool save_exr_with_float(const char fileName[], int width, int height, 
	const std::vector<float>& pixels, const int compression)
{
	Imf::Array<Imf::Rgba> pixelsdata;
	pixelsdata.resizeErase(width * height);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;

			half r = pixels[j * 3];
			half g = pixels[j * 3 + 1];
			half b = pixels[j * 3 + 2];
			pixelsdata[j] = Imf::Rgba(r, g, b);
		}
	}

	return save_exr_with_rgb(fileName, width, height, pixelsdata, compression);
}

bool save_exr_with_float_gray(const char fileName[], int width, int height, const std::vector<float>& pixels)
{
	Imf::Array<Imf::Rgba> pixelsdata;
	pixelsdata.resizeErase(width * height);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;

			half g = pixels[j];
			pixelsdata[j] = Imf::Rgba(g, g, g);
		}
	}

	return save_exr_with_rgb(fileName, width, height, pixelsdata);
}

void readPixel(const char *fileName, Imf::Rgba& data, int x, int y) {
	char s[MAX_PATH];
	strcpy(s, fileName);
	Imf::RgbaInputFile file(fileName, 1);
	Imath::Box2i dw = file.dataWindow();
	int width = dw.max.x - dw.min.x + 1;
	int height = dw.max.y - dw.min.y + 1;
	Imf::Array2D<Imf::Rgba> pixels(1, width);
	dw.min.y += y;
	file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
	file.readPixels(dw.min.y);
	data = pixels[0][x];
}

bool read_exr(const char fileName[], std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);

	Imf::Array<Imf::Rgba> pixelsdata;
	int width, height;
	{
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		int dw, dh, dx, dy;
		width = dw = dataWindow.max.x - dataWindow.min.x + 1;
		height = dh = dataWindow.max.y - dataWindow.min.y + 1;
		dx = dataWindow.min.x;
		dy = dataWindow.min.y;

		pixelsdata.resizeErase(dw * dh);
		in.setFrameBuffer(pixelsdata - dx - dy * dw, 1, dw);
		try
		{
			in.readPixels(dataWindow.min.y, dataWindow.max.y);
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}
	}

	pixels.resize(height*width*3);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;
			const Imf::Rgba &rp = pixelsdata[j];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);
		}
	}

	return true;
}

bool read_exr(const char fileName[], int& width, int& height, std::vector<float> &pixels, bool bgr)
{
	char s[MAX_PATH];
	strcpy(s, fileName);
	fstream _file;
	_file.open(s, ios::in);
	if (!_file)
	{
		cout << s << " is not exist" << endl;
		return -1;
	}
	_file.close();


	Imf::Array<Imf::Rgba> pixelsdata;
	{
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		int dw, dh, dx, dy;
		width = dw = dataWindow.max.x - dataWindow.min.x + 1;
		height = dh = dataWindow.max.y - dataWindow.min.y + 1;
		dx = dataWindow.min.x;
		dy = dataWindow.min.y;

		pixelsdata.resizeErase(dw * dh);
		in.setFrameBuffer(pixelsdata - dx - dy * dw, 1, dw);
		try
		{
			in.readPixels(dataWindow.min.y, dataWindow.max.y);
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}
	}

	pixels.resize(height*width * 3);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;
			const Imf::Rgba &rp = pixelsdata[j];

			if (bgr)
			{
				pixels[3 * j] = float(rp.b);
				pixels[3 * j + 1] = float(rp.g);
				pixels[3 * j + 2] = float(rp.r);
			}
			else
			{
				pixels[3 * j] = float(rp.r);
				pixels[3 * j + 1] = float(rp.g);
				pixels[3 * j + 2] = float(rp.b);
			}
		}
	}

	return true;
}

bool read_exr_in_row(const char fileName[], const int& start_row, const int& read_number_of_lines,
	int& width, int& height, std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);

	Imf::Array<Imf::Rgba> pixelsdata;
	int dw, dh;
	{
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		width = dw = dataWindow.max.x - dataWindow.min.x + 1;
		height = dh = dataWindow.max.y - dataWindow.min.y + 1;
		/*dx = dataWindow.min.x;
		dy = dataWindow.min.y;*/

		dataWindow.min.y += start_row;
		pixelsdata.resizeErase(dw * read_number_of_lines);
		in.setFrameBuffer(pixelsdata - dataWindow.min.x - dataWindow.min.y * dw, 1, dw);
		try
		{
			in.readPixels(dataWindow.min.y, min(dataWindow.min.y + read_number_of_lines - 1, dh));
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}
	}

	pixels.resize(read_number_of_lines*dw * 3);
	for (int y = 0; y < read_number_of_lines; y++)
	{
		for (int x = 0; x < dw; x++)
		{
			int j = y * dw + x;
			const Imf::Rgba &rp = pixelsdata[j];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);
		}
	}

	return true;
}

bool read_exr_in_row(const char fileName[], const int& start_row, const int& read_number_of_lines, std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);

	Imf::Array<Imf::Rgba> pixelsdata;
	int dw, dh;
	{
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		dw = dataWindow.max.x - dataWindow.min.x + 1;
		dh = dataWindow.max.y - dataWindow.min.y + 1;
		/*dx = dataWindow.min.x;
		dy = dataWindow.min.y;*/

		dataWindow.min.y += start_row;
		pixelsdata.resizeErase(dw * read_number_of_lines);
		in.setFrameBuffer(pixelsdata - dataWindow.min.x - dataWindow.min.y * dw, 1, dw);
		try
		{
			in.readPixels(dataWindow.min.y, min(dataWindow.min.y + read_number_of_lines - 1, dh));
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}
	}

	pixels.resize(read_number_of_lines * dw * 3);
	for (int y = 0; y < read_number_of_lines; y++)
	{
		for (int x = 0; x < dw; x++)
		{
			/*int j = y * dw + x;
			const Imf::Rgba &rp = pixelsdata[j];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);*/
			pixels[3 * (y * dw + x)] = float(pixelsdata[y * dw + x].r);
			pixels[3 * (y * dw + x) + 1] = float(pixelsdata[y * dw + x].g);
			pixels[3 * (y * dw + x) + 2] = float(pixelsdata[y * dw + x].b);
		}
	}

	return true;
}

bool read_exrs_in_row(const std::vector<std::string>& files_name, 
	const int& start_row, const int& read_number_of_lines, std::vector<std::vector<float>> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, files_name[0].c_str());
	int dw, dh;
	{
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		dw = dataWindow.max.x - dataWindow.min.x + 1;
		dh = dataWindow.max.y - dataWindow.min.y + 1;
	}

	vector<Imf::Array<Imf::Rgba>> pixelsdata(files_name.size());
#pragma omp parallel for
	for (int i = 0; i < files_name.size(); i++)
	{
		//cout << i << endl;
		char s[MAX_PATH];
		strcpy(s, files_name[i].c_str());
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();

		dataWindow.min.y += start_row;
		pixelsdata[i].resizeErase(dw * read_number_of_lines);
		in.setFrameBuffer(pixelsdata[i] - dataWindow.min.x - dataWindow.min.y * dw, 1, dw);
		/*try
		{
			in.readPixels(dataWindow.min.y, min(dataWindow.min.y + read_number_of_lines - 1, dh));
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}*/
		in.readPixels(dataWindow.min.y, min(dataWindow.min.y + read_number_of_lines - 1, dh));
	}
	
	cout << "to float ... " << endl;
	pixels.resize(files_name.size(), vector<float>(read_number_of_lines * dw * 3, -1));
#pragma omp parallel for
	for (int i = 0; i < files_name.size(); i++)
	{
		for (int y = 0; y < read_number_of_lines; y++)
		{
			for (int x = 0; x < dw; x++)
			{
				pixels[i][3 * (y * dw + x)] = float(pixelsdata[i][y * dw + x].r);
				pixels[i][3 * (y * dw + x) + 1] = float(pixelsdata[i][y * dw + x].g);
				pixels[i][3 * (y * dw + x) + 2] = float(pixelsdata[i][y * dw + x].b);
			}
		}
	}

	return true;
}

bool read_exr_in_block(const char fileName[], 
	const int& x_1, const int& x_2,
	const int& y_1, const int& y_2, std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);
	int dw, width = x_2 - x_1 + 1, height = y_2 - y_1 + 1;
	Imf::Array<Imf::Rgba> pixelsdata;
	{
		int dh;
		Imf::RgbaInputFile in(s, 1);
		Imath::Box2i dataWindow = in.dataWindow();
		dw = dataWindow.max.x - dataWindow.min.x + 1;
		dh = dataWindow.max.y - dataWindow.min.y + 1;

		dataWindow.min.y += y_1;
		pixelsdata.resizeErase(dw * (y_2 - y_1 + 1));
		in.setFrameBuffer(pixelsdata - dataWindow.min.x - dataWindow.min.y * dw, 1, dw);
		try
		{
			in.readPixels(dataWindow.min.y, min(dataWindow.min.y + y_2 - y_1, dh));
		}
		catch (const exception &e)
		{
			std::cerr << e.what() << std::endl;
			return false;
		}
	}

	pixels.resize(width * height * 3);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			/*int j = y * dw + x;
			const Imf::Rgba &rp = pixelsdata[j];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);*/
			pixels[3 * (y * width + x)] = float(pixelsdata[y * dw + x_1 + x].r);
			pixels[3 * (y * width + x) + 1] = float(pixelsdata[y * dw + x_1 + x].g);
			pixels[3 * (y * width + x) + 2] = float(pixelsdata[y * dw + x_1 + x].b);
		}
	}

	return true;
}

void write_tiled_exr(const char fileName[], int width, int height, const std::vector<float>& pixels, int tileWidth, int tileHeight)
{
	Imf::Array<Imf::Rgba> pixelsdata;
	pixelsdata.resizeErase(width * height);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;

			half r = pixels[j * 3];
			half g = pixels[j * 3 + 1];
			half b = pixels[j * 3 + 2];
			pixelsdata[j] = Imf::Rgba(r, g, b);
		}
	}

	Imf::TiledRgbaOutputFile out(fileName,
		width, height, // image size
		tileWidth, tileHeight, // tile size
		Imf::ONE_LEVEL, // level mode
		Imf::ROUND_DOWN, // rounding mode
		Imf::WRITE_RGBA); // channels in file // 1
	out.setFrameBuffer(pixelsdata, 1, width); // 2
	//cout << out.numXTiles() - 1 << " " << out.numYTiles() - 1 << endl;
	out.writeTiles(0, out.numXTiles() - 1, 0, out.numYTiles() - 1); // 3
}

void read_tiled_exr(const char fileName[], int &width, int &height, std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);
	Imf::Array2D<Imf::Rgba> pixelsdata;
	Imf::TiledRgbaInputFile in(s, 1);
	Imath::Box2i dw = in.dataWindow();
	width = dw.max.x - dw.min.x + 1;
	height = dw.max.y - dw.min.y + 1;
	int dx = dw.min.x;
	int dy = dw.min.y;
	pixelsdata.resizeErase(height, width);
	in.setFrameBuffer(&pixelsdata[-dy][-dx], 1, width);
	in.readTiles(0, in.numXTiles() - 1, 0, in.numYTiles() - 1);

	pixels.resize(height*width * 3);
	for (int y = 0; y < height; y++)
	{
		int i = y * width;
		for (int x = 0; x < width; x++)
		{
			int j = i + x;
			const Imf::Rgba &rp = pixelsdata[y][x];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);
		}
	}
}

void read_tiled_exr(const char fileName[], int &width, int &height, std::vector<float> &pixels,
	const int tile_x0, const int tile_x1, const int tile_y0, const int tile_y1)
{
	char s[MAX_PATH];
	strcpy(s, fileName);
	Imf::Array2D<Imf::Rgba> pixelsdata;
	Imf::TiledRgbaInputFile in(s, 1);
	const Imf::Header& header = in.header();
	Imf_2_2::TileDescription tile_desc = header.tileDescription();

	Imath::Box2i dw = in.dataWindow();
	width = (tile_x1 - tile_x0 + 1) * tile_desc.xSize;
	height = (tile_y1 - tile_y0 + 1) * tile_desc.ySize;
	pixelsdata.resizeErase(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);
	in.setFrameBuffer(&pixelsdata[-dw.min.y][-dw.min.x], 1, dw.max.x - dw.min.x + 1);
	in.readTiles(tile_x0, tile_x1, tile_y0, tile_y1);

	int offset_x = tile_x0 * tile_desc.ySize, offset_y = tile_y0 * tile_desc.ySize;
	pixels.resize(width * height * 3, 0);
	// need test
	/*pixels.resize(
		(((tile_x1 + 1) * tile_desc.xSize) > dw.max.x - dw.min.x + 1 ? dw.max.x - dw.min.x + 1 - offset_x : width)
		* (((tile_y1 + 1) * tile_desc.ySize) > dw.max.y - dw.min.y + 1 ? dw.max.y - dw.min.y + 1 - offset_y : height) * 3,
		0);*/

	for (int y = 0; y < height && y + offset_y < dw.max.y - dw.min.y + 1; y++)
	{
		int i = y * width;
		for (int x = 0; x < width && x + offset_x < dw.max.x - dw.min.x + 1; x++)
		{
			int j = i + x;
			const Imf::Rgba &rp = pixelsdata[y + offset_y][x + offset_x];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);
		}
	}
}

void read_tiled_exr(const char fileName[], const int tile_x0, const int tile_x1, const int tile_y0, const int tile_y1, std::vector<float> &pixels)
{
	char s[MAX_PATH];
	strcpy(s, fileName);
	Imf::Array2D<Imf::Rgba> pixelsdata;
	Imf::TiledRgbaInputFile in(s, 1);
	const Imf::Header& header = in.header();
	Imf_2_2::TileDescription tile_desc = header.tileDescription();

	Imath::Box2i dw = in.dataWindow();
	/*int width = min((tile_x1 - tile_x0 + 1) * tile_desc.xSize, dw.max.x - dw.min.x + 1);
	int height = min((tile_y1 - tile_y0 + 1) * tile_desc.ySize, dw.max.y - dw.min.y + 1);*/
	int width = (tile_x1 - tile_x0 + 1) * tile_desc.xSize;
	int height = (tile_y1 - tile_y0 + 1) * tile_desc.ySize;
	pixelsdata.resizeErase(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);
	in.setFrameBuffer(&pixelsdata[-dw.min.y][-dw.min.x], 1, dw.max.x - dw.min.x + 1);
	in.readTiles(tile_x0, tile_x1, tile_y0, tile_y1);

	int offset_x = tile_x0 * tile_desc.ySize, offset_y = tile_y0 * tile_desc.ySize;
	pixels.resize(width * height * 3, 0);
	for (int y = 0; (y < height) && (y + offset_y < dw.max.y - dw.min.y + 1); y++)
	{
		int i = y * width;
		for (int x = 0; (x < width) && (x + offset_x < dw.max.x - dw.min.x + 1); x++)
		{
			int j = i + x;
			const Imf::Rgba &rp = pixelsdata[y + offset_y][x + offset_x];

			pixels[3 * j] = float(rp.r);
			pixels[3 * j + 1] = float(rp.g);
			pixels[3 * j + 2] = float(rp.b);
		}
	}
}

bool is_file_exist(const char fileName[])
{
	fstream _file;
	_file.open(fileName, ios::in);
	if (!_file)
	{
		cout << fileName << " is not exist" << endl;
		return false;
	}
	_file.close();
	return true;
}

using namespace Imf;
#include <ImfChannelList.h>
#include <ImfDeepScanLineOutputFile.h>
#include <ImfOutputFile.h>
vector<string> channel_name{"a", "b", "c", "d", "e"};

bool write_multi_channel_exr(const char fileName[], const int channel_size, const float* gPixels,
	const int width, const int height)
{
	//Header header(width, height); // 1
	////for (int i = 0; i < pixels.size(); i++)
	//{
	//	int i = 0;
	//	header.channels().insert(channel_name[i], Channel(HALF)); // 2
	//	cout << i  << endl;
	//}
	//OutputFile file(fileName, header); // 4
	//FrameBuffer frameBuffer; // 5frameBuffer.insert(channel_name[i], // name // 6
	//Slice(PixelType::HALF, // type // 7
	//	(char *)pixels, // base // 8
	//	sizeof(*pixels) * 1, // xStride// 9
	//	sizeof(*pixels) * width); // yStride// 10
	//file.setFrameBuffer(frameBuffer); // 16
	//file.writePixels(height); // 17

	float pixels[100000];
	for (int i = 0; i < 100000; i++)
		pixels[i] = i+0.1;

	Header header(width, height); // 1
	header.channels().insert("Z", Imf::Channel(Imf::FLOAT)); // 3
	OutputFile file("test.exr", header); // 4
	FrameBuffer fb; // 5
	cout << width << " " << height <<" "<< sizeof(float) * 1 << endl;
	auto a = Slice(Imf::FLOAT, (char *)&pixels[0]);
	cout << a.base << endl;
	fb.insert("aaa", // name // 6
		a
	//	Slice(Imf::FLOAT, // type // 7
	//	(char *)&pixels[0], // base // 8
		//	sizeof(float) * 1, // xStride// 9
		//	sizeof(*pixels) * width)
	); // yStride// 10
	cout << width*height << endl;
	file.setFrameBuffer(fb); // 16
	file.writePixels(height); // 17
}

#include <ImfDeepFrameBuffer.h>
//const string DEEPSCANLINE("deepscanline");
#include <ImfPartType.h>
void writeDeepScanlineFile()
{
	try {
		/*const char filename[],
			Imath::Box2i displayWindow,
			Imath::Box2i dataWindow,
			Imf_2_2::Array2D< float* > dataZ,
			Imf_2_2::Array2D< half* > dataO,
			Imf_2_2::Array2D< unsigned int > sampleCount*/
		const char filename[] = "test.exr";
		int width = 300, height = 100;
		Imath::Box2i displayWindow({ 0,0 }, { width - 1, height - 1 });
		Imath::Box2i dataWindow({ 0,0 }, { width - 1, height - 1 });
		Imf_2_2::Array2D<float*> dataZ(width, height);
		Imf_2_2::Array2D<half*> dataO(width, height);
		Imf_2_2::Array2D<unsigned int> sampleCount(width, height);

		//int height = dataWindow.max.y - dataWindow.min.y + 1;
		//int width = dataWindow.max.x - dataWindow.min.x + 1;
		//cout << dataWindow.max.y - dataWindow.min.y + 1 <<endl;
		//cout << dataWindow.max.x - dataWindow.min.x + 1 <<endl;
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				sampleCount[i][j] = i * 2;
			}
		}

		Header header(displayWindow, dataWindow);
		header.channels().insert("Z", Channel(Imf_2_2::FLOAT));
		header.channels().insert("O", Channel(HALF));
		//header.setType(DEEPSCANLINE);
		//Imf_2_2::DeepScanLineOutputFile file(filename, header);
		//Imf_2_2::DeepFrameBuffer frameBuffer;
		//frameBuffer.insertSampleCountSlice(Slice(Imf_2_2::UINT,
		//	(char *)(&sampleCount[0][0]
		//		- dataWindow.min.x
		//		- dataWindow.min.y * width),
		//	sizeof(unsigned int) * 1, // xStride
		//	sizeof(unsigned int) * width)); // yStride
		//frameBuffer.insert("dataZ",
		//	DeepSlice(Imf_2_2::FLOAT,
		//	(char *)(&dataZ[0][0]
		//		- dataWindow.min.x
		//		- dataWindow.min.y * width),
		//		sizeof(float *) * 1, // xStride for pointer array
		//		sizeof(float *) * width, // yStride for pointer array
		//		sizeof(float) * 1)); // stride for Z data sample
		//frameBuffer.insert("dataO",
		//	DeepSlice(HALF,
		//	(char *)(&dataO[0][0]
		//		- dataWindow.min.x
		//		- dataWindow.min.y * width),
		//		sizeof(half *) * 1, // xStride for pointer array
		//		sizeof(half *) * width, // yStride for pointer array
		//		sizeof(half) * 1)); // stride for O data sample
		//file.setFrameBuffer(frameBuffer);
		//for (int i = 0; i < height; i++)
		//{
		//	for (int j = 0; j < width; j++)
		//	{
		//		//sampleCount[i][j] = getPixelSampleCount(j, i);
		//		dataZ[i][j] = new float[sampleCount[i][j]];
		//		dataO[i][j] = new half[sampleCount[i][j]];
		//		// Generate data for dataZ and dataO.
		//	}
		//	file.writePixels(1);
		//}
		//for (int i = 0; i < height; i++)
		//	for (int j = 0; j < width; j++)
		//	{
		//		delete[] dataZ[i][j];
		//		delete[] dataO[i][j];
		//	}
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}
}