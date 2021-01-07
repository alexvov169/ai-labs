#include <iostream>
#include <vector>
#include "hopfield_NN.h"
#include <iomanip>
#include <time.h>
#include <fstream>

using namespace std;

#pragma pack(push)
#pragma pack(1)
struct bmp_header
{
	uint16_t bm;
	uint32_t file_size;
	uint32_t reserved;
	uint32_t data_offset;
	uint32_t header_size;
	uint32_t width;
	uint32_t height;
	uint16_t planes;
	uint16_t bitcount;
	uint32_t compression;
	uint32_t image_size;
	uint32_t x_pix_per_meter;
	uint32_t y_pix_per_meter;
	uint32_t colors_used;
	uint32_t colors_important;
};
#pragma pack(pop)

void output_vect(image img)
{
	for (size_t i = 0; i < img.heigth; i++)
	{
		for (size_t j = 0; j < img.width; j++)
		{
			if(img.data[i * img.width + j] == 1)
				cout << setw(3) << 1 << ", ";
			else
				cout << setw(3) << 0 << ", ";
		}
		cout << endl;
	}

	cout << endl;
}

image read_file(string filename)
{
	ifstream file(filename, ios::binary);
	if (!file.is_open())
		exit(-1);
	bmp_header header;
	file.read((char*)&header, sizeof(bmp_header));

	file.seekg(header.data_offset);
	int mx3 = (3 * header.width + 3) & (-4);

	image img;
	img.heigth = header.height;
	img.width = header.width;
	img.data.resize(header.height * header.width, 0);
	int bytes = 0;
	for (size_t i = 0; i < header.height; i++)
	{
		file.seekg((double)header.data_offset + (double)i * mx3, file.beg);
		for (size_t j = 0; j < header.width; j++)
		{
			char color[3];
			file.read(color, 3);
			bytes += 3;

			if (color[0] == '\0' && color[1] == '\0' && color[2] == '\0')
			{
				img.data[(header.height - 1 - i) * header.width + j] = -1;
			}
			else
			{
				img.data[(header.height - 1 - i) * header.width + j] = 1;
			}
		}
	}
	cout << "2\n";

	file.close();
	return img;
}

vector<string> filenames = {
	"A.bmp",
	"H.bmp",
	"M.bmp",
	"T.bmp",
	"U.bmp",
	"V.bmp",
	"X.bmp"
};


int main()
{
	srand(time(NULL));
	size_t errors = 0;
	cout << "To run statistics use 0 errors." << endl;
	cout << "Errors: ";
	cin >> errors;
	vector<image> images;
	for (auto& f : filenames)
	{
		images.push_back(read_file("src/" + f));
	}
	hopfield_NN n1(images[0].data.size());

	n1.calculace_weights(images);

	if (errors > images[0].data.size())
	{
		cout << "Max error size: " << images[0].data.size();
		return 0;
	}

	cout << "Errors: " << errors << endl;

	if (errors == 0)
	{
		for (size_t err = 0; err < images[0].data.size(); err++)
		{
			int correct = 0;
			int runs = 20;
			for (int r = 0; r < runs; r++)
			{
				vector<image> stat_images = images;
				for (auto& img : stat_images)
				{
					for (size_t i = 0; i < err; i++)
					{
						int randpos = rand() % img.data.size();
						img.data[randpos] = -img.data[randpos];
					}
				}

				vector<image> results;
				for (size_t i = 0; i < stat_images.size(); i++)
				{
					results.push_back(stat_images[i]);
					results[i].data = n1.run(stat_images[i].data);
					bool is_ok = true;
					for (size_t pos = 0; pos < results[i].data.size(); pos++)
					{
						if (results[i].data[pos] != images[i].data[pos])
						{
							is_ok = false;
							break;
						}
					}
					if (is_ok)
						correct += 1;

				}
			}
			cout << "Input pixles errors: " << setw(8) << (double)err / images[0].data.size() * 100.0 << "%, Neural Net errors: " << setw(8) << (100.0 - (double)correct / images.size() * 100.0 / runs) << "%." << endl;
		}
		return 0;
	}


	for (auto& img : images)
	{
		for (size_t i = 0; i < errors; i++)
		{
			int randpos = rand() % img.data.size();
			img.data[randpos] = -img.data[randpos];
		}
	}

	//output_vect(images[0]);
	//
	vector<image> results;
	for (size_t i = 0; i < images.size(); i++)
	{
		cout << "File: " << filenames[i] << endl;
		output_vect(images[i]);

		results.push_back(images[i]);
		results[i].data = n1.run(images[i].data);

		output_vect(results[i]);
		cout << "---------------------" << endl;
	}


	//output_vect(images_test);
	//output_vect(result);

	return 0;
}
